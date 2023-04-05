/* Copyright (c) 2018–2021 SplineLib

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. */

#include "Sources/ParameterSpaces/non_zero_degree_b_spline_basis_function.hpp"

#include <utility>

#include "Sources/ParameterSpaces/knot_vector.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/numeric_operations.hpp"

namespace splinelib::sources::parameter_spaces {

NonZeroDegreeBSplineBasisFunction::NonZeroDegreeBSplineBasisFunction(
    KnotVector const& knot_vector,
    KnotSpan const& start_of_support,
    Degree degree,
    Tolerance const& tolerance)
    : Base_(knot_vector, start_of_support, std::move(degree), tolerance),
      left_lower_degree_basis_function_{CreateDynamic(knot_vector,
                                                      start_of_support,
                                                      degree_ - Degree{1},
                                                      tolerance)},
      right_lower_degree_basis_function_{
          CreateDynamic(knot_vector,
                        start_of_support + KnotSpan{1},
                        degree_ - Degree{1},
                        tolerance)} {
  Index const start_index{start_of_support.Get()};
  Degree::Type_ const& degree_value = degree_.Get();
  left_denominator_inverse_ = InvertPotentialZero(
      knot_vector[start_index + Index{degree_value}] - start_knot_,
      tolerance);
  right_denominator_inverse_ =
      InvertPotentialZero(end_knot_ - knot_vector[start_index + Index{1}],
                          tolerance);
  left_quotient_derivative_ = (degree_value * left_denominator_inverse_);
  right_quotient_derivative_ = (degree_value * right_denominator_inverse_);
}

NonZeroDegreeBSplineBasisFunction::NonZeroDegreeBSplineBasisFunction(
    KnotVector const& knot_vector,
    KnotSpan const& start_of_support,
    Degree degree,
    UniqueBSplineBasisFunctions& unique_basis_functions,
    Tolerance const& tolerance)
    : Base_(knot_vector, /* generates this_hash_ */
            start_of_support,
            std::move(degree),
            tolerance),
      left_lower_degree_basis_function_{CreateDynamic(knot_vector,
                                                      start_of_support,
                                                      degree_ - Degree{1},
                                                      unique_basis_functions,
                                                      tolerance)},
      right_lower_degree_basis_function_{
          CreateDynamic(knot_vector,
                        start_of_support + KnotSpan{1},
                        degree_ - Degree{1},
                        unique_basis_functions,
                        tolerance)} {
  Index const start_index{start_of_support.Get()};
  Degree::Type_ const& degree_value = degree_.Get();
  left_denominator_inverse_ = InvertPotentialZero(
      knot_vector[start_index + Index{degree_value}] - start_knot_,
      tolerance);
  right_denominator_inverse_ =
      InvertPotentialZero(end_knot_ - knot_vector[start_index + Index{1}],
                          tolerance);
  left_quotient_derivative_ = (degree_value * left_denominator_inverse_);
  right_quotient_derivative_ = (degree_value * right_denominator_inverse_);
}

bool IsEqual(NonZeroDegreeBSplineBasisFunction const& lhs,
             NonZeroDegreeBSplineBasisFunction const& rhs,
             Tolerance const& tolerance) {
  using Base = NonZeroDegreeBSplineBasisFunction::Base_;
  using utilities::numeric_operations::IsEqual;

#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception,
          "splinelib::sources::parameter_spaces::IsEqual::"
          "NonZeroDegreeBSplineBasisFunction");
  }
#endif
  return (IsEqual(lhs.left_denominator_inverse_,
                  rhs.left_denominator_inverse_,
                  tolerance)
          && IsEqual(lhs.right_denominator_inverse_,
                     rhs.right_denominator_inverse_,
                     tolerance)
          && IsEqual(lhs.left_quotient_derivative_,
                     rhs.left_quotient_derivative_,
                     tolerance)
          && IsEqual(lhs.right_quotient_derivative_,
                     rhs.right_quotient_derivative_,
                     tolerance)
          && IsEqual(*lhs.left_lower_degree_basis_function_,
                     *rhs.left_lower_degree_basis_function_,
                     tolerance)
          && IsEqual(*lhs.right_lower_degree_basis_function_,
                     *rhs.right_lower_degree_basis_function_,
                     tolerance)
          && IsEqual(static_cast<Base const&>(lhs),
                     static_cast<Base const&>(rhs),
                     tolerance));
}

bool operator==(NonZeroDegreeBSplineBasisFunction const& lhs,
                NonZeroDegreeBSplineBasisFunction const& rhs) {
  return IsEqual(lhs, rhs);
}

// Recurrence formula due to DeBoor, Cox, and Mansfield (see NURBS book Eq.
// (2.5)).
NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::operator()(
    ParametricCoordinate const& parametric_coordinate,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception,
          "splinelib::sources::parameter_spaces::"
          "NonZeroDegreeBSplineBasisFunction::operator()");
  }
#endif
  return IsInSupport(parametric_coordinate, tolerance) ? (
             ((parametric_coordinate - start_knot_).Get()
              * left_denominator_inverse_
              * (*left_lower_degree_basis_function_)(parametric_coordinate,
                                                     tolerance))
             + ((end_knot_ - parametric_coordinate).Get()
                * right_denominator_inverse_
                * (*right_lower_degree_basis_function_)(parametric_coordinate,
                                                        tolerance)))
                                                       : Type_{};
}

NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::ConsecutiveTopNodeEvaluation(
    ParametricCoordinate const& parametric_coordinate,
    EvaluationLookUp& evaluation_look_up,
    const int& end_support,
    const bool& is_first_support,
    const bool& check_right,
    Tolerance const& tolerance) const {

  using ReturnType = NonZeroDegreeBSplineBasisFunction::Type_;

  const auto degree_minus_one = degree_.Get() - 1;

  // first support node and its children. has very distinctive behavior
  if (is_first_support) {
    // first support just need to compute right nodes
    // look up should be allocated with same number as top level basis
    // degree. that way, each contribution can be saved at the same id
    // as source degree.
    // contribution comes from one degree lower.

    auto& right_contribution = evaluation_look_up[degree_minus_one];
    right_contribution =
        right_lower_degree_basis_function_->ConsecutiveTopNodeEvaluation(
            parametric_coordinate,
            evaluation_look_up,
            end_support,
            true,
            false,
            tolerance);

    return right_contribution * (end_knot_ - parametric_coordinate).Get()
           * right_denominator_inverse_;
  }

  // quite the same procedure for everything but the first.
  // start by initializing return value for current node contribution

  // add left support before it's overwritten by right visit
  // get reference, as we can just write here afterwards
  // at this point saved_contribution has previous tree's right contribution
  auto& saved_contribution = evaluation_look_up[degree_minus_one];
  const ReturnType left_contribution =
      (parametric_coordinate - start_knot_).Get() * left_denominator_inverse_
      * saved_contribution;

  // check if this is out of support
  if (check_right) {
    // current start of the support bigger than the end support?
    if (start_of_support_ > end_support) {
      return ReturnType{0.0};
    }

    // let's try to avoid lower visit
    if (start_of_support_ + 1 > end_support) {
      saved_contribution = 0.0;
      return left_contribution;
    }
  }

  // assign right node contribution of this tree
  saved_contribution =
      right_lower_degree_basis_function_->ConsecutiveTopNodeEvaluation(
          parametric_coordinate,
          evaluation_look_up,
          end_support,
          false,
          true,
          tolerance);

  return left_contribution
         + (saved_contribution * (end_knot_ - parametric_coordinate).Get()
            * right_denominator_inverse_);
}

// Based on recurrence formula due to DeBoor, Cox, and Mansfield (see NURBS book
// Eq. (2.9)).
NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::operator()(
    ParametricCoordinate const& parametric_coordinate,
    Derivative const& derivative,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception,
          "splinelib::sources::parameter_spaces::"
          "NonZeroDegreeBSplineBasisFunction::operator()");
  }
#endif
  if (derivative != Derivative{}) {
    if (IsInSupport(parametric_coordinate, tolerance)) {
      Derivative const lower_derivative = (derivative - Derivative{1});
      return ((left_quotient_derivative_
               * (*left_lower_degree_basis_function_)(parametric_coordinate,
                                                      lower_derivative,
                                                      tolerance))
              - (right_quotient_derivative_
                 * (*right_lower_degree_basis_function_)(parametric_coordinate,
                                                         lower_derivative,
                                                         tolerance)));
    } else {
      return Type_{};
    }
  } else {
    return operator()(parametric_coordinate, tolerance);
  }
}

NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::ConsecutiveTopNodeDerivativeEvaluation(
    ParametricCoordinate const& parametric_coordinate,
    Derivative const& derivative,
    EvaluationLookUp& derivative_look_up,
    EvaluationLookUp& evaluation_look_up,
    const int& end_support,
    const bool& is_first_support,
    const bool& check_right,
    Tolerance const& tolerance) const {

  using ReturnType = NonZeroDegreeBSplineBasisFunction::Type_;

  // zeroth step is to check if the query is actually zeroth derivative
  // because it is then equivalent to evaluation query
  if (derivative == Derivative{}) {
    // always check right
    auto zero_der = ConsecutiveTopNodeEvaluation(parametric_coordinate,
                                                 evaluation_look_up,
                                                 end_support,
                                                 is_first_support,
                                                 true,
                                                 tolerance);
    derivative_look_up[0] = zero_der;
    return zero_der;
  }

  const auto lower_derivative = derivative - Derivative{1};

  // first support node and its children have very distinctive behavior
  if (is_first_support) {
    // first support just need to compute right nodes
    // look up should be allocated with same number as derivative
    // query. that way, each contribution can be saved at the same id
    // as source derivative.
    // contribution comes from one derivative lower.
    auto& right_contribution = derivative_look_up[lower_derivative.Get()];
    right_contribution =
        right_lower_degree_basis_function_
            ->ConsecutiveTopNodeDerivativeEvaluation(parametric_coordinate,
                                                     lower_derivative,
                                                     derivative_look_up,
                                                     evaluation_look_up,
                                                     end_support,
                                                     is_first_support,
                                                     check_right,
                                                     tolerance);

    // returns left minus right, so minus here
    return -(right_quotient_derivative_ * right_contribution);
  }

  // quite the same procedure for everything but the first.
  // start by initializing return value for current node contribution

  // add left support before it's overwritten by right visit.
  // get reference, as we can just write here afterwards
  // at this point saved_contribution has previous tree's right contribution
  auto& saved_contribution = derivative_look_up[lower_derivative.Get()];
  const ReturnType left_contribution =
      left_quotient_derivative_ * saved_contribution;

  // check if this is out of support
  if (check_right) {
    // current start of the support bigger than the end support?
    if (start_of_support_ > end_support) {
      return ReturnType{0.0};
    }

    // let's try to avoid lower visit
    if (start_of_support_ + 1 > end_support) {
      saved_contribution = 0.0;
      return left_contribution;
    }
  }

  // assign right node contribution of this tree
  saved_contribution =
      right_lower_degree_basis_function_
          ->ConsecutiveTopNodeDerivativeEvaluation(parametric_coordinate,
                                                   lower_derivative,
                                                   derivative_look_up,
                                                   evaluation_look_up,
                                                   end_support,
                                                   is_first_support,
                                                   check_right,
                                                   tolerance);

  return left_contribution - (right_quotient_derivative_ * saved_contribution);
}

NonZeroDegreeBSplineBasisFunction::Type_
NonZeroDegreeBSplineBasisFunction::InvertPotentialZero(
    ParametricCoordinate const& potential_zero,
    Tolerance const& tolerance) const {
  return IsEqual(potential_zero, ParametricCoordinate{}, tolerance)
             ? Type{}
             : Type_{1.0 / potential_zero.Get()};
}

} // namespace splinelib::sources::parameter_spaces
