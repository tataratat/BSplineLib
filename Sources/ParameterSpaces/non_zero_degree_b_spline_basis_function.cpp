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
NonZeroDegreeBSplineBasisFunction::operator()(
    ParametricCoordinate const& parametric_coordinate,
    UniqueEvaluations& unique_evaluations,
    int const& tree_info,
    Tolerance const& tolerance) const {

  // First, support check - exit if not
  if (!IsInSupport(parametric_coordinate, tolerance)) {
    return Type_{};
  }

  int entry_offset{};
  const auto& id = degree_.Get();
  // evaluation is quite deterministic: it is clear when to lookup
  if (tree_info == -2) {
    // lucky, just take a look.
    // this is valid for all left branches.
    return unique_evaluations[id];

  } else if (tree_info >= 0) {
    // indicates this is top-level nodes
    // check if this has been computed
    const auto& top_level_evaluation = unique_evaluations[id + tree_info];
    // values are initialized with negative value and if computed,
    // it should be non-negative.
    if (top_level_evaluation >= 0.0) {
      return top_level_evaluation;
    }

    entry_offset = tree_info;
  }
  // it is your duty to compute.
  // this means you are either at the top-level node or right branch.
  // explicitly compute left first then right.
  // => fix for windows, as it would compute right first
  auto left_contribution =
      (parametric_coordinate - start_knot_).Get() * left_denominator_inverse_
      * (*left_lower_degree_basis_function_)(parametric_coordinate,
                                             unique_evaluations,
                                             -2,
                                             tolerance);
  auto right_contribution =
      (end_knot_ - parametric_coordinate).Get() * right_denominator_inverse_
      * (*right_lower_degree_basis_function_)(parametric_coordinate,
                                              unique_evaluations,
                                              -1,
                                              tolerance);

  // add
  auto& computed_basis = unique_evaluations[id + entry_offset];
  computed_basis = left_contribution + right_contribution;

  return computed_basis;
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
NonZeroDegreeBSplineBasisFunction::operator()(
    ParametricCoordinate const& parametric_coordinate,
    Derivative const& derivative,
    UniqueDerivatives& unique_derivatives,
    UniqueEvaluations& unique_evaluations,
    IsTopLevelComputed& top_level_computed,
    int const& tree_info,
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

  int entry_offset{}; /* will be zero, unless it is top-level node */
  if (tree_info == -2) {
    return unique_derivatives[derivative.Get()];

  } else if (tree_info >= 0) {
    // indicates this is top-level nodes.
    // check if this has been computed
    //
    // for derivatives, there are no simple rule that tells us if
    // it was evaluated.
    // so, we look at an additional info
    if (top_level_computed[tree_info]) {
      // for 0-th derivative, look at evaluations.
      if (derivative == Derivative{}) {
        throw std::runtime_error(
            String("CRITICAL ISSUE! Please help us by writing an issue at ")
            + String("github.com/tataratat/SplineLib"));
      }
      return unique_derivatives[derivative.Get() + tree_info];
    }

    // save tree_info to place freshly computed value at the right place.
    entry_offset = tree_info;

    // we are here, because it wasn't computed and in following parts,
    // it will be computed.
    // warning: premature action.
    // an alternative would be checking if entry_offset is >= 0 then set true.
    top_level_computed[tree_info] = true;
  }

  // it is your duty to compute
  // top-level-degree derivatives and all right_lower ones.
  if (derivative != Derivative{}) {
    if (!IsInSupport(parametric_coordinate, tolerance)) {
      // tschuess
      return Type_{};
    }

    Derivative const lower_derivative = (derivative - Derivative{1});
    // same idea as evaluation.
    // top-level-derivative and all right ones
    // explicitly compute left first then right.
    // => fix for windows, as it would compute right first
    const auto left_contribution =
        left_quotient_derivative_
        * (*left_lower_degree_basis_function_)(parametric_coordinate,
                                               lower_derivative,
                                               unique_derivatives,
                                               unique_evaluations,
                                               top_level_computed,
                                               -2,
                                               tolerance);
    const auto right_contribution =
        right_quotient_derivative_
        * (*right_lower_degree_basis_function_)(parametric_coordinate,
                                                lower_derivative,
                                                unique_derivatives,
                                                unique_evaluations,
                                                top_level_computed,
                                                -1,
                                                tolerance);
    auto& computed_basis = unique_derivatives[derivative.Get() + entry_offset];
    computed_basis = left_contribution - right_contribution;

    // here would be a good alternative place to set top-level computed flag.
    return computed_basis;

  } else {
    // zeroth derivative evaluation. same as normal evaluation.
    // treat it as normal evaluation <- need top level node info.
    // could use an additional input if this takes too long.
    int top_level_id{-1};
    // compiler told me `bool` will be copied anyways. Hence, no auto&
    for (const auto computed : top_level_computed) {
      if (!computed)
        break;
      ++top_level_id;
    }
    const auto evaluation = operator()(parametric_coordinate,
                                       unique_evaluations,
                                       top_level_id,
                                       tolerance);

    // Here's always 0th derivative.
    unique_derivatives[0] = evaluation;

    return evaluation;
  }
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
