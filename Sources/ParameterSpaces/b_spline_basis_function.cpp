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

#include "Sources/ParameterSpaces/b_spline_basis_function.hpp"

#include <utility>

#include "Sources/ParameterSpaces/knot_vector.hpp"
#include "Sources/ParameterSpaces/non_zero_degree_b_spline_basis_function.hpp"
#include "Sources/ParameterSpaces/zero_degree_b_spline_basis_function.hpp"
#include "Sources/Utilities/error_handling.hpp"
#include "Sources/Utilities/numeric_operations.hpp"

namespace splinelib::sources::parameter_spaces {

using std::move;

// This function also creates all basis functions of lower degree required for
// evaluation.
BSplineBasisFunction*
BSplineBasisFunction::CreateDynamic(KnotVector const& knot_vector,
                                    KnotSpan const& start_of_support,
                                    Degree degree,
                                    Tolerance const& tolerance) {
#ifndef NDEBUG
  Message const kName{"splinelib::sources::parameter_spaces::"
                      "BSplineBasisFunction::CreateDynamic"};

  try {
    KnotSpan::ThrowIfNamedIntegerIsOutOfBounds(start_of_support,
                                               (knot_vector.GetSize() - 1)
                                                   - (degree.Get() + 1));
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  if (degree == Degree{}) {
    return new ZeroDegreeBSplineBasisFunction(knot_vector,
                                              start_of_support,
                                              tolerance);
  } else {
    return new NonZeroDegreeBSplineBasisFunction(knot_vector,
                                                 start_of_support,
                                                 move(degree),
                                                 tolerance);
  }
}

// This creates a new basis function iff requested basis function does not
// exist.
// Else shared pointer of existing basis function.
SharedPointer<BSplineBasisFunction> BSplineBasisFunction::CreateDynamic(
    KnotVector const& knot_vector,
    KnotSpan const& start_of_support,
    Degree degree,
    UniqueBSplineBasisFunctions& unique_basis_functions,
    Tolerance const& tolerance) {

  // make a unique identifier str
  const std::string basis_identifier = std::to_string(start_of_support.Get())
                                       + "_" + std::to_string(degree.Get());

  // And does this exist?
  const auto& existing_basis_function =
      unique_basis_functions.find(basis_identifier);

  if (existing_basis_function != unique_basis_functions.end()) {
    // jackpot
    return existing_basis_function->second;

  } else {
    // does not exist. create - store - return
    if (degree == Degree{}) {
      auto& new_basis = unique_basis_functions[basis_identifier];
      new_basis =
          std::make_shared<ZeroDegreeBSplineBasisFunction>(knot_vector,
                                                           start_of_support,
                                                           tolerance);

      return new_basis;
    } else {
      auto& new_basis = unique_basis_functions[basis_identifier];
      new_basis = std::make_shared<NonZeroDegreeBSplineBasisFunction>(
          knot_vector,
          start_of_support,
          move(degree),
          unique_basis_functions,
          tolerance);

      return new_basis;
    }
  }
}

BSplineBasisFunction::BSplineBasisFunction(KnotVector const& knot_vector,
                                           KnotSpan const& start_of_support,
                                           Degree degree,
                                           Tolerance const& tolerance)
    : degree_(move(degree)),
      start_of_support_(start_of_support.Get()),
      start_knot_(knot_vector[Index{start_of_support_}]),
      end_knot_(knot_vector[Index{start_of_support_ + degree_.Get() + 1}]),
      end_knot_equals_last_knot_(
          knot_vector.DoesParametricCoordinateEqualBack(end_knot_, tolerance)) {
}

bool IsEqual(BSplineBasisFunction const& lhs,
             BSplineBasisFunction const& rhs,
             Tolerance const& tolerance) {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(
        exception,
        "splinelib::sources::parameter_spaces::IsEqual::BSplineBasisFunction");
  }
#endif
  return (
      (lhs.degree_ == rhs.degree_)
      && IsEqual(lhs.start_knot_, rhs.start_knot_, tolerance)
      && IsEqual(lhs.end_knot_, rhs.end_knot_, tolerance)
      && (lhs.end_knot_equals_last_knot_ == rhs.end_knot_equals_last_knot_));
}

bool operator==(BSplineBasisFunction const& lhs,
                BSplineBasisFunction const& rhs) {
  return IsEqual(lhs, rhs);
}

bool BSplineBasisFunction::IsInSupport(
    ParametricCoordinate const& parametric_coordinate,
    Tolerance const& tolerance) const {
  return ((IsGreaterOrEqual(parametric_coordinate, start_knot_, tolerance)
           && IsLess(parametric_coordinate, end_knot_, tolerance))
          || (IsEqual(parametric_coordinate, end_knot_, tolerance)
              && end_knot_equals_last_knot_));
}

} // namespace splinelib::sources::parameter_spaces
