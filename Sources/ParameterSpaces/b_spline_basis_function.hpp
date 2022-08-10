/* Copyright (c) 2018–2021 SplineLib

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. */

#ifndef SOURCES_PARAMETERSPACES_B_SPLINE_BASIS_FUNCTION_HPP_
#define SOURCES_PARAMETERSPACES_B_SPLINE_BASIS_FUNCTION_HPP_

#include "Sources/ParameterSpaces/knot_vector.hpp"
#include "Sources/Utilities/named_type.hpp"
#include "Sources/Utilities/std_container_operations.hpp"
#include "Sources/Utilities/string_operations.hpp" // for String alias.

namespace splinelib::sources::parameter_spaces {

class BSplineBasisFunction;

// proof of concept
// unorderd_map key syntax
// <start_of_support>_<degree>
using UniqueBSplineBasisFunctions =
    UnorderedMap<String, SharedPointer<BSplineBasisFunction>>;

template<int parametric_dimensionality>
using UniqueBSplineBasisFunctionsArray =
    Array<UniqueBSplineBasisFunctions, parametric_dimensionality>;

// it is just a holder to lookup
using UniqueEvaluations = Vector<ParametricCoordinate::Type_>;
// turns out they both can have the same type.
using UniqueDerivatives = UniqueEvaluations;

template<int parametric_dimensionality>
using UniqueEvaluationsArray =
    Array<UniqueEvaluations, parametric_dimensionality>;

template<int parametric_dimensionality>
using UniqueDerivativesArray =
    Array<UniqueDerivatives, parametric_dimensionality>;

using IsTopLevelComputed = Vector<bool>;

template<int parametric_dimensionality>
using IsTopLevelComputedArray = Array<IsTopLevelComputed,
                                      parametric_dimensionality>;

// BSplineBasisFunctions N_{i,p} are non-negative, piecewise polynomial functions of degree p forming a basis of the
// vector space of all piecewise polynomial functions of degree p corresponding to some knot vector.  The B-spline basis
// functions have local support only and form a partition of unity.  Actual implementations for evaluating them (as well
// as their derivatives) on their support are (Non)ZeroDegreeBSplineBasisFunctions.
//
// Example (see NURBS book Exa. 2.1):
//   KnotVector::Knot_ const k0_0{0.0}, k1_0{1.0};
//   KnotVector const kKnotVector{k0_0, k0_0, k0_0, k1_0, k1_0, k1_0};
//   ParametricCoordinate const k0_25{0.25};
//   BSplineBasisFunction const * const basis_function_2_0{CreateDynamic(kKnotVector, KnotSpan{2}, Degree{})},
//                              * const basis_function_0_2{CreateDynamic(kKnotVector, KnotSpan{}, Degree{2})};
//   BSplineBasisFunctions::Type_ const &evaluated_2_0 = (*basis_function_2_0)(k0_25),  // N_{2,0}(0.25) = 1.0.
//                                const &evaluated_0_2 = (*basis_function_0_2)(k0_25);  // N_{0,2}(0.25) = 0.5625.
class BSplineBasisFunction {
 public:
  using Type_ = ParametricCoordinate::Type_;

  static BSplineBasisFunction * CreateDynamic(
      KnotVector const &knot_vector,
      KnotSpan const &start_of_support,
      Degree degree,
      Tolerance const &tolerance = kEpsilon);

  static SharedPointer<BSplineBasisFunction> CreateDynamic(
      KnotVector const &knot_vector,
      KnotSpan const &start_of_support,
      Degree degree,
      UniqueBSplineBasisFunctions &unique_basis_functions,
      Tolerance const &tolerance = kEpsilon);

  virtual ~BSplineBasisFunction() = default;

  // Comparison based on tolerance.
  friend bool IsEqual(BSplineBasisFunction const &lhs, BSplineBasisFunction const &rhs, Tolerance const &tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==(BSplineBasisFunction const &lhs, BSplineBasisFunction const &rhs);
  virtual Type_ operator()(ParametricCoordinate const &parametric_coordinate,
                           Tolerance const &tolerance = kEpsilon) const = 0;
  // tree_info should tell you if requested basis function is ...
  //   <tree_info value> : <info>
  //                  -2 : left branch
  //                  -1 : right branch
  //         0 or bigger : top_node id. used to determine entry id for
  //                       top level evaluations
  virtual Type_ operator()(ParametricCoordinate const &parametric_coordinate,
                           UniqueEvaluations& unique_evaluations,
                           int const &tree_info,
                           Tolerance const &tolerance = kEpsilon) const = 0;
  virtual Type_ operator()(ParametricCoordinate const &parametric_coordinate,
                           Derivative const &derivative,
                           Tolerance const &tolerance = kEpsilon) const = 0;
  virtual Type_ operator()(ParametricCoordinate const &parametric_coordinate,
                           Derivative const &derivative,
                           UniqueDerivatives& unique_derivatives,
                           UniqueEvaluations& unique_evaluations,
                           IsTopLevelComputed& top_level_computed,
                           int const &tree_info,
                           Tolerance const &tolerance = kEpsilon) const = 0;

 protected:
  BSplineBasisFunction() = default;
  BSplineBasisFunction(KnotVector const &knot_vector, KnotSpan const &start_of_support, Degree degree,
                       Tolerance const &tolerance = kEpsilon);
  BSplineBasisFunction(BSplineBasisFunction const &other) = default;
  BSplineBasisFunction(BSplineBasisFunction &&other) noexcept = default;
  BSplineBasisFunction & operator=(BSplineBasisFunction const &rhs) = default;
  BSplineBasisFunction & operator=(BSplineBasisFunction &&rhs) noexcept = default;

  // The last knot is treated in a special way (cf. KnotVector::FindSpan).
  bool IsInSupport(ParametricCoordinate const &parametric_coordinate, Tolerance const &tolerance = kEpsilon) const;

  Degree degree_;
  ParametricCoordinate start_knot_, end_knot_;
  bool end_knot_equals_last_knot_;
};

bool IsEqual(BSplineBasisFunction const &lhs, BSplineBasisFunction const &rhs, Tolerance const &tolerance = kEpsilon);
bool operator==(BSplineBasisFunction const &lhs, BSplineBasisFunction const &rhs);

template<int parametric_dimensionality>
using BSplineBasisFunctions = Array<Vector<SharedPointer<BSplineBasisFunction>>, parametric_dimensionality>;

}  // namespace splinelib::sources::parameter_spaces

#endif  // SOURCES_PARAMETERSPACES_B_SPLINE_BASIS_FUNCTION_HPP_
