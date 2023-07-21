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

#ifndef SOURCES_PARAMETERSPACES_PARAMETER_SPACE_HPP_
#define SOURCES_PARAMETERSPACES_PARAMETER_SPACE_HPP_

#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>

#include "BSplineLib/ParameterSpaces/knot_vector.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/index.hpp"
#include "BSplineLib/Utilities/math_operations.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"
#include "BSplineLib/Utilities/std_container_operations.hpp"
#include "BSplineLib/Utilities/string_operations.hpp"

namespace bsplinelib::parameter_spaces {

template<int parametric_dimensionality>
class ParameterSpace;

template<int parametric_dimensionality>
bool IsEqual(ParameterSpace<parametric_dimensionality> const& lhs,
             ParameterSpace<parametric_dimensionality> const& rhs,
             Tolerance const& tolerance = kEpsilon);
template<int parametric_dimensionality>
bool operator==(ParameterSpace<parametric_dimensionality> const& lhs,
                ParameterSpace<parametric_dimensionality> const& rhs);

/// recursive combime adapted from bezman
template<std::size_t depth, typename ValueType, std::size_t array_dim>
constexpr void
RecursiveCombine_(const Array<Vector<ValueType>, array_dim>& factors,
                  Vector<ValueType>& result,
                  const ValueType& c_value) {
  static_assert(depth < array_dim,
                "Implementation error, recursion loop to deep!");

  for (std::size_t i{}; i < factors[depth].size(); ++i) {
    if constexpr (depth == 0) {
      result.push_back(c_value * factors[depth][i]);
    } else {
      RecursiveCombine_<static_cast<std::size_t>(depth - 1)>(
          factors,
          result,
          c_value * factors[depth][i]);
    }
  }
}

/// recursive combime adapted from bezman
template<typename ValueType, std::size_t array_dim>
constexpr Vector<ValueType>
RecursiveCombine(const Array<Vector<ValueType>, array_dim>& factors) {
  // Precalculate required entries
  int n_entries{1};
  std::for_each(factors.begin(),
                factors.end(),
                [&](const std::vector<double>& ii) { n_entries *= ii.size(); });

  // Init return type and reserve memory
  Vector<ValueType> result{};
  result.reserve(n_entries);

  // Start computation
  RecursiveCombine_<array_dim - 1>(factors, result, static_cast<ValueType>(1));
  return result;
}

/// recursive combime taken from bezman
template<std::size_t depth,
         typename BasisValueType,
         std::size_t array_dim,
         typename ValueType,
         typename BSplineLibIndex,
         typename CoeffType,
         typename ReturnType>
constexpr void
RecursiveCombine_(const Array<Vector<BasisValueType>, array_dim>& factors,
                  const BSplineLibIndex& index,
                  BSplineLibIndex& index_offset,
                  const CoeffType& coeffs,
                  const ValueType& c_value,
                  ReturnType& result) {
  static_assert(depth < array_dim,
                "Implementation error, recursion loop to deep!");

  for (const auto& factor : factors[depth]) {
    if constexpr (depth == 0) {
      // get basis
      const auto fac = c_value * factor;
      // get coeff
      const auto& coeff =
          coeffs[(index + index_offset.GetIndex()).GetIndex1d().Get()];
      // contribute to each dim
      int j{};
      for (auto& r : result) {
        r += typename ReturnType::value_type{static_cast<ValueType>(coeff[j])
                                             * fac};
        ++j;
      }
      ++index_offset;
    } else {
      RecursiveCombine_<static_cast<std::size_t>(depth - 1)>(factors,
                                                             index,
                                                             index_offset,
                                                             coeffs,
                                                             c_value * factor,
                                                             result);
    }
  }
}

/// recursive combime adapted from bezman
template<typename ValueType,
         std::size_t array_dim,
         typename BSplineLibIndex,
         typename CoeffType,
         typename ReturnType>
constexpr void
RecursiveCombine(const Array<Vector<ValueType>, array_dim>& factors,
                 const BSplineLibIndex& index,
                 BSplineLibIndex& index_offset,
                 const CoeffType& coeffs,
                 ReturnType& result) {

  // Start computation
  RecursiveCombine_<array_dim - 1>(factors,
                                   index,
                                   index_offset,
                                   coeffs,
                                   static_cast<ValueType>(1),
                                   result);
}

template<typename T>
struct TemporaryArray {
  T* data_;
  TemporaryArray(const int n) : data_(new T[n]) {}
  ~TemporaryArray() { delete[] data_; }
  constexpr T& operator[](const int& i) { return data_[i]; }
  constexpr const T& operator[](const int& i) const { return data_[i]; }
};
template<typename T>
struct TemporaryArray2D {
  T* data_;
  int dim_;
  TemporaryArray2D(const int n, const int d) : data_(new T[n * d]), dim_(d) {}
  ~TemporaryArray2D() { delete[] data_; }
  constexpr T* operator[](const int& i) { return &data_[i * dim_]; }
  constexpr const T* operator[](const int& i) const { return &data_[i * dim_]; }
};

// ParameterSpaces provide the B-spline basis functions corresponding to given
// knot vectors and degrees.  Only clamped knot vectors of degree p — i.e., both
// the first knot a and last knot b have multiplicity p+1 (each of them is
// repeated p times) and interior knots have multiplicities that are not greater
// than p — are allowed: U = {u_0 = a, ..., u_p = a, u_{p+1}, ..., u_{m-(p+1)},
// b = u_{m-p}, ..., b = u_m}.
//
// Example (see NURBS book Exa. 2.1):
//   using ParameterSpace2d = ParameterSpace<2>;
//   ParameterSpace2d::Knot_ const k0_0{0.0}, k1_0{1.0};
//   ParameterSpace2d::KnotVectors_::value_type const &kKnotVector =
//   std::make_shared<KnotVector>({k0_0, k0_0, k0_0, k1_0, k1_0, k1_0});
//   ParameterSpace2d parameter_space{{kKnotVector, kKnotVector}, {Degree{2},
//   Degree{2}}}}; constexpr ParametricCoordinate const k0_5{0.5};
//   ParameterSpace2d::ParametricCoordinate_ const kParametricCoordinate{k0_5,
//   k0_5};  // Set u = {0.5, 0.5}.
//   // Find the first basis function whose support contains u = {0.5, 0.5}.
//   ParameterSpace2d::Index_ const &non_zero =
//   parameter_space.FindFirstNonZeroBasisFunction(kParametricCoordinate);
//   // Evaluate N_{0,0}^{2,2}, i.e., the first basis function that is non-zero
//   for u = {0.5, 0.5}. ParameterSpace2d::Type_ const &evaluated =
//   parameter_space_.EvaluateBasisFunction(non_zero, kParametricCoordinate);
//   parameter_space.RemoveKnot(Dimension{}, ParametricCoordinate{1.0});  // Try
//   to remove u = {1.0} from U_0. parameter_space.ElevateDegree(Dimension{1});
//   // Elevate the degree p for the second parametric dimension by one.
template<int parametric_dimensionality>
class ParameterSpace {
private:
  using StringArray_ = StringArray<parametric_dimensionality>;

public:
  using BinomialRatios_ = Vector<BinomialRatio>;
  using Degrees_ = Array<Degree, parametric_dimensionality>;
  using Derivative_ = Array<Derivative, parametric_dimensionality>;
  using ElevationCoefficients_ = Vector<BinomialRatios_>;
  using ElevationInformation_ = Tuple<Index, ElevationCoefficients_>;
  using Index_ = utilities::Index<parametric_dimensionality>;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  using KnotRatios_ = Vector<KnotRatio>;
  using KnotVectors_ = KnotVectors<parametric_dimensionality>;
  using NumberOfBasisFunctions_ = IndexLength_;
  using NumberOfParametricCoordinates_ = IndexLength_;
  using OutputInformation_ =
      Tuple<KnotVectorsOutputInformation<parametric_dimensionality>,
            StringArray_,
            StringArray_>;
  using ParametricCoordinate_ =
      Array<ParametricCoordinate, parametric_dimensionality>;
  using ParametricCoordinates_ = Vector<ParametricCoordinate_>;
  using Type_ = Type;
  using InsertionCoefficients_ = Vector<KnotRatios_>;
  using InsertionInformation_ = Tuple<Index, InsertionCoefficients_>;
  using Knots_ = typename KnotVectors_::value_type::element_type::Knots_;
  using BezierInformation_ = Tuple<int, Knots_>;
  using Knot_ = typename Knots_::value_type;
  using KnotSpans_ = Array<KnotSpan, parametric_dimensionality>;

  // for evaluated basis values
  using BasisValues_ = Vector<Type_>;
  using BasisValueType_ = typename BasisValues_::value_type;
  using BasisValuesPerDimension_ =
      Array<BasisValues_, parametric_dimensionality>;

  ParameterSpace() = default;
  ParameterSpace(KnotVectors_ knot_vectors,
                 Degrees_ degrees,
                 Tolerance const& tolerance = kEpsilon);
  ParameterSpace(ParameterSpace const& other);
  ParameterSpace(ParameterSpace&& other) noexcept = default;
  ParameterSpace& operator=(ParameterSpace const& rhs);
  ParameterSpace& operator=(ParameterSpace&& rhs) noexcept = default;
  virtual ~ParameterSpace() = default;

  // Comparison based on tolerance.
  friend bool IsEqual<parametric_dimensionality>(ParameterSpace const& lhs,
                                                 ParameterSpace const& rhs,
                                                 Tolerance const& tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==<parametric_dimensionality>(ParameterSpace const& lhs,
                                                    ParameterSpace const& rhs);
  virtual Index_ First() const;
  virtual Index_ Behind() const;

  virtual NumberOfBasisFunctions_ GetNumberOfBasisFunctions() const;
  virtual int GetTotalNumberOfBasisFunctions() const;
  // Number of non-zero basis functions is equal to p+1 - see NURBS book P2.2.
  NumberOfBasisFunctions_ GetNumberOfNonZeroBasisFunctions() const;

  virtual SharedPointer<ParameterSpace<parametric_dimensionality - 1>>
  RemoveOneParametricDimension(const int parametric_dimension) const;

  virtual Index_ FindFirstNonZeroBasisFunction(
      ParametricCoordinate_ const& parametric_coordinate,
      Tolerance const& tolerance = kEpsilon) const;

  virtual KnotSpans_
  FindKnotSpans(ParametricCoordinate_ const& parametric_coordinate,
                Tolerance const& tolerance = kEpsilon) const;
  virtual BezierInformation_
  DetermineBezierExtractionKnots(Dimension const& dimension,
                                 Tolerance const& tolerance = kEpsilon) const;

  /// @brief Implements The NURBS Book A2.2
  /// @param parametric_coordinate
  /// @param tolerance
  /// @return
  virtual BasisValuesPerDimension_ EvaluateBasisValuesPerDimension(
      ParametricCoordinate_ const& parametric_coordinate,
      Tolerance const& tolerance = kEpsilon) const;

  virtual BasisValues_
  EvaluateBasisValues(ParametricCoordinate_ const& parametric_coordinate,
                      Tolerance const& tolerance = kEpsilon) const;

  /// @brief Implements The NURBS Book A2.3
  /// @param parametric_coordinate
  /// @param derivative
  /// @param tolerance
  /// @return
  virtual BasisValuesPerDimension_ EvaluateBasisDerivativeValuesPerDimension(
      ParametricCoordinate_ const& parametric_coordinate,
      Derivative_ const& derivative,
      Tolerance const& tolerance = kEpsilon) const;

  virtual BasisValues_ EvaluateBasisDerivativeValues(
      ParametricCoordinate_ const& parametric_coordinate,
      Derivative_ const& derivative,
      Tolerance const& tolerance = kEpsilon) const;

  virtual InsertionInformation_
  InsertKnot(Dimension const& dimension,
             Knot_ knot,
             Multiplicity const& multiplicity = kMultiplicity,
             Tolerance const& tolerance = kEpsilon);

  // Tries to interpret knot removal as the inverse process of knot insertion.
  virtual InsertionInformation_
  RemoveKnot(Dimension const& dimension,
             Knot_ const& knot,
             Multiplicity const& multiplicity = kMultiplicity,
             Tolerance const& tolerance = kEpsilon);

  // Elevation information refer to degree elevation of Bezier basis functions.
  virtual ElevationInformation_
  ElevateDegree(Dimension const& dimension,
                Multiplicity const& multiplicity = kMultiplicity,
                Tolerance const& tolerance = kEpsilon);

  // Tries to interpret degree reduction as the inverse process of degree
  // elevation.
  virtual ElevationInformation_
  ReduceDegree(Dimension const& dimension,
               Multiplicity const& multiplicity = kMultiplicity,
               Tolerance const& tolerance = kEpsilon);

  virtual ParametricCoordinates_
  Sample(NumberOfParametricCoordinates_ const& number_of_parametric_coordinates)
      const;
  virtual OutputInformation_
  Write(Precision const& precision = kPrecision) const;

#ifndef NDEBUG
  virtual void ThrowIfParametricCoordinateIsOutsideScope(
      Dimension const& dimension,
      ParametricCoordinate const& parametric_coordinate,
      Tolerance const& tolerance = kEpsilon) const;
#endif

  virtual const KnotVectors_& GetKnotVectors() const { return knot_vectors_; }
  virtual const Degrees_& GetDegrees() const { return degrees_; }
  virtual KnotVectors_& GetKnotVectors() { return knot_vectors_; }
  virtual Degrees_& GetDegrees() { return degrees_; }

protected:
  KnotVectors_ knot_vectors_;
  Degrees_ degrees_;

private:
  using MultiplicityType_ = Multiplicity::Type_;

  void CopyKnotVectors(KnotVectors_ const& knot_vectors);

  // Number of non-zero basis functions is equal to p+1 - see NURBS book P2.2.
  int GetNumberOfNonZeroBasisFunctions(Dimension const& dimension) const;
  // Number of basis functions n is equal to m-(p+1) - see NURBS book P2.8.
  Length GetNumberOfBasisFunctions(Dimension const& dimension) const;

  InsertionInformation_ DetermineInsertionInformation(
      Dimension const& dimension,
      ParametricCoordinate const& knot,
      Multiplicity const& multiplicity = kMultiplicity,
      Tolerance const& tolerance = kEpsilon) const;
  ElevationInformation_ DetermineElevationInformation(
      Dimension const& dimension,
      Multiplicity const& multiplicity = kMultiplicity) const;

#ifndef NDEBUG
  void
  ThrowIfBasisFunctionIndexIsInvalid(Index_ const& basis_function_index) const;
  void ThrowIfFrontOrBackKnotIsToBeInsertedOrRemoved(
      Dimension const& dimension,
      ParametricCoordinate const& knot,
      Tolerance const& tolerance) const;
#endif
};

#include "BSplineLib/ParameterSpaces/parameter_space.inc"

} // namespace bsplinelib::parameter_spaces

#endif // SOURCES_PARAMETERSPACES_PARAMETER_SPACE_HPP_
