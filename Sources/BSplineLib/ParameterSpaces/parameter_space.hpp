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
#include "BSplineLib/Utilities/containers.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/index.hpp"
#include "BSplineLib/Utilities/math_operations.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"
#include "BSplineLib/Utilities/string_operations.hpp"

namespace bsplinelib::parameter_spaces {

/// @brief ParameterSpaces provide the B-spline basis functions corresponding to
/// given knot vectors and degrees.
/// @tparam para_dim
template<int para_dim>
class ParameterSpace {
private:
  using StringArray_ = StringArray<para_dim>;

public:
  using BinomialRatios_ = Vector<BinomialRatio>;
  using Degrees_ = Array<Degree, para_dim>;
  using Derivative_ = Array<Derivative, para_dim>;
  using ElevationCoefficients_ = Vector<BinomialRatios_>;
  using ElevationInformation_ = Tuple<Index, ElevationCoefficients_>;
  using Index_ = utilities::Index<para_dim>;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  using KnotRatios_ = Vector<KnotRatio>;
  using KnotVectors_ = KnotVectors<para_dim>;
  using NumberOfBasisFunctions_ = IndexLength_;
  using NumberOfParametricCoordinates_ = IndexLength_;
  using OutputInformation_ =
      Tuple<KnotVectorsOutputInformation<para_dim>, StringArray_, StringArray_>;
  using ParametricCoordinate_ = Array<ParametricCoordinate, para_dim>;
  using Type_ = Type;
  using IntType_ = Degree; // same as derivative
  using InsertionCoefficients_ = Vector<KnotRatios_>;
  using InsertionInformation_ = Tuple<Index, InsertionCoefficients_>;
  using Knots_ = typename KnotVectors_::value_type::element_type::Knots_;
  using BezierInformation_ = Tuple<int, Knots_>;
  using Knot_ = typename Knots_::value_type;
  using KnotSpans_ = Array<KnotSpan, para_dim>;

  // for evaluated basis values
  template<typename T>
  using TemporaryData_ = bsplinelib::utilities::containers::TemporaryData<T>;
  template<typename T>
  using TemporaryData2D_ =
      bsplinelib::utilities::containers::TemporaryData2D<T>;
  template<typename T>
  using Data_ = bsplinelib::utilities::containers::Data<T>;

  using BasisValues_ = Data_<Type_>;
  using BasisValueType_ = typename BasisValues_::value_type;
  using BasisValuesPerDimension_ = Array<BasisValues_, para_dim>;

  ParameterSpace() = default;
  ParameterSpace(KnotVectors_ knot_vectors, Degrees_ degrees)
      : knot_vectors_(std::move(knot_vectors)),
        degrees_(std::move(degrees)){};
  ParameterSpace(ParameterSpace const& other);
  ParameterSpace(ParameterSpace&& other) noexcept = default;
  ParameterSpace& operator=(ParameterSpace const& rhs);
  ParameterSpace& operator=(ParameterSpace&& rhs) noexcept = default;
  virtual ~ParameterSpace() = default;

  virtual Index_ First() const;
  virtual Index_ Behind() const;

  virtual NumberOfBasisFunctions_ GetNumberOfBasisFunctions() const;
  virtual int GetTotalNumberOfBasisFunctions() const;
  // Number of non-zero basis functions is equal to p+1 - see NURBS book P2.2.
  NumberOfBasisFunctions_ GetNumberOfNonZeroBasisFunctions() const;

  virtual SharedPointer<ParameterSpace<para_dim - 1>>
  RemoveOneParametricDimension(const int parametric_dimension) const;

  virtual Index_
  FindFirstNonZeroBasisFunction(const Type_* parametric_coordinate,
                                Tolerance const& tolerance = kEpsilon) const;

  virtual KnotSpans_ FindKnotSpans(const Type_* parametric_coordinate,
                                   Tolerance const& tolerance = kEpsilon) const;
  virtual BezierInformation_
  DetermineBezierExtractionKnots(Dimension const& dimension,
                                 Tolerance const& tolerance = kEpsilon) const;

  /// @brief Return multiplicities of uninque_knots for each diemnsion
  /// @param tolerance
  virtual Vector<Vector<int>>
  KnotMultiplicities(Tolerance const& tolerance = kEpsilon) const;

  /// @brief Implements The NURBS Book A2.2
  /// @param parametric_coordinate
  /// @param tolerance
  /// @return
  virtual BasisValuesPerDimension_
  EvaluateBasisValuesPerDimension(const Type_* parametric_coordinate,
                                  Tolerance const& tolerance = kEpsilon) const;

  virtual BasisValues_
  EvaluateBasisValues(const Type_* parametric_coordinate,
                      Tolerance const& tolerance = kEpsilon) const;

  /// @brief Implements The NURBS Book A2.3
  /// @param parametric_coordinate
  /// @param derivative
  /// @param tolerance
  /// @return
  virtual BasisValuesPerDimension_ EvaluateBasisDerivativeValuesPerDimension(
      const Type_* parametric_coordinate,
      const IntType_* derivative,
      Tolerance const& tolerance = kEpsilon) const;

  virtual BasisValues_
  EvaluateBasisDerivativeValues(const Type_* parametric_coordinate,
                                const IntType_* derivative,
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

/// @brief shortcut to ParameterSpace<>::BasisValues_
using BasisValues = bsplinelib::utilities::containers::Data<Type>;

/// recursive combime adapted from bezman
template<std::size_t depth, std::size_t array_dim>
constexpr void RecursiveCombine_(const Array<BasisValues, array_dim>& factors,
                                 BasisValues& result,
                                 int& counter,
                                 const Type& c_value) {
  static_assert(depth < array_dim,
                "Implementation error, recursion loop to deep!");

  for (int i{}; i < factors[depth].size(); ++i) {
    if constexpr (depth == 0) {
      // compute result and advance counter
      result[counter++] = c_value * factors[depth][i];
    } else {
      RecursiveCombine_<static_cast<std::size_t>(depth - 1)>(
          factors,
          result,
          counter,
          c_value * factors[depth][i]);
    }
  }
}

/// recursive combime adapted from bezman
template<std::size_t array_dim>
BasisValues RecursiveCombine(const Array<BasisValues, array_dim>& factors) {
  // Precalculate required entries
  int n_entries{1};
  std::for_each(factors.begin(), factors.end(), [&](const BasisValues& ii) {
    n_entries *= ii.size();
  });

  // Init return type and reserve memory
  BasisValues result(n_entries);
  int counter{};

  // Start computation
  RecursiveCombine_<array_dim - 1>(factors, result, counter, 1.);

  assert(counter == n_entries);

  return result;
}

/// recursive combime taken from bezman
template<std::size_t depth,
         std::size_t array_dim,
         typename ValueType,
         typename BSplineLibIndex,
         typename CoeffType,
         typename ReturnType>
constexpr void RecursiveCombine_(const Array<BasisValues, array_dim>& factors,
                                 BSplineLibIndex& index,
                                 BSplineLibIndex& index_offset,
                                 const CoeffType& coeffs,
                                 const ValueType& c_value,
                                 ReturnType& result) {
  static_assert(depth < array_dim,
                "Implementation error, recursion loop to deep!");

  for (const auto& factor : factors[depth]) {
    if constexpr (depth == 0) {
      // get basis contribution
      const auto fac = c_value * factor;
      // add index offset and get beginning of corresponding control point
      const auto& index_offset_multi_index = index_offset.MultiIndex();
      index += index_offset_multi_index;
      const auto* coeff = &coeffs(index.GetIndex1d().Get(), 0);

      // contribute to each dim
      result.Add(fac, coeff);

      // setback index offset
      index -= index_offset_multi_index;
      ;

      // advance to the next index offset - this needs to be done after offset
      // setback
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
template<std::size_t array_dim,
         typename BSplineLibIndex,
         typename CoeffType,
         typename ReturnType>
constexpr void RecursiveCombine(const Array<BasisValues, array_dim>& factors,
                                BSplineLibIndex& index,
                                BSplineLibIndex& index_offset,
                                const CoeffType& coeffs,
                                ReturnType& result) {

  // Start computation
  RecursiveCombine_<array_dim - 1>(factors,
                                   index,
                                   index_offset,
                                   coeffs,
                                   1.,
                                   result);
}

#include "BSplineLib/ParameterSpaces/parameter_space.inl"

} // namespace bsplinelib::parameter_spaces

#endif // SOURCES_PARAMETERSPACES_PARAMETER_SPACE_HPP_
