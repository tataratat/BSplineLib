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

#ifndef SOURCES_VECTORSPACES_WEIGHTED_VECTOR_SPACE_HPP_
#define SOURCES_VECTORSPACES_WEIGHTED_VECTOR_SPACE_HPP_

#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>

#include "BSplineLib/Utilities/container_operations.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"
#include "BSplineLib/VectorSpaces/vector_space.hpp"

namespace bsplinelib::vector_spaces {

class WeightedVectorSpace;

bool IsEqual(WeightedVectorSpace const& lhs,
             WeightedVectorSpace const& rhs,
             Tolerance const& tolerance = kEpsilon);

bool operator==(WeightedVectorSpace const& lhs, WeightedVectorSpace const& rhs);

// WeightedVectorSpaces store coordinates and weights together using
// homogeneous, i.e., weighted coordinates.
class WeightedVectorSpace : public VectorSpace {
private:
  using VectorSpace_ = VectorSpace;

public:
  using Base_ = VectorSpace;
  using Coordinate_ = typename VectorSpace_::Coordinate_;
  using Coordinates_ = typename VectorSpace_::Coordinates_;
  using NamedCoordinate_ = Vector<Coordinate>;
  using NestedCoordinates_ = Vector<Coordinate_>;
  using HomogeneousCoordinate_ = typename Base_::Coordinate_;
  using MaximumDistanceFromOriginAndMinimumWeight_ = Tuple<Coordinate, Weight>;
  using Weights_ = Vector<double>;

  // Make Base Constructor public
  using Base_::Base_;

  WeightedVectorSpace() = default;
  WeightedVectorSpace(Coordinates_ const& coordinates,
                      Weights_ const& weights,
                      const int dimensionality);
  WeightedVectorSpace(WeightedVectorSpace const& other) = default;
  WeightedVectorSpace(WeightedVectorSpace&& other) noexcept = default;
  WeightedVectorSpace& operator=(WeightedVectorSpace const& rhs) = default;
  WeightedVectorSpace& operator=(WeightedVectorSpace&& rhs) noexcept = default;
  ~WeightedVectorSpace() override = default;

  // Comparison based on tolerance.
  friend bool IsEqual(WeightedVectorSpace const& lhs,
                      WeightedVectorSpace const& rhs,
                      Tolerance const& tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==(WeightedVectorSpace const& lhs,
                         WeightedVectorSpace const& rhs);

  virtual void SetDimensionality(const int dimensionality) override;
  virtual int GetProjectedDimensionality() const;
  // virtual Coordinate_ CreateCoordinate() const override;
  static Coordinate_
  Project(HomogeneousCoordinate_ const& homogeneous_coordinate);
  static Coordinate_ Project(ConstCoordinateView const& homogeneous_coordinate);
  virtual MaximumDistanceFromOriginAndMinimumWeight_
  DetermineMaximumDistanceFromOriginAndMinimumWeight(
      Tolerance const& tolerance = kEpsilon) const;
  virtual void SyncCoordinatesSizesAndDims(const bool to_homogeneous);
  virtual void SyncCoordinates(const bool to_homogeneous);

  // template<typename CoordinateType>
  // void Replace(Index const& coordinate_index, CoordinateType coordinate);
  // template<typename CoordinateType>
  // void Insert(Index const& coordinate_index, CoordinateType coordinate);

protected:
  VectorSpace_ unweighted_space_; // this dim - 1
  Coordinates_ weights_;

private:
  using HomogeneousCoordinates_ = typename Base_::Coordinates_;

  HomogeneousCoordinates_ HomogenizeCoordinates(Coordinates_ const& coordinates,
                                                Weights_ const& weights,
                                                const int dimensionality) const;
};

} // namespace bsplinelib::vector_spaces

#endif // SOURCES_VECTORSPACES_WEIGHTED_VECTOR_SPACE_HPP_
