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

#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"
#include "BSplineLib/Utilities/std_container_operations.hpp"
#include "BSplineLib/VectorSpaces/vector_space.hpp"

namespace bsplinelib::vector_spaces {

// WeightedVectorSpaces store coordinates and weights together using
// homogeneous, i.e., weighted coordinates.
//
class WeightedVectorSpace : public VectorSpace {
public:
  using Base_ = VectorSpace;
  using DataType_ = typename Base_::DataType_;
  using Coordinate_ = typename Base_::Coordinate_;
  using ConstCoordinate_ = typename Base_::ConstCoordinate_;
  using Coordinates_ = typename Base_::Coordinates_;
  using HomogeneousCoordinate_ = typename Base_::Coordinate_;
  using MaximumDistanceFromOriginAndMinimumWeight_ = Tuple<Coordinate, Weight>;
  using Weights_ = Base_::Array_<Weight>;

  // Make Base Constructor public
  using Base_::Base_;

  WeightedVectorSpace() = default;
  WeightedVectorSpace(Coordinates_ const& coordinates, Weights_ const& weights);
  WeightedVectorSpace(WeightedVectorSpace const& other) = default;
  WeightedVectorSpace(WeightedVectorSpace&& other) noexcept = default;
  WeightedVectorSpace& operator=(WeightedVectorSpace const& rhs) = default;
  WeightedVectorSpace& operator=(WeightedVectorSpace&& rhs) noexcept = default;
  ~WeightedVectorSpace() override = default;

  virtual MaximumDistanceFromOriginAndMinimumWeight_
  DetermineMaximumDistanceFromOriginAndMinimumWeight(
      Tolerance const& tolerance = kEpsilon) const;

private:
  constexpr void HomogenizeCoordinates(Coordinates_ const& coordinates,
                                       Weights_ const& weights) const;
};

#include "BSplineLib/VectorSpaces/weighted_vector_space.inl"

} // namespace bsplinelib::vector_spaces

#endif // SOURCES_VECTORSPACES_WEIGHTED_VECTOR_SPACE_HPP_
