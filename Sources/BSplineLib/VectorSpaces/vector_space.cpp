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

#include "BSplineLib/VectorSpaces/vector_space.hpp"

namespace bsplinelib::vector_spaces {

VectorSpace::VectorSpace(Coordinates_ coordinates, const int dimensionality)
    : coordinates_(std::move(coordinates)) {
  SetDimensionality(dimensionality);
}

bool IsEqual(VectorSpace const& lhs,
             VectorSpace const& rhs,
             Tolerance const& tolerance) {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::IsEqual::VectorSpace");
  }
#endif
  return std::equal(
      lhs.coordinates_.begin(),
      lhs.coordinates_.end(),
      rhs.coordinates_.begin(),
      rhs.coordinates_.end(),
      [&tolerance](const VectorSpace::Coordinates_::value_type& a,
                   const VectorSpace::Coordinates_::value_type& b) {
        return std::abs(a - b) < tolerance;
      });
}

bool operator==(VectorSpace const& lhs, VectorSpace const& rhs) {
  return IsEqual(lhs, rhs);
}

void VectorSpace::SetDimensionality(const int dimensionality) {
  // check if this is a valid number
  const int coord_size = GetCoordinates().size();
  const auto [quot, rem] = std::div(coord_size, dimensionality);

  if (rem != 0) {
    throw std::runtime_error("Invalid dimensionality. Size of the coordinates "
                             "must be a multiple of dimensionality.");
  }

  // good. set
  dimensionality_ = dimensionality;
}

int VectorSpace::GetDimensionality() const {
  if (dimensionality_ < 1) {
    throw std::runtime_error(
        "Dimensionality not set. Please set dimensionality. Currently set to "
        + std::to_string(dimensionality_));
  }
  return dimensionality_;
}

CoordinateView VectorSpace::operator[](const int coordinate_id) {
  const int dim = GetDimensionality();
  double* begin_ptr = GetCoordinates().data();

  return {begin_ptr + (coordinate_id * dim), dim};
}

CoordinateView VectorSpace::operator[](Index const& coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::VectorSpace::operator[]");
  }
#endif
  return operator[](coordinate.Get());
}

ConstCoordinateView VectorSpace::operator[](const int coordinate_id) const {
  const int dim = GetDimensionality();
  const double* begin_ptr = GetCoordinates().data();

  return {begin_ptr + (coordinate_id * dim), dim};
}
CoordinateView VectorSpace::DirectView(const int coordinate_id) {
  return {coordinates_.data() + (coordinate_id * dimensionality_),
          dimensionality_};
}

ConstCoordinateView VectorSpace::DirectView(const int coordinate_id) const {
  return {coordinates_.data() + (coordinate_id * dimensionality_),
          dimensionality_};
}

typename VectorSpace::Coordinate_ VectorSpace::CreateCoordinate() const {
  return Coordinate_(GetDimensionality());
}

typename VectorSpace::Coordinate_
VectorSpace::CreateCoordinate(std::initializer_list<double> coordinate) const {
  // size check
  assert(static_cast<int>(coordinate.size())
         == GetDimensionality()); // coordinate size and dim mismatch

  return Coordinate_(coordinate);
}

typename VectorSpace::Coordinates_& VectorSpace::GetCoordinates() {
  // if (coordinates_.size() == 0) {
  //   throw std::runtime_error(
  //       "Empty Coordinates. Please set coordinates first.");
  // }
  return coordinates_;
}

typename VectorSpace::Coordinates_ const& VectorSpace::GetCoordinates() const {
  if (coordinates_.size() == 0) {
    throw std::runtime_error(
        "Empty Coordinates. Please set coordinates first.");
  }
  return coordinates_;
}

int VectorSpace::GetNumberOfCoordinates() const {
  return GetCoordinates().size() / GetDimensionality();
}

void VectorSpace::Erase(Index const& coordinate_index) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::VectorSpace::Erase");
  }
#endif
  const int dim = GetDimensionality();
  const auto begin = coordinates_.begin() + (coordinate_index.Get() * dim);

  coordinates_.erase(begin, begin + dim);
}

Coordinate VectorSpace::DetermineMaximumDistanceFromOrigin(
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception,
          "bsplinelib::vector_spaces::VectorSpace::"
          "DetermineMaximumDistanceFromOrigin");
  }
#endif
  double maximum_distance{};
  for (int i{}; i < GetNumberOfCoordinates(); ++i) {
    const auto c_view = operator[](i);
    maximum_distance =
        std::max(utilities::container_operations::TwoNorm(c_view),
                 maximum_distance);
  }
  return Coordinate{maximum_distance};
}

#ifndef NDEBUG

void VectorSpace::ThrowIfIndexIsInvalid(Index const& coordinate) const {
  Index::ThrowIfNamedIntegerIsOutOfBounds(coordinate, coordinates_.size() - 1);
}
#endif

} // namespace bsplinelib::vector_spaces
