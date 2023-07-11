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

template<int dimensionality>
VectorSpace<dimensionality>::VectorSpace(Coordinates_ coordinates)
    : coordinates_(std::move(coordinates)) {}

template<int dimensionality>
bool IsEqual(VectorSpace<dimensionality> const& lhs,
             VectorSpace<dimensionality> const& rhs,
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
      std::bind(utilities::std_container_operations::DoesContainEqualValues<
                    typename VectorSpace<dimensionality>::Coordinate_>,
                std::placeholders::_1,
                std::placeholders::_2,
                tolerance));
}

template<int dimensionality>
bool operator==(VectorSpace<dimensionality> const& lhs,
                VectorSpace<dimensionality> const& rhs) {
  return IsEqual(lhs, rhs);
}

template<int dimensionality>
typename VectorSpace<dimensionality>::Coordinate_
VectorSpace<dimensionality>::CreateCoordinate() const {
  return Coordinate_(dim_);
}

template<int dimensionality>
typename VectorSpace<dimensionality>::Coordinate_
VectorSpace<dimensionality>::CreateCoordinate(std::initializer_list<double> coordinate) const {
  // size check
  assert(coordinate.size() == dim_); // coordinate size and dim mismatch
  return Coordinate_(coordinate);
}

template<int dimensionality>
CoordinateView
VectorSpace<dimensionality>::operator[](const int coordinate_id) {
  double* begin_ptr = coordinates_.data();
  return {begin_ptr + (coordinate_id * dim_), dim_};
}

template<int dimensionality>
CoordinateView //const&
VectorSpace<dimensionality>::operator[](Index const& coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate);
  } catch (OutOfRange const& exception) {
    Throw(exception,
          "bsplinelib::vector_spaces::VectorSpace::operator[]");
  }
#endif
  return operator[](coordinate.Get());
}

template<int dimensionality>
ConstCoordinateView //const&
VectorSpace<dimensionality>::operator[](const int coordinate_id) const {
  const double* begin_ptr = coordinates_.data();
  return {begin_ptr + (coordinate_id * dim_), dim_};
}

template<int dimensionality>
int VectorSpace<dimensionality>::GetNumberOfCoordinates() const {
  return coordinates_.size() / dim_;
}

template<int dimensionality>
template<typename CoordinateType>
void VectorSpace<dimensionality>::Replace(Index const& coordinate_index,
                                          CoordinateType coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::VectorSpace::Replace");
  }
#endif
  //coordinates_[coordinate_index.Get()] = std::move(coordinate);
  const auto this_coordinate = operator[](coordinate_index);
  for (int i{}; i < dim_; ++i) {
    this_coordinate[i] = coordinate[i];
  }
}

template<int dimensionality>
template<typename CoordinateType>
void VectorSpace<dimensionality>::Insert(Index const& coordinate_index,
                                         CoordinateType coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::VectorSpace::Insert");
  }
#endif

  coordinates_.insert(coordinates_.begin() + (coordinate_index.Get() * dim_),
                      coordinate.begin(), coordinate.end());
}

template<int dimensionality>
void VectorSpace<dimensionality>::Erase(Index const& coordinate_index) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::VectorSpace::Erase");
  }
#endif
  const auto begin = coordinates_.begin() + (coordinate_index.Get() * dim_);
  coordinates_.erase(begin, begin + dim_);
}

template<int dimensionality>
Coordinate VectorSpace<dimensionality>::DetermineMaximumDistanceFromOrigin(
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
    maximum_distance = std::max(
        utilities::std_container_operations::TwoNorm(c_view),
        maximum_distance
    );
  }
  return Coordinate{maximum_distance};
}

#ifndef NDEBUG
template<int dimensionality>
void VectorSpace<dimensionality>::ThrowIfIndexIsInvalid(
    Index const& coordinate) const {
  Index::ThrowIfNamedIntegerIsOutOfBounds(coordinate, coordinates_.size() - 1);
}
#endif
