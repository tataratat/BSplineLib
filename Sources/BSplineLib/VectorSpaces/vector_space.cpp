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

void VectorSpace::Replace(int const& coordinate_index,
                          const Coordinate_& coordinate) {
  std::copy_n(coordinate.begin(),
              coordinate.size(),
              &coordinates_(coordinate_index, 0));
}

void VectorSpace::ReallocateInsert(int const& coordinate_index,
                                   const Coordinate_& coordinate) {
  // This is a lot of copy
  const auto& shape = coordinates_.Shape();
  const auto& n_coord = shape[0];
  const auto& dim = shape[1];

  assert(dim == coordinate.size());

  Coordinates_ new_coordinates(n_coord + 1, dim);

  // copy (index) elements
  std::copy_n(coordinates_.begin(),
              coordinate_index * dim,
              new_coordinates.begin());

  // copy at index
  std::copy_n(coordinate.begin(), dim, &new_coordinates(coordinate_index, 0));

  // copy after index
  std::copy(&coordinates_(coordinate_index, 0),
            coordinates_.end(),
            &new_coordinates(coordinate_index + 1, 0));

  // move assign new coords as coords
  coordinates_ = std::move(new_coordinates);
}

void VectorSpace::Erase(int const& coordinate_index) {
  // we just need to "shorten" data at erase space
  std::copy(&coordinates_(coordinate_index + 1, 0),
            coordinates_.end(),
            &coordinates_(coordinate_index, 0));

  // adjust shape only.
  coordinates_.SetShape(coordinates_.Shape()[0] - 1, coordinates_.Shape()[1]);
}

typename VectorSpace::DataType_
VectorSpace::DetermineMaximumDistanceFromOrigin() const {
  Coordinate maximum_distance{};
  const auto& n_coords = coordinates_.Shape()[0];
  const auto& dim = coordinates_.Shape()[1];

  ConstCoordinate_ view;
  view.SetShape(dim);

  for (int i{}; i < n_coords; ++i) {
    view.SetData(CoordinateBegin(i));
    maximum_distance = std::max(view.NormL2(), maximum_distance);
  }

  return maximum_distance;
}

typename VectorSpace::OutputInformation_
VectorSpace::Write(Precision const& precision) const {
  // until we move iges to python, we create a type matching copy here.
  const int n = coordinates_.Shape()[0];
  const int d = coordinates_.Shape()[1];

  const DataType_* coord_ptr = coordinates_.begin();
  Vector<Vector<DataType_>> nested_coordinates(n);
  for (auto& nc : nested_coordinates) {
    nc.resize(d);
    for (auto& c : nc) {
      c = *coord_ptr++;
    }
  }

  return OutputInformation_{utilities::string_operations::Write<
      std::tuple_element_t<0, OutputInformation_>>(nested_coordinates,
                                                   precision)};
}

} // namespace bsplinelib::vector_spaces
