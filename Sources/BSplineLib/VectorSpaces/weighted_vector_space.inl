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

WeightedVectorSpace::WeightedVectorSpace(Coordinates_ const& coordinates,
                                         Weights_ const& weights)
    : {
  HomogenizeCoordinates(coordinates, weights);
}

void WeightedVectorSpace::Project(
    HomogeneousCoordinate_ const& homogeneous_coordinate,
    Coordinate_& projected) {

  assert(homogeneous_coordinate.size() == projected.size() + 1);

  Project(homogeneous_coordinate, projected.data());
}

void WeightedVectorSpace::Project(
    HomogeneousCoordinate_ const& homogeneous_coordinate,
    DataType_* projected) {

  const int dim = homogeneous_coordinate.size() - 1;

  const DataType_ w_inv =
      static_cast<DataType_>(1.0) / homogeneous_coordinate[dim];

  for (int i{}; i < dim; ++i) {
    projected[i] = homogeneous_coordinate[i] * w_inv;
  }
}

typename WeightedVectorSpace::MaximumDistanceFromOriginAndMinimumWeight_
WeightedVectorSpace::DetermineMaximumDistanceFromOriginAndMinimumWeight()
    const {
  Coordinate maximum_distance{};
  Weight minimum_weight{std::numeric_limits<Weight>::max()};

  const auto& n_coords = coordinates_.Shape()[0];
  const auto& h_dim = coordinates_.Shape()[1];

  // get a view, excluding weight
  ConstCoordinate_ view;
  view.SetShape(h_dim - 1);
  for (int i{}; i < n_coords; ++i) {
    view.SetData(Base_::CoordinateBegin(i));
    maximum_distance = std::max(view.NormL2(), maximum_distance);
    // typically end() means one after the last valid element, but in this case
    // this is still valid, because end() will point to the weight portion of
    // the homogeneous coordinate
    minimum_weight = std::min(*view.end(), minimum_weight);
  }

  return {maximum_distance, minimum_weight};
}

void WeightedVectorSpace::HomogenizeCoordinates(Coordinates_ const& coordinates,
                                                Weights_ const& weights) const {
  using std::to_string;

  const int number_of_coordinates = coordinates.Shape()[0];
  const int dim = coordinates.Shape()[1];
  const int number_of_weights = weights.size();

  // size check
  if (number_of_coordinates != number_of_weights)
    throw DomainError(
        to_string(number_of_weights) + " weights were provided but "
        + to_string(number_of_coordinates)
        + " are needed to associate each weight with a coordinate.");

  // those are homogenized coordinates
  // first set shape -> this computes size
  Base_::coordinates_.SetShape(number_of_coordinates, dim + 1);
  Base_::coordinates_.Reallocate(Base_::coordinates_.size());

  auto* h_coord = Base_::coordinates_.begin();
  const auto* coord = coordinates.begin();
  const auto* weight = weights.begin();
  for (int i{}; i < number_of_coordinates; ++i) {
    // apply weight
    const auto& w = *weight; // we deref here just to avoid (* *weight)
    for (int j{}; j < dim; ++j) {
      *h_coord++ = *coord++ * w;
    }
    // assign weight
    *h_coord++ = *weight++;
  }
}
