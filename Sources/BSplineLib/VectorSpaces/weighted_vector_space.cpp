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

#include "BSplineLib/VectorSpaces/weighted_vector_space.hpp"

namespace bsplinelib::vector_spaces {

WeightedVectorSpace::WeightedVectorSpace(Coordinates_ const& coordinates,
                                         Weights_ const& weights) {
  HomogenizeCoordinates(coordinates, weights);
}

typename WeightedVectorSpace::Coordinate_ WeightedVectorSpace::Project(
    HomogeneousCoordinate_ const& homogeneous_coordinate) {

  const int dim = homogeneous_coordinate.size() - 1;
  Coordinate_ projected(dim);

  const DataType_ w_inv =
      static_cast<DataType_>(1.0) / homogeneous_coordinate[dim];

  projected.Multiply(w_inv);

  return projected;
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

constexpr void
WeightedVectorSpace::HomogenizeCoordinates(Coordinates_ const& coordinates,
                                           Weights_ const& weights) {
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

typename WeightedVectorSpace::OutputInformation_
WeightedVectorSpace::WriteProjected(Precision const& precision) const {
  using std::tuple_element_t;
  using ProjectedCoordinatesOutput = tuple_element_t<0, OutputInformation_>;
  using utilities::string_operations::Write;

  HomogeneousCoordinates_ const& homogeneous_coordinates = Base_::coordinates_;
  const int n = homogeneous_coordinates.Shape()[0];
  const int d = homogeneous_coordinates.Shape()[1] - 1;

  ProjectedCoordinatesOutput coordinates(n);

  tuple_element_t<1, OutputInformation_> weights(n);
  for (int i{}; i < n; ++i) {
    // get beggining of the homogeneous coordinate
    const double* h_coord = &homogeneous_coordinates(i, 0);

    // get weight
    const double& w = h_coord[d];
    const double w_inv = 1. / w;

    // save weight
    weights[i] = Write(w, precision);

    // prepare projected coordinate - this is a string vector
    typename ProjectedCoordinatesOutput::value_type& projected_str =
        coordinates[i];
    projected_str.resize(d);

    // project coord and save
    for (auto& ps : projected_str) {
      ps = Write(*h_coord++ * w_inv, precision);
    }
  }

  return OutputInformation_{coordinates, weights};
}

} // namespace bsplinelib::vector_spaces
