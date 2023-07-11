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
WeightedVectorSpace<dimensionality>::WeightedVectorSpace(
    Coordinates_ const& coordinates,
    Weights_ const& weights)
    : Base_{HomogenizeCoordinates(coordinates, weights)} {}

template<int dimensionality>
bool IsEqual(WeightedVectorSpace<dimensionality> const& lhs,
             WeightedVectorSpace<dimensionality> const& rhs,
             Tolerance const& tolerance) {
  using Base = typename WeightedVectorSpace<dimensionality>::Base_;

#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::IsEqual::WeightedVectorSpace");
  }
#endif
  return IsEqual(static_cast<Base const&>(lhs),
                 static_cast<Base const&>(rhs),
                 tolerance);
}

template<int dimensionality>
bool operator==(WeightedVectorSpace<dimensionality> const& lhs,
                WeightedVectorSpace<dimensionality> const& rhs) {
  return IsEqual(lhs, rhs);
}

template<int dimensionality>
typename WeightedVectorSpace<dimensionality>::Coordinate_
WeightedVectorSpace<dimensionality>::CreateCoordinate() const {
  return Coordinate_(dimensionality);
}

template<int dimensionality>
typename WeightedVectorSpace<dimensionality>::Coordinate_
WeightedVectorSpace<dimensionality>::Project(
    HomogeneousCoordinate_ const& homogeneous_coordinate) {

  Coordinate_ coordinate;
  coordinate.reserve(dimensionality);
  const auto w_inv = 1 / homogeneous_coordinate[dimensionality];
  for (int i{}; i < dimensionality; ++i) {
    coordinate.emplace_back(homogeneous_coordinate[i] * w_inv);
  }

  return coordinate;
}

template<int dimensionality>
typename WeightedVectorSpace<dimensionality>::Coordinate_
WeightedVectorSpace<dimensionality>::Project(
    ConstCoordinateView const& homogeneous_coordinate) {
  Coordinate_ coordinate;
  coordinate.reserve(dimensionality);
  const auto w_inv = 1 / homogeneous_coordinate[dimensionality];
  for (int i{}; i < dimensionality; ++i) {
    coordinate.emplace_back(homogeneous_coordinate[i] * w_inv);
  }

  return coordinate;
}

template<int dimensionality>
typename WeightedVectorSpace<
    dimensionality>::MaximumDistanceFromOriginAndMinimumWeight_
WeightedVectorSpace<dimensionality>::
    DetermineMaximumDistanceFromOriginAndMinimumWeight(
        Tolerance const& tolerance) const {

  double max_distance{}, min_weight{};
  for (int i{}; i < Base_::GetNumberOfCoordinates(); ++i) {
    const auto hc_view = Base_::operator[](0);
    max_distance =
        std::max(utilities::std_container_operations::TwoNorm(Project(hc_view)),
                 max_distance);
    min_weight = std::min(hc_view[dimensionality], min_weight);
  }

  return typename WeightedVectorSpace<
      dimensionality>::MaximumDistanceFromOriginAndMinimumWeight_{
      Coordinate{max_distance},
      Coordinate{min_weight}};
}

template<int dimensionality>
typename WeightedVectorSpace<dimensionality>::HomogeneousCoordinates_
WeightedVectorSpace<dimensionality>::HomogenizeCoordinates(
    Coordinates_ const& coordinates,
    Weights_ const& weights) const {
  using std::to_string;

  int const& number_of_coordinates = coordinates.size() / dimensionality;
#ifndef NDEBUG
  int const& number_of_weights = weights.size();
  if (number_of_coordinates != number_of_weights)
    throw DomainError(
        to_string(number_of_weights) + " weights were provided but "
        + to_string(number_of_coordinates)
        + " are needed to associate each weight with a coordinate.");
#endif
  HomogeneousCoordinates_ homogeneous_coordinates;
  homogeneous_coordinates.reserve(number_of_coordinates * (dimensionality + 1));

  const auto* current_coord = coordinates.data();
  int c{};
  for (int i{}; i < number_of_coordinates; ++i) {
    const auto& w = weights[i];
    for (int j{}; j < dimensionality; ++j) {
      homogeneous_coordinates.emplace_back(*(current_coord++) * w);
    }

    homogeneous_coordinates.push_back(w);
    c++;
  }

  return homogeneous_coordinates;
}
