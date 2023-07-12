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
                                         Weights_ const& weights,
                                         const int dimensionality)
    : Base_(HomogenizeCoordinates(coordinates, weights, dimensionality),
            dimensionality + 1),
      unweighted_space_(coordinates, dimensionality),
      weights_(weights) {}

bool IsEqual(WeightedVectorSpace const& lhs,
             WeightedVectorSpace const& rhs,
             Tolerance const& tolerance) {
  using Base = typename WeightedVectorSpace::Base_;

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

bool operator==(WeightedVectorSpace const& lhs,
                WeightedVectorSpace const& rhs) {
  return IsEqual(lhs, rhs);
}

void WeightedVectorSpace::SetDimensionality(const int dimensionality) {
  // set projected dim
  unweighted_space_.SetDimensionality(dimensionality);
  // then normal dim
  Base_::SetDimensionality(dimensionality + 1);
}

int WeightedVectorSpace::GetProjectedDimensionality() const {
  try {
    return unweighted_space_.GetDimensionality();
  } catch (const std::runtime_error&) {
    throw std::runtime_error(
        "Projected Dimensionality not set. Please report this to "
        "github.com/tataratat/BSplineLib. Thank you!");
  }
}

typename WeightedVectorSpace::Coordinate_ WeightedVectorSpace::Project(
    HomogeneousCoordinate_ const& homogeneous_coordinate) {
  // prepare return
  Coordinate_ coordinate;
  const int p_dim = homogeneous_coordinate.size() - 1;
  coordinate.reserve(p_dim);

  // divide sum of weights
  const auto W_inv = 1 / homogeneous_coordinate[p_dim];
  for (int i{}; i < p_dim; ++i) {
    coordinate.emplace_back(homogeneous_coordinate[i] * W_inv);
  }

  return coordinate;
}

typename WeightedVectorSpace::Coordinate_ WeightedVectorSpace::Project(
    ConstCoordinateView const& homogeneous_coordinate) {
  // prepare return
  Coordinate_ coordinate;
  const int p_dim = homogeneous_coordinate.size() - 1;
  coordinate.reserve(p_dim);

  // divide sum of weights
  const auto W_inv = 1 / homogeneous_coordinate[p_dim];
  for (int i{}; i < p_dim; ++i) {
    coordinate.emplace_back(homogeneous_coordinate[i] * W_inv);
  }

  return coordinate;
}

typename WeightedVectorSpace::MaximumDistanceFromOriginAndMinimumWeight_
WeightedVectorSpace::DetermineMaximumDistanceFromOriginAndMinimumWeight(
    [[maybe_unused]] Tolerance const& tolerance) const {

  const int p_dim = GetProjectedDimensionality();

  double max_distance{}, min_weight{};
  for (int i{}; i < Base_::GetNumberOfCoordinates(); ++i) {
    const auto hc_view = Base_::operator[](0);
    max_distance =
        std::max(utilities::container_operations::TwoNorm(Project(hc_view)),
                 max_distance);
    min_weight = std::min(hc_view[p_dim], min_weight);
  }

  return
      typename WeightedVectorSpace::MaximumDistanceFromOriginAndMinimumWeight_{
          Coordinate{max_distance},
          Coordinate{min_weight}};
}

void WeightedVectorSpace::SyncCoordinatesSizesAndDims(
    const bool to_homogeneous) {
  if (to_homogeneous) {
    const int h_dim = GetDimensionality();
    const int h_size = GetNumberOfCoordinates();

    // first, weights
    weights_.resize(h_size);
    // then unweighted space
    unweighted_space_.GetCoordinates().resize(h_size * (h_dim - 1));
    unweighted_space_.SetDimensionality(h_dim - 1);
  } else {
    const int p_dim = unweighted_space_.GetDimensionality();
    const int p_size = unweighted_space_.GetNumberOfCoordinates();
    const int w_size = weights_.size();

    if (p_size != w_size) {
      throw std::runtime_error(
          "Count mismatch between unweighted coordinates and weights.");
    }

    // resize
    coordinates_.resize(p_size * (p_dim + 1));
    // set dim
    SetDimensionality(p_dim + 1);
  }
}

void WeightedVectorSpace::SyncCoordinates(const bool to_homogeneous) {
  // syncability check
  const int w_size = weights_.size();
  const int hc_size = GetNumberOfCoordinates();
  const int uwc_size = unweighted_space_.GetNumberOfCoordinates();
  if (!(w_size == hc_size == uwc_size)) {
    throw std::runtime_error("Count mismatch of weights, homogeneous "
                             "coordinates, and unweighted coordinates.");
  }

  const int projected_dim = GetProjectedDimensionality();
  if (to_homogeneous) {
    // unweighted to coordinates
    for (int i{}; i < hc_size; ++i) {
      // from
      const auto& w = weights_[i];
      const auto unweighted = unweighted_space_.DirectView(i);

      // to
      auto this_homogeneous = DirectView(i);

      // set
      for (int j{}; j < projected_dim; ++j) {
        this_homogeneous[j] = unweighted[j] * w;
      }
      this_homogeneous[projected_dim] = w;
    }
  } else {
    // homogeneous coordinates to unweighted
    for (int i{}; i < hc_size; ++i) {
      // from
      const auto this_homogeneous = DirectView(i);
      const auto& w = this_homogeneous[projected_dim];
      const auto w_inv = 1 / w;

      // to
      auto unweighted = unweighted_space_.DirectView(i);
      weights_[i] = w;

      // set
      for (int j{}; j < projected_dim; ++j) {
        unweighted[j] = this_homogeneous[j] * w_inv;
      }
    }
  }
}

// template<typename CoordinateType>
// void WeightedVectorSpace::Replace(Index const& coordinate_index,
//                           CoordinateType coordinate) {
//                             Base_::Replace(coordinate_index, coordinate);
//                           }
// template<typename CoordinateType>
// void WeightedVectorSpace::Insert(Index const& coordinate_index,
//                          CoordinateType coordinate) {
//                           Base_::Insert(coordinate_index, coordinate);
//                          }

typename WeightedVectorSpace::HomogeneousCoordinates_
WeightedVectorSpace::HomogenizeCoordinates(Coordinates_ const& coordinates,
                                           Weights_ const& weights,
                                           const int dimensionality) const {
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
  for (int i{}; i < number_of_coordinates; ++i) {

    const auto& w = weights[i];
    for (int j{}; j < dimensionality; ++j) {
      homogeneous_coordinates.emplace_back(*(current_coord++) * w);
    }
    homogeneous_coordinates.push_back(w);
  }
  return homogeneous_coordinates;
}

} // namespace bsplinelib::vector_spaces
