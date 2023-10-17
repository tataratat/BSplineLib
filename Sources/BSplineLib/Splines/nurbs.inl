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

template<int parametric_dimensionality>
Nurbs<parametric_dimensionality>::Nurbs() : Base_(true) {}

template<int parametric_dimensionality>
Nurbs<parametric_dimensionality>::Nurbs(
    SharedPointer<ParameterSpace_> parameter_space,
    SharedPointer<WeightedVectorSpace_> weighted_vector_space)
    : Base_(std::move(parameter_space), true) {
  using std::to_string;

#ifndef NDEBUG
  int const &total_number_of_basis_functions =
                Base_::parameter_space_->GetTotalNumberOfBasisFunctions(),
            &number_of_coordinates =
                weighted_vector_space->GetNumberOfCoordinates();
  if (number_of_coordinates != total_number_of_basis_functions)
    Throw(DomainError(to_string(number_of_coordinates)
                      + " coordinates were provided but "
                      + to_string(total_number_of_basis_functions)
                      + " are needed to associate each basis function "
                        "with a weighted coordinate."),
          "bsplinelib::src::spl::Nurbs::Nurbs");
#endif
  homogeneous_b_spline_ =
      std::make_shared<HomogeneousBSpline_>(Base_::parameter_space_,
                                            weighted_vector_space);
  weighted_vector_space_ = std::move(weighted_vector_space);
}

template<int parametric_dimensionality>
Nurbs<parametric_dimensionality>::Nurbs(Nurbs const& other)
    : Base_(other),
      homogeneous_b_spline_{
          std::make_shared<HomogeneousBSpline_>(*other.homogeneous_b_spline_)} {
}

template<int parametric_dimensionality>
Nurbs<parametric_dimensionality>&
Nurbs<parametric_dimensionality>::operator=(Nurbs const& rhs) {
  Base_::operator=(rhs);
  homogeneous_b_spline_ =
      std::make_shared<HomogeneousBSpline_>(*rhs.homogeneous_b_spline_);
  return *this;
}

template<int parametric_dimensionality>
void Nurbs<parametric_dimensionality>::Evaluate(
    const DataType_* parametric_coordinate,
    DataType_* evaluated) const {
  // homogeneous dimension is one higher than physical dim
  const int h_dim = weighted_vector_space_->Dim();
  const int dim = h_dim - 1;

  TemporaryArray_<DataType_> homogeneous(h_dim);
  homogeneous_b_spline_->Evaluate(parametric_coordinate, homogeneous.data());

  // project
  const DataType_ w_inv = 1. / homogeneous[dim];
  for (int i{}; i < dim; ++i) {
    evaluated[i] = w_inv * homogeneous[i];
  }
}

template<int parametric_dimensionality>
typename Spline<parametric_dimensionality>::Coordinate_
Nurbs<parametric_dimensionality>::operator()(
    ParametricCoordinate_ const& parametric_coordinate) const {
  Coordinate_ evaluated(weighted_vector_space_->Dim() - 1);
  Evaluate(parametric_coordinate.data(), evaluated.data());
  return evaluated;
}

// Adapted from bezman
template<int parametric_dimensionality>
void Nurbs<parametric_dimensionality>::EvaluateDerivative(
    const DataType_* parametric_coordinate,
    const IndexType_* derivative,
    DataType_* evaluated) const {

  // Define lambdas to switch between local indexing with coordinate style
  // indexing to global, scalar indexing
  // example: for nth (2,1,4) : (1,0,0)->1, (1,1,0)->4

  // Global (scalar) indexing to local index-system
  auto local_ids_ =
      [&derivative](
          const int& req_id) -> std::array<int, parametric_dimensionality> {
    int id{req_id};
    std::array<int, parametric_dimensionality> local_ids{};
    for (int i_pd{}; i_pd < parametric_dimensionality; i_pd++) {
      if (derivative[i_pd] == 0)
        continue;
      local_ids[i_pd] = id % (derivative[i_pd] + 1);
      id -= local_ids[i_pd];
      id /= derivative[i_pd] + 1;
    }
    return local_ids;
  };

  // Local (coordinate-style) indexing to global
  auto global_ids_ =
      [&derivative](
          const DataType_* req_derivs) -> int {
    int id{};
    int offset{1};
    for (int i_pd{}; i_pd < parametric_dimensionality; i_pd++) {
      // assert(req_derivs[i_pd] <= derivative[i_pd]);
      if (derivative[i_pd] == 0)
        continue;
      id += offset * (req_derivs[i_pd]);
      offset = derivative[i_pd] > 0 ? offset * (derivative[i_pd] + 1) : offset;
    }
    return id;
  };

  // Check if requested derivative is "subset" to current derivative
  auto is_not_subset_ =
      [](const std::array<int, parametric_dimensionality>& req_derivs_max,
         const std::array<int, parametric_dimensionality>& req_derivs) -> bool {
    for (int i_pd{}; i_pd < parametric_dimensionality; i_pd++) {
      if (req_derivs[i_pd] > req_derivs_max[i_pd])
        return true;
    }
    return false;
  };

  // Initialize return type
  const int number_of_derivs{global_ids_(derivative) + 1};
  const int h_dim = weighted_vector_space_->Dim();
  const int dim{h_dim - 1};
  // Please remember that the first derivative is not used
  TemporaryArray2D_<DataType_> derivatives(number_of_derivs, dim);
  TemporaryArray2D_<DataType_> A_derivatives(number_of_derivs, dim);
  TemporaryArray<DataType_> w_derivatives(number_of_derivs);
  // tmp array to split A_der and w_der
  TemporaryArray<DataType_> homogeneous(h_dim);

  auto assgin_der_A_w = [&derivatives, &A_derivatives, &w_derivatives, &homogeneous, &dim] (const int& i_deriv) {
    DataType* der = &derivatives(i_deriv, 0);
    DataType* A = &A_derivatives(i_deriv, 0);
    DataType* h = homogeneous.data_;
    for (int i{}; i < dim; ++i) {
      const DataType_& val = *h;
      *der++ = val;
      *A++ = val;
      ++h;
    }
    w_derivatives[i_deriv] = *h;
  }

  // Fill all polynomial spline derivatives (and values for id=0)
  for (int i_deriv{}; i_deriv < number_of_derivs; i_deriv++) {
    const auto req_derivs = local_ids_(i_deriv);
    // evaluate homogeneous
    homogeneous_b_spline_->EvaluateDerivative(parametric_coordinate, req_derivs.data(), homogeneous.data_);

    assign_der_A_w(i_deriv);
  }

  // Precompute inverse of weighted function
  const DataType_ inv_w_fact = 1. / w_derivatives[0];

  // Loop over all lower-order derivatives and assign derivatives-vector
  // Notation follows "The NURBS book" eq. 4.20 (extended for n-d splines)
  for (int i_deriv{}; i_deriv < number_of_derivs; i_deriv++) {
    // Retrieve index-wise order of the derivative for current ID
    const auto derivative_order_indexwise_LHS = local_ids_(i_deriv);
    // Assign derivative of Numerator-function
    derivatives[i_deriv] = A_derivatives[i_deriv];
    // Substract all weighted lower-order functions
    for (int j_deriv{1}; j_deriv <= i_deriv; j_deriv++) {
      // Retrieve order of current index
      const auto derivative_order_indexwise_RHS = local_ids_(j_deriv);
      // Check only subsets
      if (is_not_subset_(derivative_order_indexwise_LHS,
                         derivative_order_indexwise_RHS))
        continue;
      // Precompute Product of binomial coefficients
      int binom_fact{1};
      for (int i_pd{}; i_pd < parametric_dimensionality; i_pd++) {
        binom_fact *=
            bsplinelib::utilities::math_operations::ComputeBinomialCoefficient(
                derivative_order_indexwise_LHS[i_pd],
                derivative_order_indexwise_RHS[i_pd]);
      }
      // Substract low-order function
      derivatives[i_deriv] -=
          binom_fact * w_derivatives[j_deriv] * derivatives[i_deriv - j_deriv];
    }
    // Finalize
    derivatives[i_deriv] *= inv_w_fact;
  }
  // Return last value
  return derivatives[number_of_derivs - 1];
}

template<int parametric_dimensionality>
void Nurbs<parametric_dimensionality>::InsertKnot(
    Dimension const& dimension,
    Knot_ knot,
    Multiplicity const& multiplicity,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::splines::Nurbs::InsertKnot"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension,
                                                parametric_dimensionality - 1);
    Base_::parameter_space_->ThrowIfParametricCoordinateIsOutsideScope(
        dimension,
        knot,
        tolerance);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  }
#endif
  homogeneous_b_spline_->InsertKnot(dimension, knot, multiplicity, tolerance);
}

template<int parametric_dimensionality>
Multiplicity
Nurbs<parametric_dimensionality>::RemoveKnot(Dimension const& dimension,
                                             Knot_ const& knot,
                                             Tolerance const& tolerance_removal,
                                             Multiplicity const& multiplicity,
                                             Tolerance const& tolerance) const {
  using std::get;

#ifndef NDEBUG
  using utilities::numeric_operations::ThrowIfToleranceIsNegative;

  Message const kName{"bsplinelib::splines::Nurbs::RemoveKnot"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension,
                                                parametric_dimensionality - 1);
    ThrowIfToleranceIsNegative(tolerance_removal);
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  auto const& [maximum_distance_from_origin, minimum_weight] =
      weighted_vector_space_
          ->DetermineMaximumDistanceFromOriginAndMinimumWeight();
  // See NURBS book Eq. (5.30).
  return homogeneous_b_spline_->RemoveKnot(
      dimension,
      knot,
      tolerance_removal
          * (minimum_weight / (1.0 + maximum_distance_from_origin)),
      multiplicity,
      tolerance);
}

template<int parametric_dimensionality>
void Nurbs<parametric_dimensionality>::ElevateDegree(
    Dimension const& dimension,
    Multiplicity const& multiplicity,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::splines::Nurbs::ElevateDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension,
                                                parametric_dimensionality - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  homogeneous_b_spline_->ElevateDegree(dimension, multiplicity, tolerance);
}

template<int parametric_dimensionality>
bool Nurbs<parametric_dimensionality>::ReduceDegree(
    Dimension const& dimension,
    Tolerance const& tolerance_removal,
    Multiplicity const& multiplicity,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  using utilities::numeric_operations::ThrowIfToleranceIsNegative;

  Message const kName{"bsplinelib::splines::Nurbs::ReduceDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension,
                                                parametric_dimensionality - 1);
    ThrowIfToleranceIsNegative(tolerance_removal);
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  return homogeneous_b_spline_->ReduceDegree(dimension,
                                             tolerance_removal,
                                             multiplicity,
                                             tolerance);
}

template<int parametric_dimensionality>
Coordinate
Nurbs<parametric_dimensionality>::ComputeUpperBoundForMaximumDistanceFromOrigin(
    Tolerance const& tolerance) const {
  return std::get<0>(
      weighted_vector_space_
          ->DetermineMaximumDistanceFromOriginAndMinimumWeight(tolerance));
}
