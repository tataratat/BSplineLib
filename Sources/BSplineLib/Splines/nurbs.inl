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

template<int para_dim>
Nurbs<para_dim>::Nurbs() : Base_(true) {}

template<int para_dim>
Nurbs<para_dim>::Nurbs(
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

template<int para_dim>
Nurbs<para_dim>::Nurbs(Nurbs const& other)
    : Base_(other),
      weighted_vector_space_{std::make_shared<WeightedVectorSpace_>(
          *other.weighted_vector_space_)},
      homogeneous_b_spline_{
          std::make_shared<HomogeneousBSpline_>(Base_::parameter_space_,
                                                weighted_vector_space_)} {}

template<int para_dim>
Nurbs<para_dim>& Nurbs<para_dim>::operator=(Nurbs const& rhs) {
  Base_::operator=(rhs);
  weighted_vector_space_ =
      std::make_shared<WeightedVectorSpace_>(*rhs.weighted_vector_space_);
  homogeneous_b_spline_ =
      std::make_shared<HomogeneousBSpline_>(Base_::parameter_space_,
                                            weighted_vector_space_);
  return *this;
}

template<int para_dim>
void Nurbs<para_dim>::Evaluate(const Type_* parametric_coordinate,
                               Type_* evaluated) const {
  const int h_dim = weighted_vector_space_->Dim();
  const int dim = h_dim - 1;

  Coordinate_ homogeneous_eval(h_dim);
  homogeneous_b_spline_->Evaluate(parametric_coordinate,
                                  homogeneous_eval.data());

  const Type_ w_inv = 1. / homogeneous_eval[dim];
  const Type_* eval_ptr = homogeneous_eval.begin();
  for (int i{}; i < dim; ++i) {
    *evaluated++ = *eval_ptr++ * w_inv;
  }
}

template<int para_dim>
void Nurbs<para_dim>::EvaluateDerivative(const Type_* parametric_coordinate,
                                         const IntType_* derivative,
                                         Type_* evaluated) const {

  using Data = bsplinelib::utilities::containers::Data<double>;
  using Data2D = bsplinelib::utilities::containers::Data<double, 2>;

  // Global (scalar) indexing to local index-system
  auto multi_index_ = [&derivative](int id) -> std::array<int, para_dim> {
    std::array<int, para_dim> local_ids{};
    for (int i{}; i < para_dim; ++i) {
      if (derivative[i] == 0) {
        // local_ids[i] = 0;
        continue;
      }
      local_ids[i] = id % (derivative[i] + 1);
      id -= local_ids[i];
      id /= derivative[i] + 1;
    }
    return local_ids;
  };

  // Check if requested derivative is "subset" to current derivative
  auto is_not_subset_ =
      [](const std::array<int, para_dim>& req_derivs_max,
         const std::array<int, para_dim>& req_derivs) -> bool {
    for (int i{}; i < para_dim; ++i) {
      if (req_derivs[i] > req_derivs_max[i])
        return true;
    }
    return false;
  };

  // Initialize return type
  const int number_of_derivs = [&derivative]() {
    int id{};
    int offset{1};
    for (int i{}; i < para_dim; ++i) {
      // assert(derivative[i] <= derivative[i]);
      if (derivative[i] == 0)
        continue;
      id += offset * (derivative[i]);
      offset = derivative[i] > 0 ? offset * (derivative[i] + 1) : offset;
    }
    return id;
  }() + 1;

  // Please remember that the first derivative is not used
  const int dim = Dim();
  // evaluated derivatives and a view array for row-wise operation
  Data2D der(number_of_derivs, dim);
  Data d_row;
  d_row.SetShape(dim);

  // evaluated derivatives from homogeneous bspline and a view array for
  // row-wise operation
  Data2D homogeneous_der(number_of_derivs, dim + 1);

  // Fill all homogeneous b spline derivatives (and values for id=0)
  for (int i{}; i < number_of_derivs; ++i) {
    const auto req_derivs = multi_index_(i);
    homogeneous_b_spline_->EvaluateDerivative(parametric_coordinate,
                                              req_derivs.data(),
                                              &homogeneous_der(i, 0));
  }

  // Precompute inverse of weighted function
  const double inv_w_fact = 1. / homogeneous_der(0, dim);

  // Loop over all lower-order derivatives and assign derivatives-vector
  // Notation follows "The NURBS book" eq. 4.20 (extended for n-d splines)
  for (int i{}; i < number_of_derivs; ++i) {
    // Retrieve index-wise order of the derivative for current ID
    const auto derivative_order_indexwise_LHS = multi_index_(i);
    // copy derivative of Numerator-function
    // d_row.Copy() copies according to d_row's size (dim).
    d_row.SetData(&der(i, 0));
    d_row.Copy(&homogeneous_der(i, 0));
    // Substract all weighted lower-order functions
    for (int j{1}; j <= i; ++j) {
      // Retrieve order of current index
      const auto derivative_order_indexwise_RHS = multi_index_(j);
      // Check only subsets
      if (is_not_subset_(derivative_order_indexwise_LHS,
                         derivative_order_indexwise_RHS))
        continue;
      // Precompute Product of binomial coefficients
      int binom_fact{1};
      for (int k{}; k < para_dim; ++k) {
        binom_fact *=
            bsplinelib::utilities::math_operations::ComputeBinomialCoefficient(
                derivative_order_indexwise_LHS[k],
                derivative_order_indexwise_RHS[k]);
      }
      // Substract low-order function
      d_row.Add(-(binom_fact * homogeneous_der(j, dim)), &der(i - j, 0));
    }
    // Finalize
    d_row.Multiply(inv_w_fact);
  }

  // Return last value
  std::copy_n(&der(number_of_derivs - 1, 0), dim, evaluated);
}

template<int para_dim>
typename Spline<para_dim>::Coordinate_
Nurbs<para_dim>::operator()(const Type_* parametric_coordinate) const {
  return WeightedVectorSpace_::Project(
      (*homogeneous_b_spline_)(parametric_coordinate));
}

// See NURBS book Eq. (4.20) (cf. Generalized Leibniz rule at
// <https://en.wikipedia.org/wiki/General_Leibniz_rule>).
template<int para_dim>
typename Spline<para_dim>::Coordinate_
Nurbs<para_dim>::operator()(const Type_* parametric_coordinate,
                            const IntType_* derivative) const {
  Coordinate_ evaluated_nurbs(Dim());

  EvaluateDerivative(parametric_coordinate, derivative, evaluated_nurbs.data());

  return evaluated_nurbs;
}

template<int para_dim>
void Nurbs<para_dim>::InsertKnot(Dimension const& dimension,
                                 Knot_ knot,
                                 Multiplicity const& multiplicity,
                                 Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::splines::Nurbs::InsertKnot"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
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

template<int para_dim>
Multiplicity Nurbs<para_dim>::RemoveKnot(Dimension const& dimension,
                                         Knot_ const& knot,
                                         Tolerance const& tolerance_removal,
                                         Multiplicity const& multiplicity,
                                         Tolerance const& tolerance) const {
  using std::get;

#ifndef NDEBUG
  using utilities::numeric_operations::ThrowIfToleranceIsNegative;

  Message const kName{"bsplinelib::splines::Nurbs::RemoveKnot"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
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

template<int para_dim>
void Nurbs<para_dim>::ElevateDegree(Dimension const& dimension,
                                    Multiplicity const& multiplicity,
                                    Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::splines::Nurbs::ElevateDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  homogeneous_b_spline_->ElevateDegree(dimension, multiplicity, tolerance);
}

template<int para_dim>
bool Nurbs<para_dim>::ReduceDegree(Dimension const& dimension,
                                   Tolerance const& tolerance_removal,
                                   Multiplicity const& multiplicity,
                                   Tolerance const& tolerance) const {
#ifndef NDEBUG
  using utilities::numeric_operations::ThrowIfToleranceIsNegative;

  Message const kName{"bsplinelib::splines::Nurbs::ReduceDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
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

template<int para_dim>
Coordinate
Nurbs<para_dim>::ComputeUpperBoundForMaximumDistanceFromOrigin() const {
  return std::get<0>(
      weighted_vector_space_
          ->DetermineMaximumDistanceFromOriginAndMinimumWeight());
}

template<int para_dim>
typename Nurbs<para_dim>::OutputInformation_
Nurbs<para_dim>::Write(Precision const& precision) const {
  return OutputInformation_{Base_::parameter_space_->Write(precision),
                            weighted_vector_space_->WriteProjected(precision)};
}
