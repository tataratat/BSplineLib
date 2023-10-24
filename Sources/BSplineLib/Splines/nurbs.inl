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
    const Type_* parametric_coordinate,
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

template<int parametric_dimensionality>
void Nurbs<parametric_dimensionality>::EvaluateDerivative(
    const Type_* parametric_coordinate,
    const IntType_* derivative,
    Type_* evaluated) const {}

template<int parametric_dimensionality>
typename Spline<parametric_dimensionality>::Coordinate_
Nurbs<parametric_dimensionality>::operator()(
    const Type_* parametric_coordinate) const {
  return WeightedVectorSpace_::Project(
      (*homogeneous_b_spline_)(parametric_coordinate));
}

// See NURBS book Eq. (4.20) (cf. Generalized Leibniz rule at
// <https://en.wikipedia.org/wiki/General_Leibniz_rule>).
template<int parametric_dimensionality>
typename Spline<parametric_dimensionality>::Coordinate_
Nurbs<parametric_dimensionality>::operator()(const Type_* parametric_coordinate,
                                             const IntType_* derivative) const {
  Coordinate_ evaluated_nurbs(Dim());

  EvaluateDerivative(parametric_coordinate, derivative, evaluated_nurbs.data());

  return evaluated_nurbs;
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
Coordinate Nurbs<parametric_dimensionality>::
    ComputeUpperBoundForMaximumDistanceFromOrigin() const {
  return std::get<0>(
      weighted_vector_space_
          ->DetermineMaximumDistanceFromOriginAndMinimumWeight());
}

template<int parametric_dimensionality>
typename Nurbs<parametric_dimensionality>::OutputInformation_
Nurbs<parametric_dimensionality>::Write(Precision const& precision) const {
  return OutputInformation_{Base_::parameter_space_->Write(precision),
                            weighted_vector_space_->WriteProjected(precision)};
}
