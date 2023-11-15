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
ParameterSpace<para_dim>::ParameterSpace(ParameterSpace const& other)
    : degrees_{other.degrees_} {
  CopyKnotVectors(other.knot_vectors_);
}

template<int para_dim>
ParameterSpace<para_dim>&
ParameterSpace<para_dim>::operator=(ParameterSpace const& rhs) {
  CopyKnotVectors(rhs.knot_vectors_);
  degrees_ = rhs.degrees_;
  return *this;
}

template<int para_dim>
typename ParameterSpace<para_dim>::Index_
ParameterSpace<para_dim>::First() const {
  return Index_::First(GetNumberOfNonZeroBasisFunctions());
}

template<int para_dim>
typename ParameterSpace<para_dim>::Index_
ParameterSpace<para_dim>::Behind() const {
  return Index_::Behind(GetNumberOfNonZeroBasisFunctions());
}

template<int para_dim>
typename ParameterSpace<para_dim>::NumberOfBasisFunctions_
ParameterSpace<para_dim>::GetNumberOfBasisFunctions() const {
  NumberOfBasisFunctions_ number_of_basis_functions;
  Dimension::ForEach(0, para_dim, [&](Dimension const& dimension) {
    number_of_basis_functions[dimension.Get()] =
        GetNumberOfBasisFunctions(dimension);
  });
  return number_of_basis_functions;
}

template<int para_dim>
int ParameterSpace<para_dim>::GetTotalNumberOfBasisFunctions() const {
  NumberOfBasisFunctions_ const& number_of_basis_functions =
      GetNumberOfBasisFunctions();
  return std::reduce(number_of_basis_functions.begin(),
                     number_of_basis_functions.end(),
                     Length{1},
                     std::multiplies{})
      .Get();
}

template<int para_dim>
SharedPointer<ParameterSpace<para_dim - 1>>
ParameterSpace<para_dim>::RemoveOneParametricDimension(
    const int parametric_dimension) const {
  if constexpr (para_dim < 2) {
    return nullptr;
  } else {
    // prepare output
    auto out_pspace = std::make_shared<ParameterSpace<para_dim - 1>>();

    auto& out_degrees = out_pspace->GetDegrees();
    auto& out_knot_vectors = out_pspace->GetKnotVectors();

    int out_dim_id{};
    for (int i{}; i < para_dim; ++i) {
      if (i == parametric_dimension) {
        continue;
      }
      out_degrees[out_dim_id] = degrees_[i];
      out_knot_vectors[out_dim_id] = knot_vectors_[i];
      ++out_dim_id;
    }

    return out_pspace;
  }
}

template<int para_dim>
typename ParameterSpace<para_dim>::Index_
ParameterSpace<para_dim>::FindFirstNonZeroBasisFunction(
    const Type_* parametric_coordinate,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception,
          "bsplinelib::parameter_spaces::ParameterSpace::"
          "FindFirstNonZeroBasisFunction");
  }
#endif
  IndexValue_ first_non_zero_basis_function_index_value;
  Dimension::ForEach(0, para_dim, [&](Dimension const& dimension) {
    Dimension::Type_ const& current_dimension = dimension.Get();
    first_non_zero_basis_function_index_value[current_dimension] =
        Index{knot_vectors_[current_dimension]->FindSpan(
                  parametric_coordinate[current_dimension],
                  tolerance)
              - degrees_[current_dimension]};
  });
  return Index_{GetNumberOfBasisFunctions(),
                first_non_zero_basis_function_index_value};
}
template<int para_dim>
typename ParameterSpace<para_dim>::KnotSpans_
ParameterSpace<para_dim>::FindKnotSpans(const Type_* parametric_coordinate,
                                        Tolerance const& tolerance) const {
  KnotSpans_ out;
  for (int i{}; i < para_dim; ++i) {
    out[i] = knot_vectors_[i]->FindSpan(parametric_coordinate[i], tolerance);
  }
  return out;
}

template<int para_dim>
typename ParameterSpace<para_dim>::BezierInformation_
ParameterSpace<para_dim>::DetermineBezierExtractionKnots(
    Dimension const& dimension,
    Tolerance const& tolerance) const {
  Dimension::Type_ const& dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::ParameterSpace::"
                      "DetermineBezierExtractionKnots"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName, dimension_value);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName, dimension_value);
  }
#endif
  KnotVector const& knot_vector = *knot_vectors_[dimension_value];
  Knots_ const& unique_knots = knot_vector.GetUniqueKnots(tolerance);
  typename Knots_::const_iterator const &first_interior_knot =
                                            (unique_knots.begin() + 1),
                                        &last_knot =
                                            std::prev(unique_knots.end());
  int const& number_of_interior_knots = (last_knot - first_interior_knot);
  Knots_ bezier_extraction_knots;
  bezier_extraction_knots.reserve(number_of_interior_knots);
  std::for_each(
      first_interior_knot,
      last_knot,
      [&](ParametricCoordinate const& interior_knot) {
        std::fill_n(
            std::back_inserter(bezier_extraction_knots),
            degrees_[dimension_value]
                - knot_vector.DetermineMultiplicity(interior_knot).Get(),
            interior_knot);
      });
  return BezierInformation_{number_of_interior_knots + 1,
                            bezier_extraction_knots};
}

template<int para_dim>
Vector<Vector<int>>
ParameterSpace<para_dim>::KnotMultiplicities(Tolerance const& tolerance) const {

  Vector<Vector<int>> multiplicities;
  multiplicities.reserve(para_dim);
  for (int i{}; i < para_dim; ++i) {
    multiplicities.emplace_back(
        knot_vectors_[i]->DetermineMultiplicities(tolerance));
  }

  return multiplicities;
}

template<int para_dim>
typename ParameterSpace<para_dim>::BasisValuesPerDimension_
ParameterSpace<para_dim>::EvaluateBasisValuesPerDimension(
    const Type_* parametric_coordinate,
    Tolerance const& tolerance) const {

  // prepare output
  BasisValuesPerDimension_ output;
  for (int i{}; i < para_dim; ++i) {

    // get this dim's info
    const auto& this_dim_degree = degrees_[i];
    const auto this_dim_n_basis = this_dim_degree + 1;
    const auto& this_dim_parametric_coordinate = parametric_coordinate[i];
    const auto& this_knot_vector = *knot_vectors_[i];
    const auto& this_knots = this_knot_vector.GetKnots();
    const int this_zero_degree_support =
        this_knot_vector.FindSpan(this_dim_parametric_coordinate, tolerance)
            .Get();

    // this dim's output
    auto& this_dim_output = output[i];
    this_dim_output.Reallocate(this_dim_n_basis);
    this_dim_output[0] = 1.;

    TemporaryData_<double> left(this_dim_n_basis), right(this_dim_n_basis);

    double saved, temp;

    for (int k{1}; k < this_dim_n_basis; ++k) {
      this_dim_output[k] = 1.;

      left[k] = this_dim_parametric_coordinate
                - this_knots[this_zero_degree_support + 1 - k];
      right[k] = this_knots[this_zero_degree_support + k]
                 - this_dim_parametric_coordinate;

      saved = 0.0;
      for (int j{}; j < k; j++) {
        temp = this_dim_output[j] / (right[j + 1] + left[k - j]);
        this_dim_output[j] = saved + right[j + 1] * temp;
        saved = left[k - j] * temp;
      }
      this_dim_output[k] = saved;
    }
  }

  return output;
}

template<int para_dim>
typename ParameterSpace<para_dim>::BasisValues_
ParameterSpace<para_dim>::EvaluateBasisValues(
    const Type_* parametric_coordinate,
    Tolerance const& tolerance) const {
  return RecursiveCombine(
      EvaluateBasisValuesPerDimension(parametric_coordinate, tolerance));
}

template<int para_dim>
typename ParameterSpace<para_dim>::BasisValuesPerDimension_
ParameterSpace<para_dim>::EvaluateBasisDerivativeValuesPerDimension(
    const Type_* parametric_coordinate,
    const IntType_* derivative,
    Tolerance const& tolerance) const {

  // prepare output
  BasisValuesPerDimension_ output;

  for (int i{}; i < para_dim; ++i) {

    // get this dim's info
    const auto& this_dim_degree = degrees_[i];
    const auto& this_dim_derivative = derivative[i];
    const auto this_dim_n_basis = this_dim_degree + 1;
    // this dim's output and allocate
    auto& this_dim_output = output[i];
    this_dim_output.Reallocate(this_dim_n_basis);

    // special case for early exit - derivetive query is bigger than degree
    // all zeros.
    if (this_dim_derivative > this_dim_degree) {
      this_dim_output.Fill(0.);
      continue;
    }

    const auto& this_dim_parametric_coordinate = parametric_coordinate[i];
    const auto& this_knot_vector = *knot_vectors_[i];
    const auto& this_knots = this_knot_vector.GetKnots();
    const int this_zero_degree_support =
        this_knot_vector.FindSpan(this_dim_parametric_coordinate, tolerance)
            .Get();

    // temporary ones that we need for second special case
    TemporaryData_<double> left(this_dim_n_basis), right(this_dim_n_basis);
    double saved, temp, d;

    // special case 2 - derivative 0 query is evaluation query
    if (this_dim_derivative == 0) {
      // just raw copy of evaluation
      this_dim_output[0] = 1.;

      for (int k{1}; k < this_dim_n_basis; ++k) {
        this_dim_output[k] = 1.;

        left[k] = this_dim_parametric_coordinate
                  - this_knots[this_zero_degree_support + 1 - k];
        right[k] = this_knots[this_zero_degree_support + k]
                   - this_dim_parametric_coordinate;

        saved = 0.0;
        for (int j{}; j < k; j++) {
          temp = this_dim_output[j] / (right[j + 1] + left[k - j]);
          this_dim_output[j] = saved + right[j + 1] * temp;
          saved = left[k - j] * temp;
        }
        this_dim_output[k] = saved;
      }
      continue;
    }

    // here, proper derivative query.
    // more temporary variables
    TemporaryData2D_<double> a(2, this_dim_n_basis),
        ndu(this_dim_n_basis, this_dim_n_basis);
    int j1, j2;

    ndu(0, 0) = 1.;
    for (int j{1}; j < this_dim_n_basis; ++j) {
      left[j] = this_dim_parametric_coordinate
                - this_knots[this_zero_degree_support + 1 - j];
      right[j] = this_knots[this_zero_degree_support + j]
                 - this_dim_parametric_coordinate;

      saved = 0.0;
      for (int r{}; r < j; ++r) {
        ndu(j, r) = right[r + 1] + left[j - r];
        temp = ndu(r, j - 1) / ndu(j, r);
        ndu(r, j) = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      ndu(j, j) = saved;
    }

    if (this_dim_derivative == 0) {
      for (int j{}; j < this_dim_n_basis; ++j) {
        this_dim_output[j] = ndu(j, this_dim_degree);
      }
      continue;
    }

    for (int r{}; r < this_dim_n_basis; ++r) {
      int s1{}, s2{1}, j;
      a(0, 0) = 1.0;
      for (int k{1}; k < this_dim_derivative + 1; ++k) {
        d = 0.0;
        const int rk = r - k;
        const int pk = this_dim_degree - k;
        if (r >= k) {
          a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
          d = a(s2, 0) * ndu(rk, pk);
        }

        if (rk >= -1) {
          j1 = 1;
        } else {
          j1 = -rk;
        }

        if (r - 1 <= pk) {
          j2 = k - 1;
        } else {
          j2 = this_dim_degree - r;
        }

        for (j = j1; j < j2 + 1; ++j) {
          a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
          d += a(s2, j) * ndu(rk + j, pk);
        }

        if (r <= pk) {
          a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
          d += a(s2, j) * ndu(rk + j, pk);
        }
        // here, we could also gather lower derivative queries
        // for now, just rewrite at the same place
        this_dim_output[r] = d;

        j = s1;
        s1 = s2;
        s2 = j;
      }
    }

    temp = this_dim_degree;
    for (int k{1}; k < this_dim_derivative; ++k) {
      temp *= (this_dim_degree - k);
    }
    for (int j{}; j < this_dim_n_basis; ++j) {
      this_dim_output[j] *= temp;
    }
  }

  return output;
}

template<int para_dim>
typename ParameterSpace<para_dim>::BasisValues_
ParameterSpace<para_dim>::EvaluateBasisDerivativeValues(
    const Type_* parametric_coordinate,
    const IntType_* derivative,
    Tolerance const& tolerance) const {
  return RecursiveCombine(
      EvaluateBasisDerivativeValuesPerDimension(parametric_coordinate,
                                                derivative,
                                                tolerance));
}

template<int para_dim>
typename ParameterSpace<para_dim>::InsertionInformation_
ParameterSpace<para_dim>::InsertKnot(Dimension const& dimension,
                                     Knot_ knot,
                                     Multiplicity const& multiplicity,
                                     Tolerance const& tolerance) {
  Dimension::Type_ const& dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName(
      "bsplinelib::parameter_spaces::ParameterSpace::InsertKnot");

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
    knot_vectors_[dimension_value]->ThrowIfParametricCoordinateIsOutsideScope(
        knot,
        tolerance);
    ThrowIfFrontOrBackKnotIsToBeInsertedOrRemoved(dimension, knot, tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName, dimension_value);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName, dimension_value);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName, dimension_value);
  }
#endif
  KnotVector& knot_vector = *knot_vectors_[dimension_value];
  // instead of (degree - current_multiplicity),
  // clip at (degree - current_multiplicity) + 1 is to allow C^-1
  Multiplicity const insertion{std::min(
      multiplicity.Get(),
      degrees_[dimension_value]
          - knot_vector.DetermineMultiplicity(knot, tolerance).Get() + 1)};
  InsertionInformation_ const& insertion_information =
      DetermineInsertionInformation(dimension, knot, insertion, tolerance);
  knot_vector.Insert(std::move(knot), insertion, tolerance);
  return insertion_information;
}

template<int para_dim>
typename ParameterSpace<para_dim>::InsertionInformation_
ParameterSpace<para_dim>::RemoveKnot(Dimension const& dimension,
                                     Knot_ const& knot,
                                     Multiplicity const& multiplicity,
                                     Tolerance const& tolerance) {
  Dimension::Type_ const& dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName{
      "bsplinelib::parameter_spaces::ParameterSpace::RemoveKnot"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
    ThrowIfFrontOrBackKnotIsToBeInsertedOrRemoved(dimension, knot, tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName, dimension_value);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName, dimension_value);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName, dimension_value);
  }
#endif
  if (Multiplicity const& removals =
          knot_vectors_[dimension_value]->Remove(knot, multiplicity, tolerance);
      removals != Multiplicity{}) {
    return DetermineInsertionInformation(dimension, knot, removals, tolerance);
  } else {
    return InsertionInformation_{};
  }
}

// Cf. NURBS book Eq. (5.36).
template<int para_dim>
typename ParameterSpace<para_dim>::ElevationInformation_
ParameterSpace<para_dim>::ElevateDegree(Dimension const& dimension,
                                        Multiplicity const& multiplicity,
                                        Tolerance const& tolerance) {
  Dimension::Type_ const& dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName{
      "bsplinelib::parameter_spaces::ParameterSpace::ElevateDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName, dimension_value);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName, dimension_value);
  }
#endif
  ElevationInformation_ const& bezier_information =
      DetermineElevationInformation(dimension, multiplicity);
  knot_vectors_[dimension_value]->IncreaseMultiplicities(multiplicity);
  degrees_[dimension_value] += Degree{multiplicity.Get()};
  return bezier_information;
}

template<int para_dim>
typename ParameterSpace<para_dim>::ElevationInformation_
ParameterSpace<para_dim>::ReduceDegree(Dimension const& dimension,
                                       Multiplicity const& multiplicity,
                                       Tolerance const& tolerance) {
  Dimension::Type_ const& dimension_value = dimension.Get();
#ifndef NDEBUG
  Message const kName{
      "bsplinelib::parameter_spaces::ParameterSpace::ReduceDegree"};

  try {
    Dimension::ThrowIfNamedIntegerIsOutOfBounds(dimension, para_dim - 1);
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName, dimension_value);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName, dimension_value);
  }
#endif
  Degree& degree = degrees_[dimension_value];
  Multiplicity const reduction{std::min(multiplicity.Get(), degree - 1)};
  knot_vectors_[dimension_value]->DecreaseMultiplicities(reduction);
  degree -= reduction.Get();
  return DetermineElevationInformation(dimension, reduction);
}

template<int para_dim>
typename ParameterSpace<para_dim>::OutputInformation_
ParameterSpace<para_dim>::Write(Precision const& precision) const {
  using std::get, std::tuple_element_t, utilities::string_operations::Write;

  OutputInformation_ knot_vectors_degrees_and_maximum_basis_function_indices;
  Dimension::ForEach(0, para_dim, [&](Dimension const& dimension) {
    Dimension::Type_ const& current_dimension = dimension.Get();
    get<0>(knot_vectors_degrees_and_maximum_basis_function_indices)
        [current_dimension] =
            knot_vectors_[current_dimension]->Write(precision);
  });
  get<1>(knot_vectors_degrees_and_maximum_basis_function_indices) =
      Write<tuple_element_t<1, OutputInformation_>>(degrees_, precision);
  get<2>(knot_vectors_degrees_and_maximum_basis_function_indices) =
      Write<tuple_element_t<2, OutputInformation_>>(GetNumberOfBasisFunctions(),
                                                    precision);
  return knot_vectors_degrees_and_maximum_basis_function_indices;
}

#ifndef NDEBUG
template<int para_dim>
void ParameterSpace<para_dim>::ThrowIfParametricCoordinateIsOutsideScope(
    Dimension const& dimension,
    ParametricCoordinate const& parametric_coordinate,
    Tolerance const& tolerance) const {
  knot_vectors_[dimension.Get()]->ThrowIfParametricCoordinateIsOutsideScope(
      parametric_coordinate,
      tolerance);
}
#endif

template<int para_dim>
void ParameterSpace<para_dim>::CopyKnotVectors(
    KnotVectors_ const& knot_vectors) {
  Dimension::ForEach(0, para_dim, [&](Dimension const& dimension) {
    Dimension::Type_ const& current_dimension = dimension.Get();
    knot_vectors_[current_dimension] =
        std::make_shared<KnotVector>(*knot_vectors[current_dimension]);
  });
}

template<int para_dim>
int ParameterSpace<para_dim>::GetNumberOfNonZeroBasisFunctions(
    Dimension const& dimension) const {
  return degrees_[dimension.Get()] + 1;
}

template<int para_dim>
typename ParameterSpace<para_dim>::NumberOfBasisFunctions_
ParameterSpace<para_dim>::GetNumberOfNonZeroBasisFunctions() const {
  NumberOfBasisFunctions_ number_of_non_zero_basis_functions;
  Dimension::ForEach(0, para_dim, [&](Dimension const& dimension) {
    number_of_non_zero_basis_functions[dimension.Get()] =
        Length{GetNumberOfNonZeroBasisFunctions(dimension)};
  });
  return number_of_non_zero_basis_functions;
}

template<int para_dim>
Length ParameterSpace<para_dim>::GetNumberOfBasisFunctions(
    Dimension const& dimension) const {
  return Length{knot_vectors_[dimension.Get()]->GetSize()
                - GetNumberOfNonZeroBasisFunctions(dimension)};
}

// Cf. NURBS book below Eq. (5.15).
template<int para_dim>
typename ParameterSpace<para_dim>::InsertionInformation_
ParameterSpace<para_dim>::DetermineInsertionInformation(
    Dimension const& dimension,
    ParametricCoordinate const& knot,
    Multiplicity const& multiplicity,
    Tolerance const& tolerance) const {
  Dimension::Type_ const& dimension_value = dimension.Get();
  KnotVector& knot_vector = *knot_vectors_[dimension_value];
  MultiplicityType_ const& multiplicity_value = multiplicity.Get();
  MultiplicityType_ const previous_multiplicity{
      knot_vector.DetermineMultiplicity(knot, tolerance).Get()};
  KnotSpan::Type_ const knot_span{knot_vector.FindSpan(knot, tolerance).Get()};
  InsertionCoefficients_ knot_ratios;
  knot_ratios.reserve(multiplicity_value);
  Index const end_knot{knot_span - previous_multiplicity};
  Multiplicity::ForEach(
      0,
      multiplicity_value,
      [&](Multiplicity const& insertion) {
        KnotRatios_ current_knot_ratios;
        MultiplicityType_ const remaining_mulitplicity{degrees_[dimension_value]
                                                       - insertion.Get()};
        current_knot_ratios.reserve(remaining_mulitplicity
                                    - previous_multiplicity);
        Index::ForEach(
            knot_span - remaining_mulitplicity + 1,
            end_knot.Get() + 1,
            [&](Index const& current_knot) {
              ParametricCoordinate const& lower_knot =
                  knot_vector[current_knot];
              current_knot_ratios.push_back(
                  (knot - lower_knot)
                  / (knot_vector[current_knot + Index{remaining_mulitplicity}]
                     - lower_knot));
            });
        knot_ratios.push_back(current_knot_ratios);
      });

  return InsertionInformation_{end_knot, knot_ratios};
}

template<int para_dim>
typename ParameterSpace<para_dim>::ElevationInformation_
ParameterSpace<para_dim>::DetermineElevationInformation(
    Dimension const& dimension,
    Multiplicity const& multiplicity) const {
  using IndexType = Index::Type_;
  using utilities::math_operations::ComputeBinomialCoefficient;

  Degree const& degree = degrees_[dimension.Get()];
  MultiplicityType_ const& multiplicity_value = multiplicity.Get();
  Index::Type_ const& maximum_bezier_coordinate = (degree + multiplicity_value);
  ElevationCoefficients_ bezier_coefficients;
  bezier_coefficients.reserve(maximum_bezier_coordinate - 1);
  Index::ForEach(
      1,
      maximum_bezier_coordinate,
      [&](Index const& bezier_coordinate) {
        IndexType const &bezier_coordinate_value = bezier_coordinate.Get(),
                        &current_bezier_coordinate_end{
                            std::min(degree, bezier_coordinate_value) + 1};
        // Capturing result by reference produces dangling references if maxima
        // of the parameters are temporaries.
        IndexType const current_bezier_coordinate_begin{
            std::max(0, bezier_coordinate_value - multiplicity_value)};
        BinomialRatios_ current_bezier_coefficients(
            current_bezier_coordinate_end - current_bezier_coordinate_begin,
            1.0
                / ComputeBinomialCoefficient(maximum_bezier_coordinate,
                                             bezier_coordinate_value));
        Index::ForEach(
            current_bezier_coordinate_begin,
            current_bezier_coordinate_end,
            [&](Index const& current_bezier_coordinate) {
              IndexType const& current_bezier_coordinate_value =
                  current_bezier_coordinate.Get();
              current_bezier_coefficients[current_bezier_coordinate_value
                                          - current_bezier_coordinate_begin] *=
                  (ComputeBinomialCoefficient(degree,
                                              current_bezier_coordinate_value)
                   * ComputeBinomialCoefficient(
                       multiplicity_value,
                       bezier_coordinate_value
                           - current_bezier_coordinate_value));
            });
        bezier_coefficients.push_back(current_bezier_coefficients);
      });
  return ElevationInformation_{Index{degree}, bezier_coefficients};
}

#ifndef NDEBUG
template<int para_dim>
void ParameterSpace<para_dim>::ThrowIfBasisFunctionIndexIsInvalid(
    Index_ const& basis_function_index) const {
  Dimension::ForEach(0, para_dim, [&](Dimension const& dimension) {
    try {
      Index::ThrowIfNamedIntegerIsOutOfBounds(
          basis_function_index[dimension],
          GetNumberOfBasisFunctions(dimension).Get() - 1);
    } catch (OutOfRange const& exception) {
      throw OutOfRange("for dimension " + std::to_string(dimension.Get()) + ": "
                       + exception.what());
    }
  });
}

template<int para_dim>
void ParameterSpace<para_dim>::ThrowIfFrontOrBackKnotIsToBeInsertedOrRemoved(
    Dimension const& dimension,
    ParametricCoordinate const& knot,
    Tolerance const& tolerance) const {
  using std::to_string;

  Dimension::Type_ const& dimension_value = dimension.Get();
  if (knot_vectors_[dimension_value]->DoesParametricCoordinateEqualFrontOrBack(
          knot,
          tolerance))
    throw DomainError("for dimension" + to_string(dimension_value)
                      + ": Cannot insert or remove the first or last "
                        "knot "
                      + to_string(knot)
                      + " as only clamped knot vectors are allowed.");
}
#endif
