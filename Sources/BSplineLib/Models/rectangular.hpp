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

#ifndef SOURCES_MODELS_RECTANGULAR_HPP_
#define SOURCES_MODELS_RECTANGULAR_HPP_

#include "BSplineLib/ParameterSpaces/parameter_space.hpp"
#include "BSplineLib/Splines/b_spline.hpp"
#include "BSplineLib/Splines/nurbs.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/std_container_operations.hpp"

// Generates maximally-smooth spline models of rectangular unit domains.
//
// Example:
//   Length const kNumberOfNonZeroKnotSpans{10};
//   Nurbs<2, 2> rectangle = Generate<2, 2, true>({kNumberOfNonZeroKnotSpans,
//   kNumberOfNonZeroKnotSpans});
namespace bsplinelib::models::rectangular {

template<int parametric_dimensionality>
using Degrees = typename parameter_spaces::ParameterSpace<
    parametric_dimensionality>::Degrees_;
template<int parametric_dimensionality>
using NumbersOfNonZeroKnotSpans = Array<Length, parametric_dimensionality>;

template<int parametric_dimensionality,
         int geometric_dimensionality,
         bool is_rational>
auto Generate(NumbersOfNonZeroKnotSpans<parametric_dimensionality> const&
                  numbers_of_non_zero_knot_spans,
              Degrees<parametric_dimensionality> const& degrees);

#include "BSplineLib/Models/rectangular.inc"

} // namespace bsplinelib::models::rectangular

#endif // SOURCES_MODELS_RECTANGULAR_HPP_
