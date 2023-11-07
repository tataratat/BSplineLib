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
void Spline<para_dim>::RefineKnots(Dimension const& dimension,
                                   Knots_ knots,
                                   Multiplicity const& multiplicity,
                                   Tolerance const& tolerance) const {
  std::for_each(knots.begin(), knots.end(), [&](Knot_ const& knot) {
    InsertKnot(dimension, std::move(knot), multiplicity, tolerance);
  });
}

template<int para_dim>
Multiplicity Spline<para_dim>::CoarsenKnots(Dimension const& dimension,
                                            Knots_ const& knots,
                                            Tolerance const& tolerance_removal,
                                            Multiplicity const& multiplicity,
                                            Tolerance const& tolerance) const {
  Multiplicity successful_removals{multiplicity};
  std::for_each(knots.begin(), knots.end(), [&](Knot_ const& knot) {
    successful_removals = std::min(successful_removals,
                                   RemoveKnot(dimension,
                                              knot,
                                              tolerance_removal,
                                              multiplicity,
                                              tolerance));
  });
  return successful_removals;
}

template<int para_dim>
Spline<para_dim>::Spline(bool is_rational)
    : SplineItem(para_dim, std::move(is_rational)) {}

template<int para_dim>
Spline<para_dim>::Spline(SharedPointer<ParameterSpace_> parameter_space,
                         bool is_rational)
    : Spline(std::move(is_rational)) {
  static_assert(para_dim > 0, "The parametric dimensionality must be positive");

  parameter_space_ = std::move(parameter_space);
}

template<int para_dim>
Spline<para_dim>::Spline(Spline const& other)
    : Base_(),
      parameter_space_(
          std::make_shared<ParameterSpace_>(*other.parameter_space_)) {}

template<int para_dim>
Spline<para_dim>& Spline<para_dim>::operator=(Spline const& rhs) {
  Base_::operator=(rhs);
  parameter_space_ = std::make_shared<ParameterSpace_>(*rhs.parameter_space_);
  return *this;
}
