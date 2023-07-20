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

#include "BSplineLib/Splines/spline_item.hpp"

#include <utility>

#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"

namespace bsplinelib::splines {

SplineItem::SplineItem(int parametric_dimensionality,
                       int dimensionality,
                       bool is_rational)
    : dimensionality_(std::move(dimensionality)),
      parametric_dimensionality_(std::move(parametric_dimensionality)),
      is_rational_(std::move(is_rational)) {}

bool IsEqual(SplineItem const& lhs,
             SplineItem const& rhs,
             Tolerance const& tolerance) {
#ifndef NDEBUG
  try {
    utilities::numeric_operations::ThrowIfToleranceIsNegative(tolerance);
  } catch (RuntimeError const& exception) {
    Throw(exception, "bsplinelib::splines::IsEqual::SplineItem");
  }
#endif
  return ((lhs.parametric_dimensionality_ == rhs.parametric_dimensionality_)
          && (lhs.dimensionality_ == rhs.dimensionality_)
          && (lhs.is_rational_ == rhs.is_rational_));
}

bool operator==(SplineItem const& lhs, SplineItem const& rhs) {
  return IsEqual(lhs, rhs);
}

} // namespace bsplinelib::splines
