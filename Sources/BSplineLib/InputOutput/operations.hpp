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

#ifndef SOURCES_INPUTOUTPUT_OPERATIONS_HPP_
#define SOURCES_INPUTOUTPUT_OPERATIONS_HPP_

#include "BSplineLib/Splines/b_spline.hpp"
#include "BSplineLib/Splines/nurbs.hpp"
#include "BSplineLib/Splines/spline_item.hpp"
#include "BSplineLib/Utilities/containers.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/string_operations.hpp"

// Input and output operations such as 1.) determining the file format and
// reading spline files, 2.) casting shared_pointers that own spline items to
// shared_pointers that own actual splines, and 3.) writing coordinates
// three-dimensionally.
//
// Example:
//   Splines const &b_splines = Read(b_splines.irit);
//   SharedPointer<BSpline<2, 3>> const &b_spline = CastToSpline<2, 3,
//   false>(b_splines[0]);
namespace bsplinelib::input_output::operations {

enum class FileFormat { kInvalid, kIges, kIrit, kVtk, kXml };
template<int para_dim>
using BSpline = splines::BSpline<para_dim>;
template<int para_dim>
using Nurbs = splines::Nurbs<para_dim>;
using SplineEntry = SharedPointer<splines::SplineItem>;
using Splines = Vector<SplineEntry>;

FileFormat DetermineFileFormat(String const& file_name);
Splines Read(String const& file_name);

template<int para_dim, bool is_rational>
auto CastToSpline(SplineEntry const& spline_entry);

template<typename Coordinate>
String WriteCoordinate3d(Coordinate const& coordinate, String const& delimiter);

#include "BSplineLib/InputOutput/operations.inl"

} // namespace bsplinelib::input_output::operations

namespace bsplinelib::input_output {

template<int para_dim>
using BSpline = operations::BSpline<para_dim>;
template<int para_dim>
using Nurbs = operations::Nurbs<para_dim>;
using SplineEntry = operations::SplineEntry;
using Splines = operations::Splines;

} // namespace bsplinelib::input_output

#endif // SOURCES_INPUTOUTPUT_OPERATIONS_HPP_
