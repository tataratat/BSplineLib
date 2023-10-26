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

#ifndef SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_
#define SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_

#include <algorithm>
#include <functional>
#include <utility>

#include "BSplineLib/Utilities/containers.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"
#include "BSplineLib/Utilities/string_operations.hpp"

namespace bsplinelib::vector_spaces {

// VectorSpaces group coordinates.
//
// Example:
//   using VectorSpace3d = VectorSpace<3>;
//   using Coordinate = VectorSpace3d::Coordinate_;
//   using ScalarCoordinate = Coordinate::value_type;
//   ScalarCoordinate const k0_0{}, k1_0{1.0};
//   VectorSpace3d const vector_space{{{k0_0, k0_0, k0_0}, {k1_0, k0_0, k0_0},
//   {k0_0, k1_0, k0_0}, {k1_0, k1_0, k0_0}}}; int const &four =
//   vector_space.GetNumberOfCoordinates(); Coordinate const &coordinate =
//   vector_space[Index{1}];  // Coordinate P_1 = {1.0, 0.0, 0.0}.
//   ScalarCoordinate const &one_point_zero =
//   vector_space.DetermineMaximumDistanceFromOrigin();
class VectorSpace {
public:
  using DataType_ = Coordinate;

  template<typename T>
  using Data_ = bsplinelib::utilities::containers::Data<T, 1>;
  template<typename T>
  using Data2D_ = bsplinelib::utilities::containers::Data<T, 2>;
  template<typename T>
  using Vector_ =
      bsplinelib::utilities::containers::DefaultInitializationVector<T>;

  using Coordinate_ = Data_<DataType_>;
  using ConstCoordinate_ = Data_<const DataType_>;
  using Coordinates_ = Data2D_<DataType_>;
  using OutputInformation_ = Tuple<Vector<StringVector>>;

  VectorSpace() = default;

  /// @brief coordinate copy ctor
  /// @param coordinates
  explicit VectorSpace(const Coordinates_& coordinates) {
    // copy
    coordinates_ = coordinates;
  }

  /// @brief coordinate move ctor
  /// @param coordinates
  explicit VectorSpace(Coordinates_&& coordinates) {
    coordinates_ = std::move(coordinates);
  }

  /// @brief data pointer ctor
  /// @param data
  /// @param shape0
  /// @param shape1
  explicit VectorSpace(Coordinate* data, const int& shape0, const int& shape1) {
    // take data and make a 2d view
    coordinates_.SetData(data);
    coordinates_.SetShape(shape0, shape1);
  }

  VectorSpace(VectorSpace const& other) = default;
  VectorSpace(VectorSpace&& other) noexcept = default;
  VectorSpace& operator=(VectorSpace const& rhs) = default;
  VectorSpace& operator=(VectorSpace&& rhs) noexcept = default;
  virtual ~VectorSpace() = default;

  /// @brief dim - shape[1] of coordinates_
  /// @return
  virtual int Dim() const { return coordinates_.Shape()[1]; }

  virtual Coordinate_ operator[](Index const& i) {
    return Coordinate_(&coordinates_(i.Get(), 0), coordinates_.Shape()[1]);
  }
  virtual ConstCoordinate_ operator[](Index const& i) const {
    return ConstCoordinate_(&coordinates_(i.Get(), 0), coordinates_.Shape()[1]);
  }
  virtual Coordinate_ operator[](const int& i) {
    return Coordinate_(&coordinates_(i, 0), coordinates_.Shape()[1]);
  }
  virtual ConstCoordinate_ operator[](const int& i) const {
    return ConstCoordinate_(&coordinates_(i, 0), coordinates_.Shape()[1]);
  }

  virtual DataType_* CoordinateBegin(const int& i) {
    return &coordinates_(i, 0);
  }

  virtual const DataType_* CoordinateBegin(const int& i) const {
    return &coordinates_(i, 0);
  }

  /// @brief coordinates getter - if you change size, call
  /// @return
  virtual Coordinates_& GetCoordinates() { return coordinates_; }

  /// @brief const coordinates getter
  /// @return
  virtual Coordinates_ const& GetCoordinates() const { return coordinates_; }

  /// @brief shape[0] of coordinates array
  /// @return
  virtual int GetNumberOfCoordinates() const { return coordinates_.Shape()[0]; }

  /// @brief Replace coordinate value
  /// @param coordinate_index
  /// @param coordinate
  virtual void Replace(int const& coordinate_index,
                       const Coordinate_& coordinate);

  /// @brief Inserts a coordinate. This will invalidate any iterator / pointer
  /// to existing coordinates
  /// @param coordinate_index
  /// @param coordinate
  virtual void ReallocateInsert(int const& coordinate_index,
                                const Coordinate_& coordinate);

  /// @brief Erases a coordinate. This will invalidate any iterator / pointer to
  /// existing coords. Doesn't reallocate.
  /// @param coordinate_index
  virtual void Erase(int const& coordinate_index);

  /// @brief Computes max norm of the coordinate.
  /// @param tolerance
  /// @return
  virtual DataType_ DetermineMaximumDistanceFromOrigin() const;

  /// @brief IO operation helper
  /// @param precision
  /// @return
  virtual OutputInformation_
  Write(Precision const& precision = kPrecision) const;

protected:
  /// 2D, contiguous array. For Insert and Erase, you need to own the data
  Coordinates_ coordinates_;
};

} // namespace bsplinelib::vector_spaces

#endif // SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_
