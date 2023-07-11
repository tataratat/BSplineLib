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
#include <cassert>
#include <functional>
#include <initializer_list>
#include <utility>

#include "BSplineLib/Utilities/container_operations.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"

namespace bsplinelib::vector_spaces {

template<int dimensionality>
class VectorSpace;

template<int dimensionality>
bool IsEqual(VectorSpace<dimensionality> const& lhs,
             VectorSpace<dimensionality> const& rhs,
             Tolerance const& tolerance = kEpsilon);
template<int dimensionality>
bool operator==(VectorSpace<dimensionality> const& lhs,
                VectorSpace<dimensionality> const& rhs);

struct CoordinateView {
  using value_type = double;
  static constexpr const bool is_coordinate_view = true;
  // constant pointer
  double* const begin_;
  const int& dim_;

  double* begin() const { return begin_; }
  double* end() const { return begin_ + dim_; }
  const int& size() const { return dim_; }
  double& operator[](const int& i) const { return *(begin_ + i); }
  Vector<double> as_vector() const {
    return Vector<double>(begin_, begin_ + dim_);
  };
};

struct ConstCoordinateView {
  using value_type = double;
  static constexpr const bool is_coordinate_view = true;
  // constant pointer
  const double* begin_;
  const int& dim_;

  const double* begin() const { return begin_; }
  const double* end() const { return begin_ + dim_; }
  const int& size() const { return dim_; }
  const double& operator[](const int& i) const { return *(begin_ + i); }
  Vector<double> as_vector() const {
    return Vector<double>(begin_, begin_ + dim_);
  };
};

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
template<int dimensionality>
class VectorSpace {
public:
  using Coordinate_ = Vector<double>;
  using Coordinates_ = Vector<double>;
  using NamedCoordinate_ = Vector<Coordinate>;
  using NestedCoordinates_ = Vector<Coordinate_>;

  VectorSpace() = default;
  explicit VectorSpace(Coordinates_ coordinates);
  VectorSpace(VectorSpace const& other) = default;
  VectorSpace(VectorSpace&& other) noexcept = default;
  VectorSpace& operator=(VectorSpace const& rhs) = default;
  VectorSpace& operator=(VectorSpace&& rhs) noexcept = default;
  virtual ~VectorSpace() = default;

  // Comparison based on tolerance.
  friend bool IsEqual<dimensionality>(VectorSpace const& lhs,
                                      VectorSpace const& rhs,
                                      Tolerance const& tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==<dimensionality>(VectorSpace const& lhs,
                                         VectorSpace const& rhs);
  virtual CoordinateView operator[](const int coordinate_id);
  virtual CoordinateView operator[](Index const& coordinate);
  virtual ConstCoordinateView operator[](const int coordinate_id) const;
  virtual Coordinate_ CreateCoordinate() const;
  virtual Coordinate_
  CreateCoordinate(std::initializer_list<double> coordinate) const;
  virtual Coordinates_& GetCoordinates() { return coordinates_; }
  virtual Coordinates_ const& GetCoordinates() const { return coordinates_; }

  virtual int GetNumberOfCoordinates() const;
  template<typename CoordinateType>
  // virtual void Replace(Index const& coordinate_index, CoordinateType
  // coordinate);
  void Replace(Index const& coordinate_index, CoordinateType coordinate);
  template<typename CoordinateType>
  // virtual void Insert(Index const& coordinate_index, CoordinateType
  // coordinate);
  void Insert(Index const& coordinate_index, CoordinateType coordinate);
  virtual void Erase(Index const& coordinate_index);

  virtual Coordinate DetermineMaximumDistanceFromOrigin(
      Tolerance const& tolerance = kEpsilon) const;

protected:
  int dim_ = dimensionality;
  Coordinates_ coordinates_;

private:
#ifndef NDEBUG
  void ThrowIfIndexIsInvalid(Index const& coordinate) const;
#endif
};

#include "BSplineLib/VectorSpaces/vector_space.inc"

} // namespace bsplinelib::vector_spaces

#endif // SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_
