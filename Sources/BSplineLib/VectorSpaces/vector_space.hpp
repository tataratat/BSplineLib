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

class VectorSpace;

bool IsEqual(VectorSpace const& lhs,
             VectorSpace const& rhs,
             Tolerance const& tolerance = kEpsilon);

bool operator==(VectorSpace const& lhs, VectorSpace const& rhs);

struct CoordinateView {
  using value_type = double;
  static constexpr const bool is_coordinate_view = true;
  // constant pointer
  double* const begin_;
  const int dim_;

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
  const int dim_;

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
class VectorSpace {
public:
  using Coordinate_ = Vector<double>;
  using Coordinates_ = Vector<double>;
  using NamedCoordinate_ = Vector<Coordinate>;
  using NestedCoordinates_ = Vector<Coordinate_>;

  VectorSpace() = default;
  explicit VectorSpace(Coordinates_ coordinates, const int dimensionality);
  VectorSpace(VectorSpace const& other) = default;
  VectorSpace(VectorSpace&& other) noexcept = default;
  VectorSpace& operator=(VectorSpace const& rhs) = default;
  VectorSpace& operator=(VectorSpace&& rhs) noexcept = default;
  virtual ~VectorSpace() = default;

  // Comparison based on tolerance.
  friend bool IsEqual(VectorSpace const& lhs,
                      VectorSpace const& rhs,
                      Tolerance const& tolerance);
  // Comparison based on numeric_operations::GetEpsilon<Tolerance>().
  friend bool operator==(VectorSpace const& lhs, VectorSpace const& rhs);

  virtual void SetDimensionality(const int dimensionality);
  virtual int GetDimensionality() const;
  virtual CoordinateView operator[](const int coordinate_id);
  virtual CoordinateView operator[](Index const& coordinate);
  virtual ConstCoordinateView operator[](const int coordinate_id) const;
  // similar to operator[], but without checks.
  // safe to use iff you haven't got an error with GetDimensionality() and
  // GetCoordinates() GetNumberOfCoordinates() uses both.
  virtual CoordinateView DirectView(const int coordinate_id);
  virtual ConstCoordinateView DirectView(const int coordinate_id) const;
  virtual Coordinate_ CreateCoordinate() const;
  virtual Coordinate_
  CreateCoordinate(std::initializer_list<double> coordinate) const;
  virtual Coordinates_& GetCoordinates();
  virtual Coordinates_ const& GetCoordinates() const;

  virtual int GetNumberOfCoordinates() const;
  template<typename CoordinateType>
  void Replace(Index const& coordinate_index, CoordinateType coordinate);
  template<typename CoordinateType>
  void Insert(Index const& coordinate_index, CoordinateType coordinate);
  virtual void Erase(Index const& coordinate_index);

  virtual Coordinate DetermineMaximumDistanceFromOrigin(
      Tolerance const& tolerance = kEpsilon) const;

protected:
  int dimensionality_{-1};
  Coordinates_ coordinates_;

private:
#ifndef NDEBUG
  void ThrowIfIndexIsInvalid(Index const& coordinate) const;
#endif
};

// templated implementations should come here
template<typename CoordinateType>
void VectorSpace::Insert(Index const& coordinate_index,
                         CoordinateType coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::VectorSpace::Insert");
  }
#endif
  coordinates_.insert(coordinates_.begin()
                          + (coordinate_index.Get() * GetDimensionality()),
                      coordinate.begin(),
                      coordinate.end());
}

template<typename CoordinateType>
void VectorSpace::Replace(Index const& coordinate_index,
                          CoordinateType coordinate) {
#ifndef NDEBUG
  try {
    ThrowIfIndexIsInvalid(coordinate_index);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::vector_spaces::VectorSpace::Replace");
  }
#endif
  const auto this_coordinate = operator[](coordinate_index);
  for (int i{}; i < GetDimensionality(); ++i) {
    this_coordinate[i] = coordinate[i];
  }
}

// #include "BSplineLib/VectorSpaces/vector_space.inc"

} // namespace bsplinelib::vector_spaces

#endif // SOURCES_VECTORSPACES_VECTOR_SPACE_HPP_
