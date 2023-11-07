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

#ifndef SOURCES_SPLINES_NURBS_HPP_
#define SOURCES_SPLINES_NURBS_HPP_

#include <algorithm>
#include <deque>
#include <iterator>
#include <utility>

#include "BSplineLib/ParameterSpaces/parameter_space.hpp"
#include "BSplineLib/Splines/b_spline.hpp"
#include "BSplineLib/Splines/spline.hpp"
#include "BSplineLib/Utilities/containers.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/math_operations.hpp"
#include "BSplineLib/Utilities/named_type.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"
#include "BSplineLib/VectorSpaces/weighted_vector_space.hpp"

namespace bsplinelib::splines {

// NURBSs are rational B-splines.  Currently only single-patch NURBSs are
// supported.
template<int para_dim>
class Nurbs : public Spline<para_dim> {
public:
  using Base_ = Spline<para_dim>;
  using Coordinate_ = typename Base_::Coordinate_;
  using Derivative_ = typename Base_::Derivative_;
  using Knot_ = typename Base_::Knot_;
  using ParameterSpace_ = typename Base_::ParameterSpace_;
  using ParametricCoordinate_ = typename Base_::ParametricCoordinate_;
  using WeightedVectorSpace_ = vector_spaces::WeightedVectorSpace;
  using OutputInformation_ =
      Tuple<typename ParameterSpace_::OutputInformation_,
            typename WeightedVectorSpace_::OutputInformation_>;

  using Type_ = typename ParameterSpace_::Type_;
  using IntType_ = typename ParameterSpace_::IntType_;

  Nurbs();
  Nurbs(SharedPointer<ParameterSpace_> parameter_space,
        SharedPointer<WeightedVectorSpace_> weighted_vector_space);
  Nurbs(Nurbs const& other);
  Nurbs(Nurbs&& other) noexcept = default;
  Nurbs& operator=(Nurbs const& rhs);
  Nurbs& operator=(Nurbs&& rhs) noexcept = default;
  ~Nurbs() override = default;

  /// @brief nurbs keeps weighted vector space which will return dim + 1. For
  /// splines, we reverse that and return true dim.
  /// @return
  virtual int Dim() const override {
    return weighted_vector_space_->Dim() - 1;
  };

  void Evaluate(const Type_* parametric_coordinate, Type_* evaluated) const;
  void EvaluateDerivative(const Type_* parametric_coordinate,
                          const IntType_* derivative,
                          Type_* evaluated) const;

  Coordinate_ operator()(const Type_* parametric_coordinate) const final;
  Coordinate_ operator()(const Type_* parametric_coordinate,
                         const IntType_* derivative) const final;

  void InsertKnot(Dimension const& dimension,
                  Knot_ knot,
                  Multiplicity const& multiplicity = kMultiplicity,
                  Tolerance const& tolerance = kEpsilon) const final;
  Multiplicity RemoveKnot(Dimension const& dimension,
                          Knot_ const& knot,
                          Tolerance const& tolerance_removal,
                          Multiplicity const& multiplicity = kMultiplicity,
                          Tolerance const& tolerance = kEpsilon) const final;
  void ElevateDegree(Dimension const& dimension,
                     Multiplicity const& multiplicity = kMultiplicity,
                     Tolerance const& tolerance = kEpsilon) const final;
  bool ReduceDegree(Dimension const& dimension,
                    Tolerance const& tolerance_removal,
                    Multiplicity const& multiplicity = kMultiplicity,
                    Tolerance const& tolerance = kEpsilon) const final;

  Coordinate ComputeUpperBoundForMaximumDistanceFromOrigin() const final;
  OutputInformation_ Write(Precision const& precision = kPrecision) const;

  /// @brief
  /// @param from
  /// @param to
  virtual void ShareMembers(std::shared_ptr<Nurbs>& from) {
    homogeneous_b_spline_ = from->homogeneous_b_spline_;
    Base_::parameter_space_ = from->parameter_space_;
    weighted_vector_space_ = from->weighted_vector_space_;
  }

protected:
  using HomogeneousBSpline_ = BSpline<para_dim>;

  SharedPointer<HomogeneousBSpline_> homogeneous_b_spline_;
  SharedPointer<WeightedVectorSpace_> weighted_vector_space_;
};

#include "BSplineLib/Splines/nurbs.inl"

} // namespace bsplinelib::splines

#endif // SOURCES_SPLINES_NURBS_HPP_
