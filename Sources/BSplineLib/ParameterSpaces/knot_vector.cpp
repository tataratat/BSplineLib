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

#include "BSplineLib/ParameterSpaces/knot_vector.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>

#include "BSplineLib/Utilities/containers.hpp"
#include "BSplineLib/Utilities/error_handling.hpp"
#include "BSplineLib/Utilities/numeric_operations.hpp"

namespace bsplinelib::parameter_spaces {

using Knot = KnotVector::Knot_;
using std::move, std::to_string;
#ifndef NDEBUG
using utilities::numeric_operations::ThrowIfToleranceIsNegative;
#endif

KnotVector::KnotVector(Knots_ knots, Tolerance const& tolerance)
    : knots_(std::move(knots)) {
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::KnotVector::KnotVector"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  }
#endif
}

Knot const& KnotVector::operator[](Index const& index) const {
#ifndef NDEBUG
  try {
    Index::ThrowIfNamedIntegerIsOutOfBounds(index, knots_.size() - 1);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::parameter_spaces::KnotVector::operator[]");
  }
#endif
  return knots_[index.Get()];
}
Knot const& KnotVector::operator[](int const& index) const {
  return knots_[index];
}

int KnotVector::GetSize() const { return knots_.size(); }

Knot const& KnotVector::GetFront() const { return knots_[0]; }

Knot const& KnotVector::GetBack() const { return knots_[knots_.size() - 1]; }

void KnotVector::UpdateKnot(const int id, Knot const& knot) {
  if (knots_[id - 1] > knot || knots_[id + 1] < knot) {
    throw RuntimeError(
        "KnotVector::UpdateKnot - updated knot must be non-decreasing.");
  }
  knots_[id] = knot;
}

void KnotVector::Scale(Knot const& min, Knot const& max) {
  if (max < min) {
    throw RuntimeError("KnotVector::Scale - max is smaller than min.");
  }
  const auto current_min = GetFront();
  const auto scale_factor = (max - min) / (GetBack() - current_min);
  for (auto& knot : knots_) {
    knot = ((knot - current_min) * scale_factor) + min;
  }
}

bool KnotVector::DoesParametricCoordinateEqualBack(
    Knot const& parametric_coordinate,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::KnotVector::"
                      "DoesParametricCoordinateEqualBack"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  }
#endif
  return std::abs(static_cast<Knot>(parametric_coordinate) - GetBack())
         < tolerance;
}

bool KnotVector::DoesParametricCoordinateEqualFrontOrBack(
    Knot const& parametric_coordinate,
    Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::KnotVector::"
                      "DoesParametricCoordinateEqualFrontOrBack"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  }
#endif
  return (std::abs(static_cast<Knot>(parametric_coordinate) - GetFront())
          < tolerance)
         || DoesParametricCoordinateEqualBack(parametric_coordinate, tolerance);
}

KnotSpan KnotVector::FindSpan(Knot const& parametric_coordinate,
                              Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::KnotVector::FindSpan"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfParametricCoordinateIsOutsideScope(parametric_coordinate, tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  ConstIterator_ const &knots_begin = knots_.begin(), &knots_end = knots_.end();
  auto is_less = [&](Knot const& knot,
                     Knot const& parametric_coordinate_value) {
    return (knot < parametric_coordinate_value);
  };
  return KnotSpan{static_cast<int>(
      std::distance(
          knots_begin,
          DoesParametricCoordinateEqualBack(parametric_coordinate, tolerance)
              ? std::lower_bound(knots_begin,
                                 knots_end,
                                 parametric_coordinate,
                                 is_less)
              : std::upper_bound(knots_begin,
                                 knots_end,
                                 parametric_coordinate,
                                 is_less))
      - 1)};
}

Multiplicity
KnotVector::DetermineMultiplicity(Knot const& parametric_coordinate,
                                  Tolerance const& tolerance) const {
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::KnotVector::"
                      "DetermineMultiplicity"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  }
#endif
  auto knot_iter = knots_.cbegin();

  // initialize multiplicity and loop around
  int multiplicity{};
  for (; knot_iter != knots_.cend();) {
    if (std::abs(parametric_coordinate - *knot_iter++) < tolerance) {
      ++multiplicity;
      while (std::abs(parametric_coordinate - *knot_iter++) < tolerance) {
        ++multiplicity;
      }
      break;
    }
  }

  return Multiplicity{multiplicity};
}

Vector<int> KnotVector::DetermineMultiplicities(const Knot_* knot_vector_data,
                                                const int knot_vector_size,
                                                Tolerance const& tolerance) {
  // a block to sort can come here

  // create return
  Vector<int> multiplicities;
  multiplicities.reserve(knot_vector_size);

  // initialize unique_knot - alternative is to have an iterator
  double unique_knot{knot_vector_data[0]};
  int multiplicity{};

  // again, assumed sorted kv
  for (int i{}; i < knot_vector_size; ++i) {
    const Knot_& knot = knot_vector_data[i];
    if (std::abs(knot - unique_knot) < tolerance) {
      ++multiplicity;
    } else {
      multiplicities.push_back(multiplicity);
      // not same -> update unique_knot and reset counter
      unique_knot = knot;
      multiplicity = 1;
    }
  }

  // Add last knot if repeated (e.g. open knot vector)
  if (std::abs(knot_vector_data[knot_vector_size - 1] - unique_knot)
      < tolerance) {
    multiplicities.push_back(multiplicity);
  }

  return multiplicities;
}

KnotVector::Knots_
KnotVector::GetUniqueKnots(Tolerance const& tolerance) const {
#ifndef NDEBUG
  try {
    ThrowIfToleranceIsNegative(tolerance);
  } catch (InvalidArgument const& exception) {
    Throw(exception,
          "bsplinelib::parameter_spaces::KnotVector::GetUniqueKnots");
  }
#endif
  Knots_ unique_knots;
  std::unique_copy(knots_.begin(),
                   knots_.end(),
                   std::back_inserter(unique_knots),
                   [&](Knot const& lhs_knot, Knot const& rhs_knot) {
                     return std::abs(lhs_knot - rhs_knot) < tolerance;
                   });
  return unique_knots;
}

void KnotVector::Insert(Knot knot,
                        Multiplicity const& multiplicity,
                        Tolerance const& tolerance) {
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::KnotVector::Insert"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfParametricCoordinateIsOutsideScope(ParametricCoordinate{knot},
                                              tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  knots_.insert(knots_.begin() + FindSpan(knot, tolerance).Get() + 1,
                multiplicity.Get(),
                std::move(knot));
}

Multiplicity KnotVector::Remove(Knot const& knot,
                                Multiplicity const& multiplicity,
                                Tolerance const& tolerance) {
#ifndef NDEBUG
  Message const kName{"bsplinelib::parameter_spaces::KnotVector::Remove"};

  try {
    ThrowIfToleranceIsNegative(tolerance);
    ThrowIfTooSmallOrNotNonDecreasing(tolerance);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (InvalidArgument const& exception) {
    Throw(exception, kName);
  }
#endif
  if (Multiplicity::Type_ const number_of_removals{
          std::min(multiplicity, DetermineMultiplicity(knot, tolerance)).Get()};
      number_of_removals != 0) {
    KnotSpan const& knot_span = FindSpan(knot, tolerance);
    if (DoesParametricCoordinateEqualBack(knot, tolerance)) {
      ConstIterator_ const& end = knots_.end();
      knots_.erase(end - number_of_removals, end);
    } else {
      ConstIterator_ const& first_knot = (knots_.begin() + knot_span.Get());
      knots_.erase(first_knot - (number_of_removals - 1), first_knot + 1);
    }
    return Multiplicity{number_of_removals};
  } else {
    return Multiplicity{};
  }
}

void KnotVector::IncreaseMultiplicities(Multiplicity const& multiplicity,
                                        Tolerance const& tolerance) {
  for (Knot const& knot : GetUniqueKnots(tolerance))
    Insert(knot, multiplicity, tolerance);
}

void KnotVector::DecreaseMultiplicities(Multiplicity const& multiplicity,
                                        Tolerance const& tolerance) {
  if (GetSize() > 2)
    for (Knot const& knot : GetUniqueKnots(tolerance))
      Remove(knot, multiplicity, tolerance);
}

typename KnotVector::OutputInformation_
KnotVector::Write(Precision const& precision) const {
  return utilities::string_operations::Write<OutputInformation_>(knots_,
                                                                 precision);
}

#ifndef NDEBUG
void KnotVector::ThrowIfParametricCoordinateIsOutsideScope(
    Knot const& parametric_coordinate,
    Tolerance const& tolerance) const {
  Knot const &first_knot = GetFront(), &last_knot = GetBack();
  if ((parametric_coordinate + tolerance) < first_knot
      || (parametric_coordinate - tolerance) > last_knot)
    throw OutOfRange(
        "The parametric coordinate " + to_string(parametric_coordinate)
        + " is outside of the "
          "knot vector's scope ["
        + to_string(first_knot) + "," + to_string(last_knot) + "].");
}
#endif

void KnotVector::ThrowIfTooSmallOrNotNonDecreasing(
    Tolerance const& tolerance) const {
  int const& number_of_knots = knots_.size();
  if (number_of_knots < 2)
    throw DomainError(
        "The knot vector has to contain at least 2 knots but only contains "
        + to_string(number_of_knots) + ".");
  Index::ForEach(1, number_of_knots, [&](Index const& knot) {
    Index::Type_ const& index = knot.Get();
    Knot const &current_knot = knots_[index],
               &previous_knot = knots_[index - 1];
    if ((current_knot + tolerance) < previous_knot)
      throw DomainError("The knot vector has to be a non-decreasing sequence "
                        "of real numbers but the knot "
                        + to_string(current_knot) + " at index "
                        + to_string(index)
                        + " is less than the "
                          "previous knot "
                        + to_string(previous_knot) + ".");
  });
}

} // namespace bsplinelib::parameter_spaces
