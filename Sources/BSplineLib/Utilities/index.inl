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

template<int size>
Index<size>::Index(const Length_& length,
                   const Value_& value,
                   const bool invalid)
    : invalid_{invalid},
      length_{length},
      value_{value} {

#ifndef NDEBUG
  try {
    ThrowIfValueIsInvalid(length, value);
  } catch (OutOfRange const& exception) {
    Throw(exception, "bsplinelib::utilities::Index::Index");
  }
#endif
}

template<int size>
Index<size>& Index<size>::operator+=(Value_ const& value) {
  containers::AddAndAssignToFirst(value_, value);
#ifndef NDEBUG
  Message const kName{"bsplinelib::utilities::Index::operator+="};

  try {
    ThrowIfValueIsInvalid(length_, value_);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  return *this;
}

template<int size>
Index<size>& Index<size>::operator-=(Value_ const& value) {
  containers::SubtractAndAssignToFirst(value_, value);
#ifndef NDEBUG
  Message const kName{"bsplinelib::utilities::Index::operator-="};

  try {
    ThrowIfValueIsInvalid(length_, value_);
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  return *this;
}

template<int size>
constexpr Index<size>& Index<size>::operator++() {
  for (Dimension dimension{}; dimension < size; ++dimension) {
    Length const& length = length_[dimension];
    Index_& value = value_[dimension];
    if (length == 0) {
      continue;
    } else if (value == (length - 1)) {
      value = 0;
    } else {
      ++value;
      break;
    }
  }
  if (value_ == Value_{})
    invalid_ = true; // Invalidate the index if it was incremented from the
                     // maximum index.
  return *this;
}

template<int size>
Index<size>& Index<size>::Increment(Dimension const& dimension) {
#ifndef NDEBUG
  Message const kName{"bsplinelib::utilities::Index::Increment"};
  DimensionBoundCheck(kName, dimension);
#endif
  Length const& length = length_[dimension];
  Index_& value = value_[dimension];
  if (length != 0) {
    if (value == (length - 1)) {
      value = 0;
    } else {
      ++value;
    }
  }
  return *this;
}

template<int size>
constexpr Index<size>& Index<size>::operator--() {
  if (value_ == Value_{}) { // Invalidate the index if it is decremented from
                            // the zero index.
    invalid_ = true;
  } else {
    for (Dimension dimension{}; dimension < size; ++dimension) {
      int const& length = length_[dimension];
      int& value = value_[dimension];
      if (length == 0) {
        continue;
      } else if (value == 0) {
        value = (length - 1);
      } else {
        --value;
        break;
      }
    }
  }
  return *this;
}

template<int size>
Index<size>& Index<size>::Decrement(Dimension const& dimension) {
#ifndef NDEBUG
  Message const kName{"bsplinelib::utilities::Index::Decrement"};
  DimensionBoundCheck(kName, dimension);
#endif
  int const& length = length_[dimension];
  int& value = value_[dimension];
  if (length != 0) {
    if (value == 0) {
      value = (length - 1);
    } else {
      --value;
    }
  }
  return *this;
}

template<int size>
Index<size> operator+(Index<size> const& lhs,
                      typename Index<size>::Value_ const& rhs) {
  Index sum{lhs};
#ifndef NDEBUG
  Message const kName{"bsplinelib::utilities::operator+::Index"};

  try {
#endif
    sum += rhs;
#ifndef NDEBUG
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  return sum;
}

template<int size>
Index<size> operator-(Index<size> const& lhs,
                      typename Index<size>::Value_ const& rhs) {
  Index difference{lhs};
#ifndef NDEBUG
  Message const kName{"bsplinelib::utilities::operator-::Index"};

  try {
#endif
    difference -= rhs;
#ifndef NDEBUG
  } catch (DomainError const& exception) {
    Throw(exception, kName);
  } catch (OutOfRange const& exception) {
    Throw(exception, kName);
  }
#endif
  return difference;
}

template<int size>
bool operator==(Index<size> const& lhs, Index<size> const& rhs) {
  return ((lhs.invalid_ == rhs.invalid_) && (lhs.length_ == rhs.length_)
          && (lhs.value_ == rhs.value_));
}

template<int size>
bool operator!=(Index<size> const& lhs, Index<size> const& rhs) {
  return (!(lhs == rhs));
}

template<int size>
constexpr bsplinelib::Index const&
Index<size>::operator[](Dimension const& dimension) const {
#ifndef NDEBUG
  DimensionBoundCheck("bsplinelib::utilities::Index::operator[]", dimension);
#endif
  return value_[dimension];
}

template<int size>
constexpr Index<size> Index<size>::First(const Length_& length) {
  return Index{length};
}

template<int size>
constexpr Index<size> Index<size>::Last(const Length_& length) {
  Value_ value;

  auto v_iter = value.begin();
  for (const int& len : length) {
    if (len != 0) {
      *v_iter = len - 1;
    } else {
      *v_iter = 0;
    }
    ++v_iter;
  }

  return Index{length, value};
}

template<int size>
constexpr Index<size> Index<size>::Behind(const Length_& length) {
  return ++Last(length);
}

template<int size>
constexpr Index<size> Index<size>::Before(const Length_& length) {
  return --First(length);
}

template<int size>
constexpr int Index<size>::GetIndex1d(const Length_& length,
                                      const Value_& value) {
  auto stride = [&length](const int& dim) {
    int s{1};
    for (int i{}; i < dim; ++i) {
      const int len = length[i];
      if (len != 0) {
        s *= len;
      }
    }
    return s;
  };

  int id{};
  for (int i{}; i < size; ++i) {
    id += value[i] * stride(i);
  }
  return id;
}

template<int size>
constexpr int Index<size>::GetTotalNumberOfIndices() const {
  return DetermineStride(length_, size);
}

template<int size>
constexpr typename Index<size>::Value_ Index<size>::GetIndex() const {
  return value_;
}

template<int size>
constexpr bsplinelib::Index Index<size>::GetIndex1d() const {
  Index_ index_1d{};
  for (int i{}; i < size; ++i) {
    index_1d += (value_[i] * DetermineStride(length_, i));
  };
  return index_1d;
}

template<int size>
constexpr int Index<size>::DetermineStride(Length_ const& length,
                                           Dimension const& dimension) const {
  int stride{1};
  for (int i{}; i < dimension; ++i) {
    const int& len = length[i];
    if (len != 0) {
      stride *= len;
    }
  }

  return stride;
}

#ifndef NDEBUG
template<int size>
void Index<size>::ThrowIfValueIsInvalid(Length_ const& length,
                                        Value_ const& value) {
  for (int i{}; i < size; ++i) {
    Message const dimension_string{"for dimension " + std::to_string(i) + ": "};
    if (value[i] < 0 || value[i] > std::max(0, length[i] - 1)) {
      throw OutOfRange(dimension_string);
    }
  }
}
#endif
