#include <array>
#include <vector>

#include "BSplineLib/Utilities/error_handling.hpp"

namespace bsplinelib::utilities::containers {

/// @brief lightweight self deleting array. Meant to be used for simple tasks,
/// where you want to avoid value initialization of every element, unlike
/// std::vector does.
/// @tparam T
template<typename T>
struct TemporaryArray {
  T* data_;
  TemporaryArray(const int n) : data_(new T[n]) {}
  ~TemporaryArray() { delete[] data_; }
  constexpr T& operator[](const int& i) { return data_[i]; }
  constexpr const T& operator[](const int& i) const { return data_[i]; }
};

/// @brief lightweight self deleting 2D array. Meant to be used for simple
/// tasks, where you want to avoid value initialization of every element, unlike
/// std::vector.
/// @tparam T
template<typename T>
struct TemporaryArray2D {
  T* data_;
  int dim_;
  TemporaryArray2D(const int n, const int d) : data_(new T[n * d]), dim_(d) {}
  ~TemporaryArray2D() { delete[] data_; }
  constexpr T& operator[](const int& i, const int& j) {
    return data_[i * dim_ + j];
  }
  constexpr const T& operator[](const int& i, const int& j) const {
    return data_[i * dim_ + j];
  }
};

/// @brief Fully dynamic array that can view another data. Equipped with basic
/// math operations.
/// @tparam DataType
/// @tparam dim
template<typename DataType, int dim = 1, typename IndexType = int>
class Array {
  static_assert(dim > 0, "dim needs to be positive value bigger than zero.");
  static_assert(std::is_integral_v<IndexType>,
                "IndexType should be an integral type");

public:
  using ShapeType_ = std::array<IndexType, dim>;
  using StridesType_ = std::array<IndexType, dim - 1>;
  static constexpr const IndexType kDim = static_cast<IndexType>(dim);

  using DataType_ = DataType;
  using IndexType_ = IndexType;

  // std container like types
  using value_type = DataType;
  using size_type = IndexType;

protected:
  bool own_data_{false};
  /// @brief beginning of the array pointer that this object manages
  DataType* data_{nullptr};

  /// @brief size of this Array
  IndexType size_{};

  /// @brief strides in case this is a higher dim. last entry should be the same
  /// as size_
  StridesType_ strides_;

  /// @brief shape of the array
  ShapeType_ shape_;

  /// @brief helper to find id offset
  /// @tparam ...Ts
  /// @param index
  /// @param counter
  /// @param ...id
  template<typename... Ts>
  constexpr void ComputeRemainingOffset(IndexType& index,
                                        IndexType& counter,
                                        const IndexType& id0,
                                        const Ts&... id) const {
    // last element should be just added
    if constexpr (sizeof...(Ts) == 0) {
      index += id0;
      // intentionally skipping
      // ++counter;
      return;
    }

    // add strides * id
    index += strides_[counter++] * id0;

    // do the rest
    if constexpr (sizeof...(Ts) > 0) {
      ComputeRemainingOffset(index, counter, id...);
    }
  }

public:
  constexpr DataType* data() {
    assert(data_);
    return data_;
  }

  constexpr const DataType* data() const {
    assert(data_);
    return data_;
  }

  constexpr DataType* begin() {
    assert(data_);
    return data_;
  }

  constexpr DataType* end() {
    assert(data_);
    assert(size_ > 0);
    return data_ + size_;
  }

  constexpr IndexType size() const { return size_; }

  constexpr void DestroyData() {
    if (own_data_ && data_) {
      delete[] data_;
    }

    data_ = nullptr;
    own_data_ = false;
  }

  constexpr void SetData(DataType* data_pointer) {
    DestroyData();
    data_ = data_pointer;
  }

  constexpr DataType* GetData() {
    assert(data_);
    return data_;
  }

  constexpr const DataType* GetData() const {
    assert(data_);
    return data_;
  }

  /// @brief reallocates space. After this call, own_data_ should be true
  /// @param size
  constexpr void Reallocate(const IndexType& size) {
    // destroy and reallocate space
    DestroyData();
    data_ = new DataType[size];
    own_data_ = true;

    // set size - don't forget to set shape in case this is multi-dim array
    size_ = size;

    // set shape if this is a 1d array
    if constexpr (dim == 1) {
      shape_[0] = size;
    }
  }

  template<typename... Ts>
  constexpr void SetShape(const IndexType& shape0, const Ts&... shape) {
    static_assert(sizeof...(Ts) == static_cast<std::size_t>(dim - 1),
                  "shape should match the dim");

    // early exit in case of dim 1
    if constexpr (sizeof...(Ts) == 0) {
      size_ = shape0;
      shape_[0] = shape0;
      return;
    }

    // (re)initialize shape_
    shape_ = {shape0, shape...};

    // (re)initialize  strides_
    strides_ = {shape...};

    // reverse iter and accumulate
    IndexType sub_total{1};
    typename StridesType_::reverse_iterator r_stride_iter = strides_.rbegin();
    for (; r_stride_iter != strides_.rend(); ++r_stride_iter) {
      // accumulate
      sub_total *= *r_stride_iter;
      // assign
      *r_stride_iter = sub_total;
    }

    // at this point, sub total has accumulated shape values except shape0
    // set size
    size_ = sub_total * shape0;
  }

  /// @brief default. use SetData() and SetShape<false>
  Array() = default;

  /// @brief basic array ctor
  /// @tparam ...Ts
  /// @param ...shape
  template<typename... Ts>
  Array(const Ts&... shape) {
    SetShape(shape...);
    Reallocate(size_);
  }

  /// @brief data wrapping array
  /// @tparam ...Ts
  /// @param data_pointer
  /// @param ...shape
  template<typename... Ts>
  Array(DataType* data_pointer, const Ts&... shape) {
    SetData(data_pointer);
    SetShape(shape...);
  }

  /// @brief copy ctor
  /// @param other
  Array(const Array& other) {
    // memory alloc
    Reallocate(other.size_);
    // copy data
    std::copy_n(other.data(), other.size_, data_);
    // copy shape
    shape_ = other.shape_;
    // copy strides
    strides_ = other.strides_;
  }

  /// @brief move ctor
  /// @param other
  Array(Array&& other) {
    own_data_ = true;
    data_ = std::move(other.data_);
    size_ = std::move(other.size_);
    strides_ = std::move(other.strides_);
    shape_ = std::move(other.shape_);
    other.own_data_ = false;
  }

  /// @brief copy assignment - need to make sure that you have enough space
  /// @param rhs
  /// @return
  constexpr Array& operator=(const Array& rhs) {
    // size check is crucial, so runtime check
    if (size_ != rhs.size_) {
      throw bsplinelib::utilities::error_handling::RuntimeError(
          "Array::operator=(const Array&) - size mismatch between rhs");
    }
    std::copy_n(rhs.data(), rhs.size(), data_);

    // copy shape and stride info
    shape_ = rhs.Shape();
    strides_ = rhs.Strides();

    return *this;
  }

  /// @brief move assignment. currently same as ctor.
  /// @param rhs
  /// @return
  constexpr Array& operator=(Array&& rhs) {
    own_data_ = true;
    data_ = std::move(rhs.data_);
    size_ = std::move(rhs.size_);
    strides_ = std::move(rhs.strides_);
    shape_ = std::move(rhs.shape_);
    rhs.own_data_ = false;
    return *this;
  }

  /// @brief
  ~Array() { DestroyData(); }

  /// @brief Returns const shape object
  /// @return
  constexpr const ShapeType_& Shape() const { return shape_; }

  /// @brief Returns const strides object
  /// @return
  constexpr const StridesType_& Strides() const { return strides_; }

  /// @brief flat indexed array access
  /// @param index
  /// @return
  constexpr DataType& operator[](const IndexType& index) {
    assert(index > -1 && index < size_);

    return data_[index];
  }

  /// @brief flat indexed array access
  /// @param index
  /// @return
  constexpr const DataType& operator[](const IndexType& index) const {
    assert(index > -1 && index < size_);

    return data_[index];
  }

  /// @brief multi-indexed array access
  /// @tparam ...Ts
  /// @param id0
  /// @param ...id
  /// @return
  template<typename... Ts>
  constexpr DataType& operator()(const IndexType& id0, const Ts&... id) {
    static_assert(sizeof...(Ts) == static_cast<std::size_t>(dim - 1),
                  "number of index should match dim");

    if constexpr (dim == 1) {
      return operator[](id0);
    }

    IndexType final_id{id0 * strides_[0]}, i{1};
    if constexpr (sizeof...(Ts) > 0) {
      ComputeRemainingOffset(final_id, i, id...);
    }

    assert(final_id > -1 && final_id < size_);

    return data_[final_id];
  }

  /// @brief multi-indexed array access
  /// @tparam ...Ts
  /// @param id0
  /// @param ...id
  /// @return
  template<typename... Ts>
  constexpr const DataType& operator()(const IndexType& id0,
                                       const Ts&... id) const {
    static_assert(sizeof...(Ts) == static_cast<std::size_t>(dim - 1),
                  "number of index should match dim");

    if constexpr (dim == 1) {
      return operator[](id0);
    }

    IndexType final_id{id0 * strides_[0]}, i{1};
    if constexpr (sizeof...(Ts) > 0) {

      ComputeRemainingOffset(final_id, i, id...);
    }

    assert(final_id > -1 && final_id < size_);

    return data_[final_id];
  }

  /// @brief this[:] = v
  /// @param v
  constexpr Array& Fill(const DataType& v) {
    assert(data_);

    std::fill_n(data_, size_, v);

    return *this;
  }

  /// @brief this[i] = a[i]
  /// @param a
  /// @return
  constexpr Array& Copy(const DataType* a) {
    assert(data_);

    std::copy_n(a, size_, data_);
  }

  /// @brief this[i] += a[i]
  /// @tparam Iterable
  /// @param a
  template<typename Iterable>
  constexpr Array& Add(const Iterable& a) {
    assert(a.size() == size_);
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] += a[i];
    }
    return *this;
  }

  /// @brief this[i] += a * b[i]
  /// @param a 
  /// @param b 
  /// @return 
  constexpr Array& Add(const DataType_& a, const DataType* b) {
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] += a * b[i];
    }
    return *this;
  }

  /// @brief this[i] *= a * this[i]
  /// @param a
  constexpr Array& Multiply(const DataType& a) {
    assert(data_);

    for (DataType* d = begin(); d != end();) {
      *d++ *= a;
    }

    return *this;
  }

  /// @brief this[i] = a * b[i]. Must set size and data beforehand.
  /// @param a
  /// @param v
  constexpr Array& MultiplyAssign(const DataType& a, const DataType* b) {
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] = a * v[i];
    }

    return *this;
  }

  /// @brief this[i] -= a[i]
  /// @tparam Iterable
  /// @param a
  template<typename Iterable>
  constexpr Array& Subtract(const Iterable& a) {
    assert(a.size() == size_);
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] -= a[i];
    }

    return *this;
  }

  /// @brief this[i] -= a[i]
  /// @tparam Iterable
  /// @param a
  constexpr Array& Subtract(const DataType* a) {
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] -= a[i];
    }

    return *this;
  }

  /// @brief this[i] = a[i] - this[i]. Must set size and data beforehand.
  /// @param a
  /// @param b
  constexpr Array& FlipSubtract(const DataType* a) {
    assert(data_);

    for (IndexType i{}; i < size_; ++i) {
      data_[i] = a[i] - data_[i];
    }

    return *this;
  }

  /// @brief c = (this) dot (a).
  /// @tparam Iterable
  /// @param a
  /// @return
  template<typename Iterable>
  constexpr DataType InnerProduct(const Iterable& a) const {
    static_assert(dim == 1,
                  "inner product is only applicable for dim=1 Array.");
    assert(a.size() == size_);
    assert(data_);

    DataType dot{};
    for (IndexType i{}; i < size_; ++i) {
      dot += a[i] * data_[i];
    }

    return dot;
  }

  constexpr void AAt(Array& aa_t) const {
    static_assert(dim == 2, "AAt is only applicable for dim=2 array.");

    const auto& height = shape_[0];
    const auto& width = shape_[1];

    for (int i{}; i < height; ++i) {
      // compute upper triangle
      for (int j{i}; j < height; ++j) {

        // element to fill
        DataType& ij = aa_t(i, j);
        // init
        ij = 0.0;
        // get beginning of i and j row
        const DataType* i_ptr = &operator()(i, 0);
        const DataType* j_ptr = &operator()(j, 0);

        // inner product
        for (int k{}; k < width; ++k) {
          ij += i_ptr[k] * j_ptr[k];
        }

        // fill symmetric part - this will write diagonal twice. it's okay.
        aa_t(j, i) = ij;
      }
    }
  }

  /// @brief given upper and lower bounds, clips data values inplace.
  /// @tparam Iterable
  /// @tparam IntegerType
  /// @param lower_bound
  /// @param upper_bound
  /// @param clipped
  template<typename Iterable, typename IntegerType>
  constexpr void Clip(const Iterable& lower_bound,
                      const Iterable& upper_bound,
                      Array<IntegerType, dim, IndexType>& clipped) {
    static_assert(std::is_integral_v<IntegerType>,
                  "clipped array should be an integral type");
    static_assert(std::is_signed_v<IntegerType>,
                  "clipped array should be a signed type");

    assert(lower_bound.size() == size_);
    assert(upper_bound.size() == size_);

    for (IndexType i{}; i < size_; ++i) {

      // deref values to write
      DataType& data = data_[i];
      IntegerType& clip = clipped[i];

      if (const DataType& ub = upper_bound[i]; data > ub) {
        // upper bound check
        clip = 1;
        data = ub;
      } else if (const DataType& lb = lower_bound[i]; data < lb) {
        // lower bound check
        clip = -1;
        data = lb;
      } else {
        // within the bounds
        clip = 0;
      }
    }
  }

  /// @brief L2 norm
  /// @return
  constexpr DataType NormL2() {
    DataType norm{};
    for (IndexType i{}; i < size_; ++i) {
      const DataType& data_i = data_[i];
      norm += data_i * data_i;
    }
    return std::sqrt(norm);
  }

  /// @brief returns number of non (exact) zero elements.
  /// @return
  constexpr IndexType NonZeros() const {
    IndexType n{};
    for (IndexType i{}; i < size_; ++i) {
      if (n != static_cast<DataType>(0)) {
        ++n;
      }
    }
    return n;
  }
};

/// @brief x[i] = a * v[i] + b * w[i].
/// @tparam DataType
/// @tparam ArrayType
/// @param a
/// @param v
/// @param b
/// @param w
/// @param x
template<typename DataType, typename ArrayType>
inline void Add(const DataType& a,
                const DataType* v,
                const DataType& b,
                const DataType* w,
                ArrayType& x) {
  const auto size = e.size();
  for (std::remove_const_t<delctype(size)> i{}; i < size; ++i) {
    x[i] = a * v[i] + c * w[i];
  }
}

/// Adapted from a post from Casey
/// http://stackoverflow.com/a/21028912/273767
/// mentioned in `Note` at
/// http://en.cppreference.com/w/cpp/container/vector/resize
///
/// comments, also from the post:
/// Allocator adaptor that interposes construct() calls to
/// convert value initialization into default initialization.
template<typename Type, typename BaseAllocator = std::allocator<Type>>
class DefaultInitializationAllocator : public BaseAllocator {
  using AllocatorTraits_ = std::allocator_traits<BaseAllocator>;

public:
  template<typename U>
  /// @brief Rebind
  struct rebind {
    using other = DefaultInitializationAllocator<
        U,
        typename AllocatorTraits_::template rebind_alloc<U>>;
  };

  using BaseAllocator::BaseAllocator;

  /// @brief Construct
  /// @tparam U
  /// @param ptr
  template<typename U>
  void
  construct(U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value) {
    ::new (static_cast<void*>(ptr)) U;
  }
  /// @brief Construct
  /// @tparam U
  /// @tparam ...Args
  /// @param ptr
  /// @param ...args
  template<typename U, typename... Args>
  void construct(U* ptr, Args&&... args) {
    AllocatorTraits_::construct(static_cast<BaseAllocator&>(*this),
                                ptr,
                                std::forward<Args>(args)...);
  }
};

/// @brief short-cut to vector that default initializes
/// @tparam Type
template<typename Type>
using DefaultInitializationVector =
    std::vector<Type, DefaultInitializationAllocator<Type>>;

} // namespace bsplinelib::utilities::containers
