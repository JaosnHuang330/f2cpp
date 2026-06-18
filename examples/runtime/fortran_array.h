#ifndef FORTRAN_ARRAY_H
#define FORTRAN_ARRAY_H

#include <vector>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <algorithm>

// Non-owning Fortran-style array views.  These are used for assumed-shape
// dummy arguments and array sections such as a(:, :, k).  They preserve
// Fortran lower bounds and column-major strides without copying data.
template<typename T>
class FortranArray1DRef {
public:
    FortranArray1DRef() : data_(nullptr), n1_(0), lb1_(1), s1_(1) {}
    FortranArray1DRef(T* data, int n1, int lb1 = 1, int s1 = 1)
        : data_(data), n1_(n1), lb1_(lb1), s1_(s1) {}

    T& operator()(int i) { return data_[static_cast<size_t>(i - lb1_) * s1_]; }
    const T& operator()(int i) const { return data_[static_cast<size_t>(i - lb1_) * s1_]; }
    T& operator[](int i) { return (*this)(i); }
    const T& operator[](int i) const { return (*this)(i); }
    FortranArray1DRef& operator=(const T& val) {
        for (int i = lb1_; i <= ubound(1); ++i) (*this)(i) = val;
        return *this;
    }

    T* data() { return data_; }
    const T* data() const { return data_; }
    int size() const { return n1_; }
    int size(int) const { return n1_; }
    int extent(int dim) const { return size(dim); }
    size_t total_size() const { return static_cast<size_t>(n1_); }
    int lbound(int) const { return lb1_; }
    int ubound(int) const { return lb1_ + n1_ - 1; }
    int stride(int) const { return s1_; }
    bool is_contiguous() const { return n1_ <= 1 || s1_ == 1; }
    bool allocated() const { return data_ != nullptr; }
    FortranArray1DRef<T> slice(int lb, int ub) {
        T* base = data_ + static_cast<size_t>(lb - lb1_) * s1_;
        return FortranArray1DRef<T>(base, ub - lb + 1, lb, s1_);
    }

private:
    T* data_;
    int n1_;
    int lb1_;
    int s1_;
};

// LOGICAL arrays are stored as int values rather than std::vector<bool>.
// This matching view type keeps array sections and dummy arguments on the
// same storage model; using the generic bool* view would corrupt memory.
template<>
class FortranArray1DRef<bool> {
public:
    FortranArray1DRef() : data_(nullptr), n1_(0), lb1_(1), s1_(1) {}
    FortranArray1DRef(int* data, int n1, int lb1 = 1, int s1 = 1)
        : data_(data), n1_(n1), lb1_(lb1), s1_(s1) {}

    int& operator()(int i) { return data_[static_cast<size_t>(i - lb1_) * s1_]; }
    const int& operator()(int i) const { return data_[static_cast<size_t>(i - lb1_) * s1_]; }
    int& operator[](int i) { return (*this)(i); }
    const int& operator[](int i) const { return (*this)(i); }
    FortranArray1DRef& operator=(bool val) {
        for (int i = lb1_; i <= ubound(1); ++i) (*this)(i) = val ? 1 : 0;
        return *this;
    }

    int* data() { return data_; }
    const int* data() const { return data_; }
    int size() const { return n1_; }
    int size(int) const { return n1_; }
    int extent(int dim) const { return size(dim); }
    size_t total_size() const { return static_cast<size_t>(n1_); }
    int lbound(int) const { return lb1_; }
    int ubound(int) const { return lb1_ + n1_ - 1; }
    int stride(int) const { return s1_; }
    bool is_contiguous() const { return n1_ <= 1 || s1_ == 1; }
    bool allocated() const { return data_ != nullptr; }
    FortranArray1DRef<bool> slice(int lb, int ub) {
        int* base = data_ + static_cast<size_t>(lb - lb1_) * s1_;
        return FortranArray1DRef<bool>(base, ub - lb + 1, lb, s1_);
    }

private:
    int* data_;
    int n1_;
    int lb1_;
    int s1_;
};

template<typename T>
class FortranArray2DRef {
public:
    FortranArray2DRef() : data_(nullptr), n1_(0), n2_(0), lb1_(1), lb2_(1), s1_(1), s2_(0) {}
    FortranArray2DRef(T* data, int n1, int n2,
                      int lb1 = 1, int lb2 = 1,
                      int s1 = 1, int s2 = 0)
        : data_(data), n1_(n1), n2_(n2), lb1_(lb1), lb2_(lb2),
          s1_(s1), s2_(s2 ? s2 : n1) {}

    T& operator()(int i, int j) {
        return data_[static_cast<size_t>(i - lb1_) * s1_ +
                     static_cast<size_t>(j - lb2_) * s2_];
    }
    const T& operator()(int i, int j) const {
        return data_[static_cast<size_t>(i - lb1_) * s1_ +
                     static_cast<size_t>(j - lb2_) * s2_];
    }
    FortranArray2DRef& operator=(const T& val) {
        for (int j = lb2_; j <= ubound(2); ++j)
            for (int i = lb1_; i <= ubound(1); ++i)
                (*this)(i, j) = val;
        return *this;
    }

    T* data() { return data_; }
    const T* data() const { return data_; }
    int size1() const { return n1_; }
    int size2() const { return n2_; }
    int size(int dim) const { return dim == 1 ? n1_ : n2_; }
    int extent(int dim) const { return size(dim); }
    int lbound(int dim) const { return dim == 1 ? lb1_ : lb2_; }
    int ubound(int dim) const { return dim == 1 ? lb1_ + n1_ - 1 : lb2_ + n2_ - 1; }
    size_t total_size() const { return static_cast<size_t>(n1_) * n2_; }
    int stride(int dim) const { return dim == 1 ? s1_ : s2_; }
    bool is_contiguous() const {
        return n1_ == 0 || n2_ == 0 || ((n1_ <= 1 || s1_ == 1) && (n2_ <= 1 || s2_ == n1_));
    }
    bool allocated() const { return data_ != nullptr; }
    FortranArray2DRef<T> slice(int lb1, int ub1, int lb2, int ub2) {
        T* base = data_ + static_cast<size_t>(lb1 - lb1_) * s1_
                       + static_cast<size_t>(lb2 - lb2_) * s2_;
        return FortranArray2DRef<T>(base, ub1 - lb1 + 1, ub2 - lb2 + 1,
                                    lb1, lb2, s1_, s2_);
    }

private:
    T* data_;
    int n1_, n2_;
    int lb1_, lb2_;
    int s1_, s2_;
};

template<typename T>
class FortranArray3DRef {
public:
    FortranArray3DRef()
        : data_(nullptr), n1_(0), n2_(0), n3_(0),
          lb1_(1), lb2_(1), lb3_(1), s1_(1), s2_(0), s3_(0) {}
    FortranArray3DRef(T* data, int n1, int n2, int n3,
                      int lb1 = 1, int lb2 = 1, int lb3 = 1,
                      int s1 = 1, int s2 = 0, int s3 = 0)
        : data_(data), n1_(n1), n2_(n2), n3_(n3),
          lb1_(lb1), lb2_(lb2), lb3_(lb3),
          s1_(s1), s2_(s2 ? s2 : n1), s3_(s3 ? s3 : n1 * n2) {}

    T& operator()(int i, int j, int k) {
        return data_[static_cast<size_t>(i - lb1_) * s1_ +
                     static_cast<size_t>(j - lb2_) * s2_ +
                     static_cast<size_t>(k - lb3_) * s3_];
    }
    const T& operator()(int i, int j, int k) const {
        return data_[static_cast<size_t>(i - lb1_) * s1_ +
                     static_cast<size_t>(j - lb2_) * s2_ +
                     static_cast<size_t>(k - lb3_) * s3_];
    }
    FortranArray3DRef& operator=(const T& val) {
        for (int k = lb3_; k <= ubound(3); ++k)
            for (int j = lb2_; j <= ubound(2); ++j)
                for (int i = lb1_; i <= ubound(1); ++i)
                    (*this)(i, j, k) = val;
        return *this;
    }

    T* data() { return data_; }
    const T* data() const { return data_; }
    int size1() const { return n1_; }
    int size2() const { return n2_; }
    int size3() const { return n3_; }
    int size(int dim) const {
        if (dim == 1) return n1_;
        if (dim == 2) return n2_;
        return n3_;
    }
    int extent(int dim) const { return size(dim); }
    int lbound(int dim) const {
        if (dim == 1) return lb1_;
        if (dim == 2) return lb2_;
        return lb3_;
    }
    int ubound(int dim) const {
        if (dim == 1) return lb1_ + n1_ - 1;
        if (dim == 2) return lb2_ + n2_ - 1;
        return lb3_ + n3_ - 1;
    }
    size_t total_size() const { return static_cast<size_t>(n1_) * n2_ * n3_; }
    int stride(int dim) const {
        if (dim == 1) return s1_;
        if (dim == 2) return s2_;
        return s3_;
    }
    bool is_contiguous() const {
        return n1_ == 0 || n2_ == 0 || n3_ == 0 ||
               ((n1_ <= 1 || s1_ == 1) &&
                (n2_ <= 1 || s2_ == n1_) &&
                (n3_ <= 1 || s3_ == n1_ * n2_));
    }
    bool allocated() const { return data_ != nullptr; }
    FortranArray3DRef<T> slice(int lb1, int ub1, int lb2, int ub2, int lb3, int ub3) {
        T* base = data_ + static_cast<size_t>(lb1 - lb1_) * s1_
                       + static_cast<size_t>(lb2 - lb2_) * s2_
                       + static_cast<size_t>(lb3 - lb3_) * s3_;
        return FortranArray3DRef<T>(base, ub1 - lb1 + 1, ub2 - lb2 + 1, ub3 - lb3 + 1,
                                    lb1, lb2, lb3, s1_, s2_, s3_);
    }

private:
    T* data_;
    int n1_, n2_, n3_;
    int lb1_, lb2_, lb3_;
    int s1_, s2_, s3_;
};

template<typename T>
class FortranArray4DRef {
public:
    FortranArray4DRef()
        : data_(nullptr), n1_(0), n2_(0), n3_(0), n4_(0),
          lb1_(1), lb2_(1), lb3_(1), lb4_(1), s1_(1), s2_(0), s3_(0), s4_(0) {}
    FortranArray4DRef(T* data, int n1, int n2, int n3, int n4,
                      int lb1 = 1, int lb2 = 1, int lb3 = 1, int lb4 = 1,
                      int s1 = 1, int s2 = 0, int s3 = 0, int s4 = 0)
        : data_(data), n1_(n1), n2_(n2), n3_(n3), n4_(n4),
          lb1_(lb1), lb2_(lb2), lb3_(lb3), lb4_(lb4),
          s1_(s1), s2_(s2 ? s2 : n1), s3_(s3 ? s3 : n1 * n2),
          s4_(s4 ? s4 : n1 * n2 * n3) {}

    T& operator()(int i, int j, int k, int l) {
        return data_[static_cast<size_t>(i - lb1_) * s1_ +
                     static_cast<size_t>(j - lb2_) * s2_ +
                     static_cast<size_t>(k - lb3_) * s3_ +
                     static_cast<size_t>(l - lb4_) * s4_];
    }
    const T& operator()(int i, int j, int k, int l) const {
        return data_[static_cast<size_t>(i - lb1_) * s1_ +
                     static_cast<size_t>(j - lb2_) * s2_ +
                     static_cast<size_t>(k - lb3_) * s3_ +
                     static_cast<size_t>(l - lb4_) * s4_];
    }
    FortranArray4DRef& operator=(const T& val) {
        for (int l = lb4_; l <= ubound(4); ++l)
            for (int k = lb3_; k <= ubound(3); ++k)
                for (int j = lb2_; j <= ubound(2); ++j)
                    for (int i = lb1_; i <= ubound(1); ++i)
                        (*this)(i, j, k, l) = val;
        return *this;
    }

    T* data() { return data_; }
    const T* data() const { return data_; }
    int size(int dim) const {
        if (dim == 1) return n1_;
        if (dim == 2) return n2_;
        if (dim == 3) return n3_;
        return n4_;
    }
    int extent(int dim) const { return size(dim); }
    int lbound(int dim) const {
        if (dim == 1) return lb1_;
        if (dim == 2) return lb2_;
        if (dim == 3) return lb3_;
        return lb4_;
    }
    int ubound(int dim) const {
        if (dim == 1) return lb1_ + n1_ - 1;
        if (dim == 2) return lb2_ + n2_ - 1;
        if (dim == 3) return lb3_ + n3_ - 1;
        return lb4_ + n4_ - 1;
    }
    size_t total_size() const { return static_cast<size_t>(n1_) * n2_ * n3_ * n4_; }
    int stride(int dim) const {
        if (dim == 1) return s1_;
        if (dim == 2) return s2_;
        if (dim == 3) return s3_;
        return s4_;
    }
    bool is_contiguous() const {
        return n1_ == 0 || n2_ == 0 || n3_ == 0 || n4_ == 0 ||
               ((n1_ <= 1 || s1_ == 1) &&
                (n2_ <= 1 || s2_ == n1_) &&
                (n3_ <= 1 || s3_ == n1_ * n2_) &&
                (n4_ <= 1 || s4_ == n1_ * n2_ * n3_));
    }
    bool allocated() const { return data_ != nullptr; }
    FortranArray4DRef<T> slice(int lb1, int ub1, int lb2, int ub2,
                               int lb3, int ub3, int lb4, int ub4) {
        T* base = data_ + static_cast<size_t>(lb1 - lb1_) * s1_
                       + static_cast<size_t>(lb2 - lb2_) * s2_
                       + static_cast<size_t>(lb3 - lb3_) * s3_
                       + static_cast<size_t>(lb4 - lb4_) * s4_;
        return FortranArray4DRef<T>(base, ub1 - lb1 + 1, ub2 - lb2 + 1,
                                    ub3 - lb3 + 1, ub4 - lb4 + 1,
                                    lb1, lb2, lb3, lb4, s1_, s2_, s3_, s4_);
    }

private:
    T* data_;
    int n1_, n2_, n3_, n4_;
    int lb1_, lb2_, lb3_, lb4_;
    int s1_, s2_, s3_, s4_;
};

// Allocatable Fortran scalar, including allocatable derived types.
template<typename T>
class FortranScalar {
public:
    FortranScalar() : allocated_(false), value_{} {}

    void allocate() { allocated_ = true; value_ = T{}; }
    void allocate(int) { allocate(); }
    void deallocate() { allocated_ = false; value_ = T{}; }
    bool allocated() const { return allocated_; }

    T& operator()() { allocated_ = true; return value_; }
    const T& operator()() const { return value_; }
    T& operator()(int) { allocated_ = true; return value_; }
    const T& operator()(int) const { return value_; }

    operator T&() { allocated_ = true; return value_; }
    operator const T&() const { return value_; }
    T* operator&() { allocated_ = true; return &value_; }
    const T* operator&() const { return &value_; }
    T* data() { allocated_ = true; return &value_; }
    const T* data() const { return &value_; }

private:
    bool allocated_;
    T value_;
};

// 1D Fortran-style array with 1-based indexing (default) or custom lower bound
template<typename T>
class FortranArray1D {
public:
    FortranArray1D() : n1_(0), lb_(1) {}
    explicit FortranArray1D(int n1, T init = T{})
        : n1_(n1), lb_(1), data_(static_cast<size_t>(n1), init) {}

    void allocate(int n1) {
        n1_ = n1;
        lb_ = 1;
        data_.assign(static_cast<size_t>(n1), T{});
    }
    void allocate(int n1, T init) {
        n1_ = n1;
        lb_ = 1;
        data_.assign(static_cast<size_t>(n1), init);
    }
    // Allocate with custom lower bound: allocate_range(lb, ub) gives indices lb..ub
    void allocate_range(int lb, int ub) {
        lb_ = lb;
        n1_ = ub - lb + 1;
        data_.assign(static_cast<size_t>(n1_), T{});
    }
    void deallocate() { data_.clear(); n1_ = 0; lb_ = 1; }
    bool allocated() const { return !data_.empty(); }

    T& operator()(int i) { return data_[i - lb_]; }
    const T& operator()(int i) const { return data_[i - lb_]; }
    T& operator[](int i) { return (*this)(i); }
    const T& operator[](int i) const { return (*this)(i); }
    // Scalar assignment: arr = val fills all elements
    FortranArray1D& operator=(T val) { for (auto& x : data_) x = val; return *this; }

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }
    int size() const { return n1_; }
    int size(int) const { return n1_; }
    int lbound() const { return lb_; }
    int lbound(int) const { return lb_; }
    int ubound() const { return lb_ + n1_ - 1; }
    int ubound(int) const { return lb_ + n1_ - 1; }
    size_t total_size() const { return data_.size(); }
    FortranArray1DRef<T> ref() { return FortranArray1DRef<T>(data(), n1_, lb_, 1); }
    operator FortranArray1DRef<T>() { return ref(); }
    // Fortran array section: arr(lb:ub) → arr.slice(lb, ub)
    FortranArray1DRef<T> slice(int lb, int ub) {
        return FortranArray1DRef<T>(data() + (lb - lb_), ub - lb + 1, lb, 1);
    }

private:
    int n1_;
    int lb_;
    std::vector<T> data_;
};

// Specialization for bool to avoid std::vector<bool> issues
// ARPACK and LAPACK use int* for logical arrays, so we store as int
template<>
class FortranArray1D<bool> {
public:
    FortranArray1D() : n1_(0), lb_(1) {}
    explicit FortranArray1D(int n1, bool init = false)
        : n1_(n1), lb_(1), data_(static_cast<size_t>(n1), init ? 1 : 0) {}

    void allocate(int n1, bool init = false) {
        n1_ = n1;
        lb_ = 1;
        data_.assign(static_cast<size_t>(n1), init ? 1 : 0);
    }
    void allocate_range(int lb, int ub, bool init = false) {
        lb_ = lb;
        n1_ = ub - lb + 1;
        data_.assign(static_cast<size_t>(n1_), init ? 1 : 0);
    }
    void deallocate() { data_.clear(); n1_ = 0; lb_ = 1; }
    bool allocated() const { return !data_.empty(); }

    int& operator()(int i) { return data_[i - lb_]; }
    const int& operator()(int i) const { return data_[i - lb_]; }
    int& operator[](int i) { return (*this)(i); }
    const int& operator[](int i) const { return (*this)(i); }
    FortranArray1D& operator=(bool val) {
        for (auto& x : data_) x = val ? 1 : 0;
        return *this;
    }

    int* data() { return data_.data(); }
    const int* data() const { return data_.data(); }
    int size() const { return n1_; }
    int size(int) const { return n1_; }
    int lbound() const { return lb_; }
    int lbound(int) const { return lb_; }
    int ubound() const { return lb_ + n1_ - 1; }
    int ubound(int) const { return lb_ + n1_ - 1; }
    size_t total_size() const { return data_.size(); }
    FortranArray1DRef<bool> ref() { return FortranArray1DRef<bool>(data(), n1_, lb_, 1); }
    operator FortranArray1DRef<bool>() { return ref(); }
    FortranArray1DRef<bool> slice(int lb, int ub) {
        return FortranArray1DRef<bool>(data() + (lb - lb_), ub - lb + 1, lb, 1);
    }

private:
    int n1_;
    int lb_;
    std::vector<int> data_;
};

// 2D Fortran-style array with 1-based indexing, column-major storage
template<typename T>
class FortranArray2D {
public:
    FortranArray2D() : n1_(0), n2_(0) {}
    FortranArray2D(int n1, int n2, T init = T{})
        : n1_(n1), n2_(n2), data_(static_cast<size_t>(n1)*n2, init) {}
    // Range constructor for non-1 lower bounds
    FortranArray2D(int lb1, int ub1, int lb2, int ub2, T init = T{})
        : n1_(ub1-lb1+1), n2_(ub2-lb2+1), lb1_(lb1), lb2_(lb2),
          data_(static_cast<size_t>(n1_)*n2_, init) {}

    void allocate(int n1, int n2, T init = T{}) {
        n1_ = n1; n2_ = n2; lb1_ = 1; lb2_ = 1;
        data_.assign(static_cast<size_t>(n1)*n2, init);
    }
    void allocate_range(int lb1, int ub1, int lb2, int ub2, T init = T{}) {
        lb1_ = lb1; lb2_ = lb2;
        n1_ = ub1 - lb1 + 1; n2_ = ub2 - lb2 + 1;
        data_.assign(static_cast<size_t>(n1_)*n2_, init);
    }
    void deallocate() { data_.clear(); n1_ = n2_ = 0; lb1_ = lb2_ = 1; }
    bool allocated() const { return !data_.empty(); }

    T& operator()(int i, int j) {
        return data_[static_cast<size_t>(i-lb1_) + static_cast<size_t>(j-lb2_)*n1_];
    }
    const T& operator()(int i, int j) const {
        return data_[static_cast<size_t>(i-lb1_) + static_cast<size_t>(j-lb2_)*n1_];
    }
    FortranArray2D& operator=(T val) { for (auto& x : data_) x = val; return *this; }

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }
    int size1() const { return n1_; }
    int size2() const { return n2_; }
    int size(int dim) const { return dim == 1 ? n1_ : n2_; }
    int extent(int dim) const { return size(dim); }
    int lbound(int dim) const { return dim == 1 ? lb1_ : lb2_; }
    int ubound(int dim) const { return dim == 1 ? lb1_+n1_-1 : lb2_+n2_-1; }
    size_t total_size() const { return data_.size(); }
    FortranArray2DRef<T> ref() {
        return FortranArray2DRef<T>(data(), n1_, n2_, lb1_, lb2_, 1, n1_);
    }
    operator FortranArray2DRef<T>() { return ref(); }
    FortranArray1DRef<T> slice_dim1(int i) {
        T* base = data() + static_cast<size_t>(i - lb1_);
        return FortranArray1DRef<T>(base, n2_, lb2_, n1_);
    }
    FortranArray1DRef<T> slice_dim2(int j) {
        T* base = data() + static_cast<size_t>(j - lb2_) * n1_;
        return FortranArray1DRef<T>(base, n1_, lb1_, 1);
    }
    FortranArray2DRef<T> slice(int lb1, int ub1, int lb2, int ub2) {
        T* base = data() + static_cast<size_t>(lb1 - lb1_)
                         + static_cast<size_t>(lb2 - lb2_) * n1_;
        return FortranArray2DRef<T>(base, ub1 - lb1 + 1, ub2 - lb2 + 1,
                                    lb1, lb2, 1, n1_);
    }

private:
    int n1_, n2_;
    int lb1_ = 1, lb2_ = 1;
    std::vector<T> data_;
};

// 3D Fortran-style array with 1-based indexing, column-major storage
template<typename T>
class FortranArray3D {
public:
    FortranArray3D() : n1_(0), n2_(0), n3_(0) {}
    FortranArray3D(int n1, int n2, int n3, T init = T{})
        : n1_(n1), n2_(n2), n3_(n3),
          data_(static_cast<size_t>(n1)*n2*n3, init) {}
    // Range constructor for non-1 lower bounds (ghost cells, etc.)
    FortranArray3D(int lb1, int ub1, int lb2, int ub2, int lb3, int ub3, T init = T{})
        : n1_(ub1-lb1+1), n2_(ub2-lb2+1), n3_(ub3-lb3+1),
          lb1_(lb1), lb2_(lb2), lb3_(lb3),
          data_(static_cast<size_t>(n1_)*n2_*n3_, init) {}

    void allocate(int n1, int n2, int n3, T init = T{}) {
        n1_ = n1; n2_ = n2; n3_ = n3; lb1_ = 1; lb2_ = 1; lb3_ = 1;
        data_.assign(static_cast<size_t>(n1)*n2*n3, init);
    }
    void allocate_range(int lb1, int ub1, int lb2, int ub2, int lb3, int ub3, T init = T{}) {
        lb1_ = lb1; lb2_ = lb2; lb3_ = lb3;
        n1_ = ub1 - lb1 + 1; n2_ = ub2 - lb2 + 1; n3_ = ub3 - lb3 + 1;
        data_.assign(static_cast<size_t>(n1_)*n2_*n3_, init);
    }
    void deallocate() { data_.clear(); n1_ = n2_ = n3_ = 0; lb1_ = lb2_ = lb3_ = 1; }
    bool allocated() const { return !data_.empty(); }

    // Column-major with custom lower bounds
    T& operator()(int i, int j, int k) {
        return data_[static_cast<size_t>(i-lb1_) +
                     static_cast<size_t>(j-lb2_)*n1_ +
                     static_cast<size_t>(k-lb3_)*n1_*n2_];
    }
    const T& operator()(int i, int j, int k) const {
        return data_[static_cast<size_t>(i-lb1_) +
                     static_cast<size_t>(j-lb2_)*n1_ +
                     static_cast<size_t>(k-lb3_)*n1_*n2_];
    }
    FortranArray3D& operator=(T val) { for (auto& x : data_) x = val; return *this; }

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }
    int size1() const { return n1_; }
    int size2() const { return n2_; }
    int size3() const { return n3_; }
    int size(int dim) const {
        if (dim == 1) return n1_;
        if (dim == 2) return n2_;
        return n3_;
    }
    int extent(int dim) const { return size(dim); }
    int lbound(int dim) const {
        if (dim == 1) return lb1_;
        if (dim == 2) return lb2_;
        return lb3_;
    }
    int ubound(int dim) const {
        if (dim == 1) return lb1_ + n1_ - 1;
        if (dim == 2) return lb2_ + n2_ - 1;
        return lb3_ + n3_ - 1;
    }
    size_t total_size() const { return data_.size(); }
    FortranArray3DRef<T> ref() {
        return FortranArray3DRef<T>(data(), n1_, n2_, n3_,
                                    lb1_, lb2_, lb3_, 1, n1_, n1_ * n2_);
    }
    operator FortranArray3DRef<T>() { return ref(); }

    FortranArray2DRef<T> slice_dim1(int i) {
        T* base = data() + static_cast<size_t>(i - lb1_);
        return FortranArray2DRef<T>(base, n2_, n3_, lb2_, lb3_, n1_, n1_ * n2_);
    }
    FortranArray2DRef<T> slice_dim2(int j) {
        T* base = data() + static_cast<size_t>(j - lb2_) * n1_;
        return FortranArray2DRef<T>(base, n1_, n3_, lb1_, lb3_, 1, n1_ * n2_);
    }
    FortranArray2DRef<T> slice_dim3(int k) {
        T* base = data() + static_cast<size_t>(k - lb3_) * n1_ * n2_;
        return FortranArray2DRef<T>(base, n1_, n2_, lb1_, lb2_, 1, n1_);
    }
    FortranArray3DRef<T> slice(int lb1, int ub1, int lb2, int ub2, int lb3, int ub3) {
        T* base = data() + static_cast<size_t>(lb1 - lb1_)
                         + static_cast<size_t>(lb2 - lb2_) * n1_
                         + static_cast<size_t>(lb3 - lb3_) * n1_ * n2_;
        return FortranArray3DRef<T>(base, ub1 - lb1 + 1, ub2 - lb2 + 1, ub3 - lb3 + 1,
                                    lb1, lb2, lb3, 1, n1_, n1_ * n2_);
    }

    // Fortran-style slice: copy sub-array along dim 2
    FortranArray3D<T> slice_dim2(int j1, int j2) const {
        int nj = j2 - j1 + 1;
        FortranArray3D<T> result(n1_, nj, n3_);
        for (int k = lb3_; k <= lb3_+n3_-1; ++k)
            for (int j = j1; j <= j2; ++j)
                for (int i = lb1_; i <= lb1_+n1_-1; ++i)
                    result(i-lb1_+1, j-j1+1, k-lb3_+1) = (*this)(i, j, k);
        return result;
    }

    void set_slice_dim2(int j1, int j2, const FortranArray3D<T>& src) {
        for (int k = lb3_; k <= lb3_+n3_-1; ++k)
            for (int j = j1; j <= j2; ++j)
                for (int i = lb1_; i <= lb1_+n1_-1; ++i)
                    (*this)(i, j, k) = src(i-lb1_+1, j-j1+1, k-lb3_+1);
    }

private:
    int n1_, n2_, n3_;
    int lb1_ = 1, lb2_ = 1, lb3_ = 1;
    std::vector<T> data_;
};

// 4D Fortran-style array with 1-based indexing, column-major storage
template<typename T>
class FortranArray4D {
public:
    FortranArray4D() : n1_(0), n2_(0), n3_(0), n4_(0) {}
    FortranArray4D(int n1, int n2, int n3, int n4, T init = T{})
        : n1_(n1), n2_(n2), n3_(n3), n4_(n4),
          data_(static_cast<size_t>(n1)*n2*n3*n4, init) {}

    void allocate(int n1, int n2, int n3, int n4, T init = T{}) {
        n1_ = n1; n2_ = n2; n3_ = n3; n4_ = n4;
        lb1_ = 1; lb2_ = 1; lb3_ = 1; lb4_ = 1;
        data_.assign(static_cast<size_t>(n1)*n2*n3*n4, init);
    }
    void allocate_range(int lb1, int ub1, int lb2, int ub2,
                        int lb3, int ub3, int lb4, int ub4, T init = T{}) {
        lb1_ = lb1; lb2_ = lb2; lb3_ = lb3; lb4_ = lb4;
        n1_ = ub1 - lb1 + 1; n2_ = ub2 - lb2 + 1;
        n3_ = ub3 - lb3 + 1; n4_ = ub4 - lb4 + 1;
        data_.assign(static_cast<size_t>(n1_)*n2_*n3_*n4_, init);
    }
    void deallocate() { data_.clear(); n1_ = n2_ = n3_ = n4_ = 0; lb1_ = lb2_ = lb3_ = lb4_ = 1; }
    bool allocated() const { return !data_.empty(); }

    T& operator()(int i, int j, int k, int l) {
        return data_[static_cast<size_t>(i-lb1_) +
                     static_cast<size_t>(j-lb2_)*n1_ +
                     static_cast<size_t>(k-lb3_)*n1_*n2_ +
                     static_cast<size_t>(l-lb4_)*n1_*n2_*n3_];
    }
    const T& operator()(int i, int j, int k, int l) const {
        return data_[static_cast<size_t>(i-lb1_) +
                     static_cast<size_t>(j-lb2_)*n1_ +
                     static_cast<size_t>(k-lb3_)*n1_*n2_ +
                     static_cast<size_t>(l-lb4_)*n1_*n2_*n3_];
    }
    FortranArray4D& operator=(T val) { for (auto& x : data_) x = val; return *this; }

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }
    int size1() const { return n1_; }
    int size2() const { return n2_; }
    int size3() const { return n3_; }
    int size4() const { return n4_; }
    int size(int dim) const {
        if (dim == 1) return n1_;
        if (dim == 2) return n2_;
        if (dim == 3) return n3_;
        return n4_;
    }
    int lbound(int dim) const {
        if (dim == 1) return lb1_;
        if (dim == 2) return lb2_;
        if (dim == 3) return lb3_;
        return lb4_;
    }
    int ubound(int dim) const {
        if (dim == 1) return lb1_ + n1_ - 1;
        if (dim == 2) return lb2_ + n2_ - 1;
        if (dim == 3) return lb3_ + n3_ - 1;
        return lb4_ + n4_ - 1;
    }
    size_t total_size() const { return data_.size(); }
    FortranArray4DRef<T> ref() {
        return FortranArray4DRef<T>(data(), n1_, n2_, n3_, n4_,
                                    lb1_, lb2_, lb3_, lb4_, 1, n1_, n1_ * n2_, n1_ * n2_ * n3_);
    }
    operator FortranArray4DRef<T>() { return ref(); }
    FortranArray3DRef<T> slice_dim1(int i) {
        T* base = data() + static_cast<size_t>(i - lb1_);
        return FortranArray3DRef<T>(base, n2_, n3_, n4_, lb2_, lb3_, lb4_, n1_, n1_ * n2_, n1_ * n2_ * n3_);
    }
    FortranArray3DRef<T> slice_dim2(int j) {
        T* base = data() + static_cast<size_t>(j - lb2_) * n1_;
        return FortranArray3DRef<T>(base, n1_, n3_, n4_, lb1_, lb3_, lb4_, 1, n1_ * n2_, n1_ * n2_ * n3_);
    }
    FortranArray3DRef<T> slice_dim3(int k) {
        T* base = data() + static_cast<size_t>(k - lb3_) * n1_ * n2_;
        return FortranArray3DRef<T>(base, n1_, n2_, n4_, lb1_, lb2_, lb4_, 1, n1_, n1_ * n2_ * n3_);
    }
    FortranArray3DRef<T> slice_dim4(int l) {
        T* base = data() + static_cast<size_t>(l - lb4_) * n1_ * n2_ * n3_;
        return FortranArray3DRef<T>(base, n1_, n2_, n3_, lb1_, lb2_, lb3_, 1, n1_, n1_ * n2_);
    }
    FortranArray4DRef<T> slice(int lb1, int ub1, int lb2, int ub2,
                               int lb3, int ub3, int lb4, int ub4) {
        T* base = data() + static_cast<size_t>(lb1 - lb1_)
                         + static_cast<size_t>(lb2 - lb2_) * n1_
                         + static_cast<size_t>(lb3 - lb3_) * n1_ * n2_
                         + static_cast<size_t>(lb4 - lb4_) * n1_ * n2_ * n3_;
        return FortranArray4DRef<T>(base, ub1 - lb1 + 1, ub2 - lb2 + 1,
                                    ub3 - lb3 + 1, ub4 - lb4 + 1,
                                    lb1, lb2, lb3, lb4, 1, n1_, n1_ * n2_,
                                    n1_ * n2_ * n3_);
    }

private:
    int n1_, n2_, n3_, n4_;
    int lb1_ = 1, lb2_ = 1, lb3_ = 1, lb4_ = 1;
    std::vector<T> data_;
};

template<typename T>
FortranArray1D<T> operator*(FortranArray1DRef<T> arr, const T& scalar) {
    FortranArray1D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1));
    for (int i = arr.lbound(1); i <= arr.ubound(1); ++i) out(i) = arr(i) * scalar;
    return out;
}
template<typename T>
FortranArray1D<T> operator*(const T& scalar, FortranArray1DRef<T> arr) { return arr * scalar; }
template<typename T>
FortranArray1D<T> operator/(FortranArray1DRef<T> arr, const T& scalar) {
    FortranArray1D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1));
    for (int i = arr.lbound(1); i <= arr.ubound(1); ++i) out(i) = arr(i) / scalar;
    return out;
}
template<typename T>
FortranArray1D<T> operator*(FortranArray1D<T>& arr, const T& scalar) { return arr.ref() * scalar; }
template<typename T>
FortranArray1D<T> operator*(const T& scalar, FortranArray1D<T>& arr) { return arr.ref() * scalar; }
template<typename T>
FortranArray1D<T> operator/(FortranArray1D<T>& arr, const T& scalar) { return arr.ref() / scalar; }

template<typename T>
FortranArray2D<T> operator*(FortranArray2DRef<T> arr, const T& scalar) {
    FortranArray2D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1), arr.lbound(2), arr.ubound(2));
    for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
        for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
            out(i, j) = arr(i, j) * scalar;
    return out;
}
template<typename T>
FortranArray2D<T> operator*(const T& scalar, FortranArray2DRef<T> arr) { return arr * scalar; }
template<typename T>
FortranArray2D<T> operator/(FortranArray2DRef<T> arr, const T& scalar) {
    FortranArray2D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1), arr.lbound(2), arr.ubound(2));
    for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
        for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
            out(i, j) = arr(i, j) / scalar;
    return out;
}
template<typename T>
FortranArray2D<T> operator*(FortranArray2D<T>& arr, const T& scalar) { return arr.ref() * scalar; }
template<typename T>
FortranArray2D<T> operator*(const T& scalar, FortranArray2D<T>& arr) { return arr.ref() * scalar; }
template<typename T>
FortranArray2D<T> operator/(FortranArray2D<T>& arr, const T& scalar) { return arr.ref() / scalar; }

template<typename T>
FortranArray3D<T> operator*(FortranArray3DRef<T> arr, const T& scalar) {
    FortranArray3D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1),
                       arr.lbound(2), arr.ubound(2),
                       arr.lbound(3), arr.ubound(3));
    for (int k = arr.lbound(3); k <= arr.ubound(3); ++k)
        for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
            for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
                out(i, j, k) = arr(i, j, k) * scalar;
    return out;
}
template<typename T>
FortranArray3D<T> operator*(const T& scalar, FortranArray3DRef<T> arr) { return arr * scalar; }
template<typename T>
FortranArray3D<T> operator/(FortranArray3DRef<T> arr, const T& scalar) {
    FortranArray3D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1),
                       arr.lbound(2), arr.ubound(2),
                       arr.lbound(3), arr.ubound(3));
    for (int k = arr.lbound(3); k <= arr.ubound(3); ++k)
        for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
            for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
                out(i, j, k) = arr(i, j, k) / scalar;
    return out;
}
template<typename T>
FortranArray3D<T> operator*(FortranArray3D<T>& arr, const T& scalar) { return arr.ref() * scalar; }
template<typename T>
FortranArray3D<T> operator*(const T& scalar, FortranArray3D<T>& arr) { return arr.ref() * scalar; }
template<typename T>
FortranArray3D<T> operator/(FortranArray3D<T>& arr, const T& scalar) { return arr.ref() / scalar; }

template<typename T>
FortranArray4D<T> operator*(FortranArray4DRef<T> arr, const T& scalar) {
    FortranArray4D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1),
                       arr.lbound(2), arr.ubound(2),
                       arr.lbound(3), arr.ubound(3),
                       arr.lbound(4), arr.ubound(4));
    for (int l = arr.lbound(4); l <= arr.ubound(4); ++l)
        for (int k = arr.lbound(3); k <= arr.ubound(3); ++k)
            for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
                for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
                    out(i, j, k, l) = arr(i, j, k, l) * scalar;
    return out;
}
template<typename T>
FortranArray4D<T> operator*(const T& scalar, FortranArray4DRef<T> arr) { return arr * scalar; }
template<typename T>
FortranArray4D<T> operator/(FortranArray4DRef<T> arr, const T& scalar) {
    FortranArray4D<T> out;
    out.allocate_range(arr.lbound(1), arr.ubound(1),
                       arr.lbound(2), arr.ubound(2),
                       arr.lbound(3), arr.ubound(3),
                       arr.lbound(4), arr.ubound(4));
    for (int l = arr.lbound(4); l <= arr.ubound(4); ++l)
        for (int k = arr.lbound(3); k <= arr.ubound(3); ++k)
            for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
                for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
                    out(i, j, k, l) = arr(i, j, k, l) / scalar;
    return out;
}
template<typename T>
FortranArray4D<T> operator*(FortranArray4D<T>& arr, const T& scalar) { return arr.ref() * scalar; }
template<typename T>
FortranArray4D<T> operator*(const T& scalar, FortranArray4D<T>& arr) { return arr.ref() * scalar; }
template<typename T>
FortranArray4D<T> operator/(FortranArray4D<T>& arr, const T& scalar) { return arr.ref() / scalar; }

template<typename T>
std::vector<T> fortran_pack(FortranArray1DRef<T> arr) {
    std::vector<T> out(arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(arr.data(), arr.total_size(), out.data());
        return out;
    }
    size_t p = 0;
    for (int i = arr.lbound(1); i <= arr.ubound(1); ++i) out[p++] = arr(i);
    return out;
}

template<typename T>
std::vector<T> fortran_pack(FortranArray2DRef<T> arr) {
    std::vector<T> out(arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(arr.data(), arr.total_size(), out.data());
        return out;
    }
    size_t p = 0;
    for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
        for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
            out[p++] = arr(i, j);
    return out;
}

template<typename T>
std::vector<T> fortran_pack(FortranArray3DRef<T> arr) {
    std::vector<T> out(arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(arr.data(), arr.total_size(), out.data());
        return out;
    }
    size_t p = 0;
    for (int k = arr.lbound(3); k <= arr.ubound(3); ++k)
        for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
            for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
                out[p++] = arr(i, j, k);
    return out;
}

template<typename T>
std::vector<T> fortran_pack(FortranArray4DRef<T> arr) {
    std::vector<T> out(arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(arr.data(), arr.total_size(), out.data());
        return out;
    }
    size_t p = 0;
    for (int l = arr.lbound(4); l <= arr.ubound(4); ++l)
        for (int k = arr.lbound(3); k <= arr.ubound(3); ++k)
            for (int j = arr.lbound(2); j <= arr.ubound(2); ++j)
                for (int i = arr.lbound(1); i <= arr.ubound(1); ++i)
                    out[p++] = arr(i, j, k, l);
    return out;
}

template<typename T>
std::vector<T> fortran_pack(FortranArray1D<T>& arr) { return fortran_pack(arr.ref()); }
template<typename T>
std::vector<T> fortran_pack(FortranArray2D<T>& arr) { return fortran_pack(arr.ref()); }
template<typename T>
std::vector<T> fortran_pack(FortranArray3D<T>& arr) { return fortran_pack(arr.ref()); }
template<typename T>
std::vector<T> fortran_pack(FortranArray4D<T>& arr) { return fortran_pack(arr.ref()); }

template<typename T>
std::vector<T> fortran_pack(const FortranArray1D<T>& arr) {
    return std::vector<T>(arr.data(), arr.data() + arr.total_size());
}
template<typename T>
std::vector<T> fortran_pack(const FortranArray2D<T>& arr) {
    return std::vector<T>(arr.data(), arr.data() + arr.total_size());
}
template<typename T>
std::vector<T> fortran_pack(const FortranArray3D<T>& arr) {
    return std::vector<T>(arr.data(), arr.data() + arr.total_size());
}
template<typename T>
std::vector<T> fortran_pack(const FortranArray4D<T>& arr) {
    return std::vector<T>(arr.data(), arr.data() + arr.total_size());
}

template<typename T>
void fortran_unpack(FortranArray1DRef<T> arr, const T* src, size_t count) {
    const size_t n = std::min(count, arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(src, n, arr.data());
        return;
    }
    size_t p = 0;
    for (int i = arr.lbound(1); i <= arr.ubound(1) && p < n; ++i) arr(i) = src[p++];
}

template<typename T>
void fortran_unpack(FortranArray2DRef<T> arr, const T* src, size_t count) {
    const size_t n = std::min(count, arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(src, n, arr.data());
        return;
    }
    size_t p = 0;
    for (int j = arr.lbound(2); j <= arr.ubound(2) && p < n; ++j)
        for (int i = arr.lbound(1); i <= arr.ubound(1) && p < n; ++i)
            arr(i, j) = src[p++];
}

template<typename T>
void fortran_unpack(FortranArray3DRef<T> arr, const T* src, size_t count) {
    const size_t n = std::min(count, arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(src, n, arr.data());
        return;
    }
    size_t p = 0;
    for (int k = arr.lbound(3); k <= arr.ubound(3) && p < n; ++k)
        for (int j = arr.lbound(2); j <= arr.ubound(2) && p < n; ++j)
            for (int i = arr.lbound(1); i <= arr.ubound(1) && p < n; ++i)
                arr(i, j, k) = src[p++];
}

template<typename T>
void fortran_unpack(FortranArray4DRef<T> arr, const T* src, size_t count) {
    const size_t n = std::min(count, arr.total_size());
    if (arr.is_contiguous()) {
        std::copy_n(src, n, arr.data());
        return;
    }
    size_t p = 0;
    for (int l = arr.lbound(4); l <= arr.ubound(4) && p < n; ++l)
        for (int k = arr.lbound(3); k <= arr.ubound(3) && p < n; ++k)
            for (int j = arr.lbound(2); j <= arr.ubound(2) && p < n; ++j)
                for (int i = arr.lbound(1); i <= arr.ubound(1) && p < n; ++i)
                    arr(i, j, k, l) = src[p++];
}

template<typename T>
void fortran_unpack(FortranArray1DRef<T> arr, const std::vector<T>& src) {
    fortran_unpack(arr, src.data(), src.size());
}
template<typename T>
void fortran_unpack(FortranArray2DRef<T> arr, const std::vector<T>& src) {
    fortran_unpack(arr, src.data(), src.size());
}
template<typename T>
void fortran_unpack(FortranArray3DRef<T> arr, const std::vector<T>& src) {
    fortran_unpack(arr, src.data(), src.size());
}
template<typename T>
void fortran_unpack(FortranArray4DRef<T> arr, const std::vector<T>& src) {
    fortran_unpack(arr, src.data(), src.size());
}

template<typename T>
void fortran_unpack(FortranArray1D<T>& arr, const T* src, size_t count) { fortran_unpack(arr.ref(), src, count); }
template<typename T>
void fortran_unpack(FortranArray2D<T>& arr, const T* src, size_t count) { fortran_unpack(arr.ref(), src, count); }
template<typename T>
void fortran_unpack(FortranArray3D<T>& arr, const T* src, size_t count) { fortran_unpack(arr.ref(), src, count); }
template<typename T>
void fortran_unpack(FortranArray4D<T>& arr, const T* src, size_t count) { fortran_unpack(arr.ref(), src, count); }

#endif // FORTRAN_ARRAY_H
