/*
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * A CRTP base class for coordinate objects, providing partial specializations for 2D and 3D.
 */
#ifndef LSST_GEOM_COORDINATEBASE_H
#define LSST_GEOM_COORDINATEBASE_H

#include <iostream>
#include <type_traits>
#include <tuple>
#include <utility>

#include "Eigen/Core"

namespace lsst {
namespace geom {

template <typename T, int N = 2>
class Point;
template <typename T, int N = 2>
class Extent;

/// Test that a type is nothrow-copy-convertible from U to T.
template <typename T, typename U>
bool constexpr IS_NOTHROW_CONVERTIBLE =
        std::is_nothrow_copy_constructible<T>::value&& noexcept(static_cast<T>(std::declval<U>()));

/**
 *  A CRTP base class for coordinate objects.
 *
 *  CoordinateBase has partial specializations for 2 and 3 dimensions so its subclasses don't have to.
 */
template <typename Derived, typename T, int N>
class CoordinateBase {
public:
    static_assert(N > 0, "CoordinateBase must have a positive length.");
    typedef T Element;
    static int const dimensions = N;
    typedef Eigen::Matrix<T, N, 1, Eigen::DontAlign> EigenVector;
    static bool constexpr IS_ELEMENT_NOTHROW_COPYABLE = std::is_nothrow_copy_constructible<T>::value;
    static bool constexpr IS_ELEMENT_NOTHROW_ASSIGNABLE = std::is_nothrow_copy_assignable<T>::value;

    // Can't use both =default and noexcept until Eigen supports noexcept
    CoordinateBase(CoordinateBase const& other) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(other._vector) {}
    CoordinateBase(CoordinateBase&& other) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(std::move(other._vector)) {}
    CoordinateBase& operator=(CoordinateBase const& other) noexcept(IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        _vector = other._vector;
        return *this;
    }
    CoordinateBase& operator=(CoordinateBase&& other) noexcept(IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        _vector = std::move(other._vector);
        return *this;
    }
    ~CoordinateBase() noexcept = default;

    T& operator[](int n) { return _vector[n]; }
    T const& operator[](int n) const { return const_cast<EigenVector&>(_vector)[n]; }
    T& coeffRef(int n) { return _vector.coeffRef(n); }
    T const& coeffRef(int n) const { return const_cast<EigenVector&>(_vector).coeffRef(n); }

    /**
     *  Return a fixed-size Eigen representation of the coordinate object.
     *
     *  The fact that this returns by const reference rather than by value should not be considered
     *  part of the API; this is merely an optimization enabled by the implementation.
     */
    EigenVector const& asEigen() const noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { return _vector; }

protected:
    /**
     *  Initialize all elements to a scalar.
     *
     *  A public constructor with the same signature is expected for subclasses.
     */
    explicit CoordinateBase(T val = static_cast<T>(0)) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(EigenVector::Constant(val)) {}

    /**
     *  Initialize all elements from an N-d Eigen vector.
     *
     *  A public constructor with the same signature is expected for subclasses.
     */
    template <typename Vector>
    explicit CoordinateBase(Eigen::MatrixBase<Vector> const& vector) : _vector(vector) {}

    void _swap(CoordinateBase& other) noexcept { _vector.swap(other._vector); }
    EigenVector _vector;
};

/**
 *  Floating-point comparison with tolerance.
 *
 *  Interface, naming, and default tolerances matches Numpy.
 *
 *  @relatesalso CoordinateBase
 */
template <typename Derived, typename T, int N>
bool allclose(CoordinateBase<Derived, T, N> const& a, CoordinateBase<Derived, T, N> const& b,
              T rtol = static_cast<T>(1E-5),
              T atol = static_cast<T>(1E-8)) noexcept(std::is_nothrow_copy_constructible<T>::value&&
                                                              std::is_nothrow_copy_assignable<T>::value);

/**
 *  Specialization of CoordinateBase for 2 dimensions.
 */
template <typename Derived, typename T>
class CoordinateBase<Derived, T, 2> {
public:
    typedef T Element;
    static int const dimensions = 2;
    typedef Eigen::Matrix<T, 2, 1, Eigen::DontAlign> EigenVector;
    static bool constexpr IS_ELEMENT_NOTHROW_COPYABLE = std::is_nothrow_copy_constructible<T>::value;
    static bool constexpr IS_ELEMENT_NOTHROW_ASSIGNABLE = std::is_nothrow_copy_assignable<T>::value;

    // Can't use both =default and noexcept until Eigen supports noexcept
    CoordinateBase(CoordinateBase const& other) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(other._vector) {}
    CoordinateBase(CoordinateBase&& other) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(std::move(other._vector)) {}
    CoordinateBase& operator=(CoordinateBase const& other) noexcept(IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        _vector = other._vector;
        return *this;
    }
    CoordinateBase& operator=(CoordinateBase&& other) noexcept(IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        _vector = std::move(other._vector);
        return *this;
    }
    ~CoordinateBase() noexcept = default;

    T& operator[](int n) { return _vector[n]; }
    T const& operator[](int n) const { return const_cast<EigenVector&>(_vector)[n]; }
    T& coeffRef(int n) { return _vector.coeffRef(n); }
    T const& coeffRef(int n) const { return const_cast<EigenVector&>(_vector).coeffRef(n); }

    /**
     *  Return a fixed-size Eigen representation of the coordinate object.
     *
     *  The fact that this returns by const reference rather than by value should not be considered
     *  part of the API; this is merely an optimization enabled by the implementation.
     */
    EigenVector const& asEigen() const noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { return _vector; }

    T const& getX() const noexcept { return _vector.x(); }
    T const& getY() const noexcept { return _vector.y(); }
    T& getX() noexcept { return _vector.x(); }
    T& getY() noexcept { return _vector.y(); }
    void setX(T x) noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { _vector.x() = x; }
    void setY(T y) noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { _vector.y() = y; }

    /// Return a std::pair representation of the coordinate object.
    std::pair<T, T> asPair() const noexcept(IS_ELEMENT_NOTHROW_COPYABLE) {
        return std::make_pair(_vector.x(), _vector.y());
    }

    /// Return a std::tuple representation of the coordinate object.
    std::tuple<T, T> asTuple() const noexcept(IS_ELEMENT_NOTHROW_COPYABLE) {
        return std::make_tuple(_vector.x(), _vector.y());
    }

protected:
    explicit CoordinateBase(T val = static_cast<T>(0)) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(EigenVector::Constant(val)) {}

    template <typename Vector>
    explicit CoordinateBase(Eigen::MatrixBase<Vector> const& vector) : _vector(vector) {}
    void _swap(CoordinateBase& other) noexcept { _vector.swap(other._vector); }
    EigenVector _vector;
};

/**
 *  Specialization of CoordinateBase for 3 dimensions.
 */
template <typename Derived, typename T>
class CoordinateBase<Derived, T, 3> {
public:
    typedef T Element;
    static int const dimensions = 3;
    typedef Eigen::Matrix<T, 3, 1, Eigen::DontAlign> EigenVector;
    static bool constexpr IS_ELEMENT_NOTHROW_COPYABLE = std::is_nothrow_copy_constructible<T>::value;
    static bool constexpr IS_ELEMENT_NOTHROW_ASSIGNABLE = std::is_nothrow_copy_assignable<T>::value;

    // Can't use both =default and noexcept until Eigen supports noexcept
    CoordinateBase(CoordinateBase const& other) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(other._vector) {}
    CoordinateBase(CoordinateBase&& other) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(std::move(other._vector)) {}
    CoordinateBase& operator=(CoordinateBase const& other) noexcept(IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        _vector = other._vector;
        return *this;
    }
    CoordinateBase& operator=(CoordinateBase&& other) noexcept(IS_ELEMENT_NOTHROW_ASSIGNABLE) {
        _vector = std::move(other._vector);
        return *this;
    }
    ~CoordinateBase() noexcept = default;

    T& operator[](int n) { return _vector[n]; }
    T const& operator[](int n) const { return const_cast<EigenVector&>(_vector)[n]; }
    T& coeffRef(int n) { return _vector.coeffRef(n); }
    T const& coeffRef(int n) const { return const_cast<EigenVector&>(_vector).coeffRef(n); }

    /**
     *  Return a fixed-size Eigen representation of the coordinate object.
     *
     *  The fact that this returns by const reference rather than by value should not be considered
     *  part of the API; this is merely an optimization enabled by the implementation.
     */
    EigenVector const& asEigen() const noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { return _vector; }

    T const& getX() const noexcept { return _vector.x(); }
    T const& getY() const noexcept { return _vector.y(); }
    T const& getZ() const noexcept { return _vector.z(); }
    T& getX() noexcept { return _vector.x(); }
    T& getY() noexcept { return _vector.y(); }
    T& getZ() noexcept { return _vector.z(); }
    void setX(T x) noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { _vector.x() = x; }
    void setY(T y) noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { _vector.y() = y; }
    void setZ(T z) noexcept(IS_ELEMENT_NOTHROW_COPYABLE) { _vector.z() = z; }

    /// Return a std::tuple representation of the coordinate object.
    std::tuple<T, T, T> asTuple() const noexcept(IS_ELEMENT_NOTHROW_COPYABLE) {
        return std::make_tuple(_vector.x(), _vector.y(), _vector.z());
    }

protected:
    explicit CoordinateBase(T val = static_cast<T>(0)) noexcept(IS_ELEMENT_NOTHROW_COPYABLE)
            : _vector(EigenVector::Constant(val)) {}

    template <typename Vector>
    explicit CoordinateBase(Eigen::MatrixBase<Vector> const& vector) : _vector(vector) {}
    void _swap(CoordinateBase& other) noexcept { _vector.swap(other._vector); }
    EigenVector _vector;
};

template <typename Derived, typename T, int N>
std::ostream& operator<<(std::ostream& os, CoordinateBase<Derived, T, N> const& coordinate) {
    os << "(" << coordinate[0];
    for (int n = 1; n < N; ++n) os << ", " << coordinate[n];
    return os << ")";
}

}  // namespace geom
}  // namespace lsst

#endif
