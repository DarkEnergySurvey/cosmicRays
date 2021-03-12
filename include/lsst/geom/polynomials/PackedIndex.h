// -*- LSST-C++ -*-
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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PackedIndex_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PackedIndex_h_INCLUDED

namespace lsst { namespace geom { namespace polynomials {

/// Enum defining the packing orders used to order 2-d polynomial coefficients.
enum class PackingOrder {

    /**
     *  A pair of indices @f$(nx, ny)@f$ is mapped to the flattened position
     *  @f$i = (nx+ny)(nx+ny+1)/2 + nx@f$, which yields the (nx, ny) ordering
     *  ```
     *      (0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0), ...
     *  ```
     *
     *  In other words, if @f$a_i@f$ are the coefficients of a 2-d polynomial
     *  expansion with order=2, the full polynomial is:
     *  @f[
     *     a_0 + a_1 y + a_2 x + a_3 y^2 + a_4 x y + a_5 x^2
     *  @f]
     */
    YX,

    /**
     *  A pair of indices @f$(nx, ny)@f$ is mapped to the flattened position
     *  @f$i = (nx+ny)(nx+ny+1)/2 + ny@f$, which yields the (nx, ny) ordering
     *  ```
     *      (0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), ...
     *  ```
     *
     *  In other words, if @f$a_i@f$ are the coefficients of a 2-d polynomial
     *  expansion with order=2, the full polynomial is:
     *  @f[
     *     a_0 + a_1 x + a_2 y + a_3 x^2 + a_4 x y + a_5 y^2
     *  @f]
     */
    XY
};

/**
 *  A custom tuple that relates the indices of two 1-d functions for x and y
 *  to the flattened index for the 2-d function they form.
 *
 *  The packing algorithm is not specified by Index2d itself; it is intended
 *  to be a common data structure for other classes that define possibly
 *  different algorithms (e.g. PackedIndexIterator).
 */
struct Index2d {

    /// Construct an index with zero entries.
    constexpr Index2d() noexcept : flat(0), nx(0), ny(0) {}

    /// Construct with the provided values.
    constexpr Index2d(std::size_t flat_, std::size_t nx_, std::size_t ny_) noexcept :
        flat(flat_), nx(nx_), ny(ny_)
    {}

    /// Equality comparison.
    constexpr bool operator==(Index2d const & other) const noexcept {
        return flat == other.flat && nx == other.nx && ny == other.ny;
    }

    /// Inequality comparison.
    constexpr bool operator!=(Index2d const & other) const noexcept {
        return !(*this == other);
    }

    std::size_t flat; ///< Index into the flattened 2-d function.
    std::size_t nx;    ///< Index into the 1-d function for nx.
    std::size_t ny;    ///< Index into the 1-d functoin for ny.
};

namespace detail {

// Specialization of all PackedIndexIterator/PackedIndexRange logic
// that isn't common across PackingOrders.
template <PackingOrder packing> struct PackingOrderTraits;

template <>
struct PackingOrderTraits<PackingOrder::YX> {

    // Return the offset of the given coefficient after the (nx + ny) offset is subtracted off.
    static std::size_t computeInnerIndex(std::size_t nx, std::size_t ny) { return nx; }

    static void increment(Index2d & index) {
        if (index.ny == 0) {
            index.ny = index.nx + 1;
            index.nx = 0;
        } else {
            --index.ny;
            ++index.nx;
        }
    }

    // Return the nx and ny values appropriate for the end iterator of a range of the given order.
    static std::size_t getEndX(std::size_t order) { return 0; }
    static std::size_t getEndY(std::size_t order) { return order + 1; }

};

template <>
struct PackingOrderTraits<PackingOrder::XY> {

    // Return the offset of the given coefficient after the (nx + ny) offset is subtracted off.
    static std::size_t computeInnerIndex(std::size_t nx, std::size_t ny) { return ny; }

    static void increment(Index2d & index) {
        if (index.nx == 0) {
            index.nx = index.ny + 1;
            index.ny = 0;
        } else {
            --index.nx;
            ++index.ny;
        }
    }

    // Return the nx and ny values appropriate for the end iterator of a range of the given order.
    static std::size_t getEndX(std::size_t order) { return order + 1; }
    static std::size_t getEndY(std::size_t order) { return 0; }

};

} // namespace detail

/**
 *  An iterator for traversing "packed" triangular 2-d series expansions,
 *  in which two 1-d expansions are truncated according to the sum of their
 *  orders and all values for one order are stored before any values of the
 *  subsequent order.
 *
 *  PackedIndexIterator dereferences to Index2d.  Typical usage is via a
 *  PackedIndexRange.
 *
 *  This packing ensures that the coefficients for an nth-order expansion are
 *  a contiguous subset of the coefficients for an (n+1)th-order expansion.
 *
 *  The packing within each order is set by the `packing` template parameter.
 *
 *  PackedIndexIterator is an STL input iterator, *not* a forward iterator;
 *  incrementing the iterator invalidates all previously-dereferenced values.
 */
template <PackingOrder packing>
class PackedIndexIterator {
    using Traits = detail::PackingOrderTraits<packing>;
public:

    using difference_type = std::ptrdiff_t;
    using value_type = Index2d;
    using pointer = Index2d const *;
    using reference = Index2d const &;
    using iterator_category = std::input_iterator_tag;

    /// Return the flattened offset to the start of the given order.
    static constexpr std::size_t computeOffset(std::size_t order) noexcept {
        return order*(order + 1)/2;
    }

    /// Return the flattened size of an expansion with the given maximum order (inclusive).
    static constexpr std::size_t computeSize(std::size_t order) noexcept {
        return computeOffset(order + 1);
    }

    /// Return the flattened index for the element with the given x and y orders.
    static constexpr std::size_t computeIndex(std::size_t nx, std::size_t ny) noexcept {
        return computeOffset(nx + ny) + Traits::computeInnerIndex(nx, ny);
    }

    /// Construct an iterator one past the end of an expansion with the given order.
    static constexpr PackedIndexIterator makeEnd(std::size_t order) noexcept {
        return PackedIndexIterator(order);
    }

    /// Construct an iterator at the beginning of an expansion of any order.
    constexpr PackedIndexIterator() noexcept : _index() {}

    /// Construct an iterator pointing to the element with the given x and y orders.
    constexpr PackedIndexIterator(std::size_t nx, std::size_t ny) noexcept :
        _index(computeIndex(nx, ny), nx, ny)
    {}

    /// Dereference the iterator, yielding a Index2d const reference.
    constexpr reference operator*() const noexcept { return _index; }

    /// Dereference the iterator, yielding a Index2d const pointer.
    constexpr pointer operator->() const noexcept { return &_index; }

    /// Move to the next element in the packed array and return the iterator.
    PackedIndexIterator & operator++() noexcept {
        ++_index.flat;
        Traits::increment(_index);
        return *this;
    }

    /// Move to the next element in the packed array and return a copy of the iterator before the move.
    PackedIndexIterator operator++(int) noexcept {
        PackedIndexIterator r(*this);
        ++(*this);
        return r;
    }

    /// Equality comparison.
    constexpr bool operator==(PackedIndexIterator const & other) const noexcept {
        return _index == other._index;
    }

    /// Inequality comparison.
    constexpr bool operator!=(PackedIndexIterator const & other) const noexcept {
        return !(*this == other);
    }

private:

    constexpr PackedIndexIterator(std::size_t order) noexcept :
        _index(computeOffset(order + 1), Traits::getEndX(order), Traits::getEndY(order))
    {}

    Index2d _index;
};

/**
 *  A specialized iterator range class for PackedIndexIterator, providing
 *  size calculation, comparison, and range-based `for` support.
 *
 *  @see PackedIndexIterator for information on the packing algorithm.
 */
template <PackingOrder packing>
class PackedIndexRange {
public:

    using iterator = PackedIndexIterator<packing>;
    using const_iterator = iterator;
    using value_type = typename iterator::value_type;
    using reference = typename iterator::reference;
    using pointer = typename iterator::pointer;
    using difference_type = typename iterator::difference_type;
    using size_type = std::size_t;

    /// Return the flattened offset to the start of the given order.
    static constexpr std::size_t computeOffset(std::size_t order) noexcept {
        return iterator::computeOffset(order);
    }

    /// Return the flattened size of an expansion with the given maximum order (inclusive).
    static constexpr std::size_t computeSize(std::size_t order) noexcept {
        return iterator::computeSize(order);
    }

    /// Return the flattened index for the element with the given x and y orders.
    static constexpr std::size_t computeIndex(std::size_t nx, std::size_t ny) noexcept {
        return iterator::computeIndex(nx, ny);
    }

    /// Construct from begin and end iterators.
    constexpr PackedIndexRange(iterator first, iterator last) noexcept :
        _begin(first),
        _end(last)
    {}

    /// Return an iterator to the start of the range.
    constexpr iterator begin() const noexcept { return _begin; }

    /// Return an iterator to the start of the range.
    constexpr iterator cbegin() const noexcept { return _begin; }

    /// Return an iterator to one past the end of the range.
    constexpr iterator end() const noexcept { return _end; }

    /// Return an iterator to one past the end of the range.
    constexpr iterator cend() const noexcept { return _end; }

    /// Return the number of elements in the flattened expansion.
    constexpr std::size_t size() const noexcept { return _end->flat - _begin->flat; }

    /// Return true if the number of elements in the flattened expansion is zero.
    constexpr bool empty() const noexcept { return size() == 0u; }

    /// Equality comparison.
    constexpr bool operator==(PackedIndexRange const & other) const noexcept {
        return _begin == other._begin && _end == other._end;
    }

    /// Inequality comparison.
    constexpr bool operator!=(PackedIndexRange const & other) const noexcept {
        return !(*this == other);
    }

private:
    iterator _begin;
    iterator _end;
};

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PackedIndex_h_INCLUDED
