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

#ifndef LSST_GEOM_INTERVAL_H
#define LSST_GEOM_INTERVAL_H

#include <vector>

#include "ndarray.h"

#include "lsst/geom/Point.h"
#include "lsst/geom/Extent.h"

namespace lsst {
namespace geom {

class IntervalD;

/**
 *  A 1-d integer coordinate range.
 *
 *  IntervalI is an inclusive range that represents region of pixels.  An
 *  IntervalI never has negative size; the empty interval is defined to have
 *  zero size, and is treated as though it does not have a well-defined
 *  position (regardless of the return value of getMin() or getMax() for an
 *  empty interval).
 *
 *  All IntervalI methods that return a new instance (and are not marked
 *  `noexcept`) throw `OverflowError` if the lower bound, upper bound, or the
 *  size would be too large to fit in `int`.
 */
class IntervalI final {
public:
    using Element = int;

    /**
     *  Enum used to indicate how to handle conversions from floating-point to
     *  integer intervals.
     *
     *  Note that the floating point bounds for some pixel `x` in an
     *  integer interval are `[x - 0.5, x + 0.5]` (because the pixel has
     *  unit size, and integer coordinates correspond to the centers of
     *  pixels).  This means interval conversions involve more than just
     *  rounding the bounds inward or outward.
     */
    enum class EdgeHandlingEnum {
        /**
         *  Include all pixels that overlap the floating-point interval at all.
         */
        EXPAND,
        /**
         *  Include only pixels that are wholly contained by the floating-point
         *  interval.
         */
        SHRINK
    };

    /// Construct an empty interval.
    IntervalI() noexcept : _min(0), _size(0) {}

    /**
     *  Construct an interval that contains all of the given points.
     *
     *  @param[in] first  Iterator to the beginning of a sequence of integers.
     *  @param[in] last   Iterator to one-past-the-end of a sequence of
     *                    integers.
     *
     *  An empty interval is returned if `first == last`.
     *
     *  @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *      interval would overflow.
     */
    template <typename Iter>
    static IntervalI fromSpannedPoints(Iter first, Iter last) {
        IntervalI result;
        for (auto i = first; i != last; ++i) {
            result = result.expandedTo(*i);
        }
        return result;
    }

    //@{
    /**
     *  Construct an interval that contains all of the given points.
     *
     *  @param[in] elements  Points (integer values) to include in the
     *                       interval.
     *
     *  An empty interval is returned if the given container has no elements.
     *
     *  @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *      interval would overflow.
     */
    static IntervalI fromSpannedPoints(std::vector<Element> const& elements) {
        return fromSpannedPoints(elements.begin(), elements.end());
    }
    static IntervalI fromSpannedPoints(ndarray::Array<Element const, 1> const& elements) {
        return fromSpannedPoints(elements.begin(), elements.end());
    }
    //@}

    /**
     *  Construct an interval from its lower and upper bounds.
     *
     *  If min > max, the resulting interval is empty.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] max       Maximum coordinate (inclusive).
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     */
    static IntervalI fromMinMax(Element min, Element max);

    /**
     *  Construct an interval from its lower bound and size.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] size      Number of pixels in interval.  Nonpositive values
     *                       will produce an empty interval.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     */
    static IntervalI fromMinSize(Element min, Element size);

    /**
     *  Construct an interval from its upper bound and size.
     *
     *  @param[in] max       Maximum coordinate (inclusive).
     *  @param[in] size      Number of pixels in interval.  Nonpositive values
     *                       will produce an empty interval.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     */
    static IntervalI fromMaxSize(Element max, Element size);

    /**
     * Create an interval centered as closely as possible on a particular
     * point.
     *
     * @param center The desired center of the interval.
     * @param size   Number of pixels in interval.  Nonpositive values
     *               will produce an empty interval.
     *
     * @returns if `size` is positive, an interval with size `size`; if `size`
     *          is zero or negative, an empty interval. If the returned
     *          interval is not empty, its center shall be within half a pixel
     *          of `center`.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if `center`
     *     is not finite and size is positive.
     */
    static IntervalI fromCenterSize(double center, Element size);

    /**
     *  Construct an integer interval from a floating-point interval.
     *
     *  Floating-point to integer interval conversion is based on the concept
     *  that a pixel is not an infinitesimal point but rather a square of unit
     *  size centered on integer-valued coordinates.  Converting a
     *  floating-point interval to an integer interval thus requires a choice
     *  on how to handle pixels which are only partially contained by the
     *  input floating-point interval.
     *
     *  @param[in] other          A floating-point interval to convert.
     *  @param[in] edgeHandling   If EXPAND, the integer interval will contain
     *                            any pixels that overlap the floating-point
     *                            interval.  If SHRINK, the integer interval
     *                            will contain only pixels completely
     *                            contained by the floating-point interval.
     *
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if `other`
     *     is not finite.
     * @throws lsst::pex::exceptions::LogicError Thrown if an invalid enum
     *     value is passed.
     */
    explicit IntervalI(IntervalD const& other, EdgeHandlingEnum edgeHandling = EdgeHandlingEnum::EXPAND);

    /// Standard copy constructor.
    IntervalI(IntervalI const&) noexcept = default;

    /// Standard move constructor.
    IntervalI(IntervalI&&) noexcept = default;

    ~IntervalI() noexcept = default;

    void swap(IntervalI& other) noexcept {
        using std::swap;
        swap(_size, other._size);
        swap(_min, other._min);
    }

    /// Standard copy assignment operator.
    IntervalI& operator=(IntervalI const&) noexcept = default;

    /// Standard move assignment operatior.
    IntervalI& operator=(IntervalI&&) noexcept = default;

    /**
     *  @name Min/Max Accessors
     *
     *  Return the minimum and maximum coordinates of the interval (inclusive).
     */
    //@{
    Element getMin() const noexcept { return _min; }
    Element getMax() const noexcept { return _min + _size - 1; }
    //@}

    /**
     *  @name Begin/End Accessors
     *
     *  Return begin (inclusive) and end (exclusive) coordinates for the
     *  interval.
     */
    //@{
    Element getBegin() const noexcept { return _min; }
    Element getEnd() const noexcept { return _min + _size; }
    //@}

    /**
     *  @name Size Accessors
     *
     *  Return the size of the interval in pixels.
     */
    Element getSize() const noexcept { return _size; }

    /**
     *  Return slice to extract the interval's region from an ndarray::Array.
     *
     *  For example:
     *
     *      auto array = ndarray::copy(ndarray::arange(5);
     *      auto interval = IntervalI::fromMinMax(2, 4);
     *      auto subarray = array[interval.getSlice()];
     *
     */
    ndarray::View<boost::fusion::vector1<ndarray::index::Range> > getSlice() const;

    /// Return true if the interval contains no points.
    bool isEmpty() const noexcept { return _size == 0; }

    /// Return true if the interval contains the point.
    bool contains(Element point) const noexcept;

    /**
     *  Return true if all points contained by other are also contained by
     *  this.
     *
     *  An empty interval is contained by every other interval, including
     *  other empty intervals.
     */
    bool contains(IntervalI const& other) const noexcept;

    //@{
    /**
     *  Return true if there are any points in both this and other.
     *
     *  Any overlap operation involving an empty interval returns false.
     */
    bool overlaps(IntervalI const& other) const noexcept;
    bool intersects(IntervalI const& other) const noexcept { return overlaps(other); }
    //@}

    /**
     *  Return true if there are no points in both `this` and `other`.
     */
    bool isDisjointFrom(IntervalI const& other) const noexcept;

    /**
     *  Increase the size of the interval by the given amount in both
     *  directions (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to eroding by `-buffer`.
     *
     *  If the final size of the interval is less than or equal to zero the
     *  new interval will be empty.
     *
     *  Empty intervals remain empty after dilation.
     */
    IntervalI dilatedBy(Element buffer) const;

    /**
     *  Decrease the size of the interval by the given amount in both
     *  directions (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to dilating by `-buffer`.
     *
     *  If the final size of the interval is less than or equal to zero the
     *  new interval will be empty.
     *
     *  Empty intervals remain empty after erosion.
     */
    IntervalI erodedBy(Element buffer) const { return dilatedBy(-buffer); }

    /**
     * Shift the position of the interval by the given offset (returning a new
     * object).)
     *
     * Empty intervals remain empty when shifted.
     */
    IntervalI shiftedBy(Element offset) const;

    /**
     * Reflect an interval about a point (returning a new object).
     *
     * Empty intervals remain empty when reflected.
     */
    IntervalI reflectedAbout(Element point) const;

    //@{
    /**
     * Expand an interval to ensure that `contains(other)` is true (returning a
     * new object).
     *
     * Expanding an empty interval with a single point yields an interval with
     * size=1 at that point; expanding an empty interval with a second interval
     * is equivalent to assignment.
     */
    IntervalI expandedTo(Element other) const;
    IntervalI expandedTo(IntervalI const& other) const;
    //@}

    /**
     * Shrink an interval to ensure that it is contained by other (returning a
     * new)
     *
     * In particular, if `other` and this interval do not overlap, the returned
     * interval will be empty.
     */
    IntervalI clippedTo(IntervalI const& other) const noexcept;

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator==(IntervalI const& other) const noexcept;

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator!=(IntervalI const& other) const noexcept;

    /// Return a hash of this object.
    std::size_t hash_value() const noexcept;

    std::string toString() const;

private:
    template <typename T>
    static IntervalI _fromMinMaxChecked(T min, T max);

    IntervalI(Element min, Element max);

    /*
     *  IntervalI internally stores its minimum point and size, because we
     *  expect these will be the most commonly accessed quantities.
     *
     *  We set the minimum point to the origin for an empty interval, and
     *  use -1 for the maximum point in that case, but this is an internal
     *  detail - all the API guarantees is that for an empty interval,
     *  size == 0.
     */

    Element _min;
    Element _size;
};

/**
 *  A floating-point coordinate rectangle geometry.
 *
 *  IntervalD is closed (its bounds are considered included in the interval).
 *  An interval never has negative size; the empty interval is defined to
 *  zero-size size and its minimum and maximum values are set to NaN.
 *  Non-empty intervals representing infinitesimal points may also have zero
 *  size, but are not considered empty.
 *
 *  The existence of zero-size, non-empty intervals is an important,
 *  intentional difference between IntervalI and IntervalD, related to the fact
 *  that IntervalI models a discrete set while IntervalD (imperfectly, due to
 *  floating-point limitations) models a continuous set.
 *
 *  @internal
 *
 *  IntervalD internally stores its minimum point and maximum point, instead of
 *  minimum point and size, to ensure roundoff error does not affect whether
 *  points are contained by the interval.
 */
class IntervalD final {
public:
    using Element = double;

    /// Construct an empty interval.
    IntervalD() noexcept;

    /**
     *  Construct an interval that contains all of the given points.
     *
     *  @param[in] first  Iterator to the beginning of a sequence of
     *                    floating-point values.
     *  @param[in] last   Iterator to one-past-the-end of a sequence of
     *                    floating-point values.
     *
     *  An empty interval is returned if `first == last`.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError Thrown if any
     *      input point is non-finite.
     */
    template <typename Iter>
    static IntervalD fromSpannedPoints(Iter first, Iter last) {
        IntervalD result;
        for (auto i = first; i != last; ++i) {
            result = result.expandedTo(*i);
        }
        return result;
    }

    //@{
    /**
     *  Construct an interval that contains all of the given points.
     *
     *  @param[in] elements  Points (floating-point values) to include in the
     *                       interval.
     *
     *  An empty interval is returned if the given container has no elements.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError Thrown if any
     *      input point is non-finite.
     */
    static IntervalD fromSpannedPoints(std::vector<Element> const& elements) {
        return fromSpannedPoints(elements.begin(), elements.end());
    }
    static IntervalD fromSpannedPoints(ndarray::Array<Element const, 1> const& elements) {
        return fromSpannedPoints(elements.begin(), elements.end());
    }
    //@}

    /**
     *  Construct an interval from its lower and upper bounds.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] max       Maximum coordinate (inclusive).
     *
     *  Bounds may be non-finite, with some restrictions.
     *   - If `min` is `-inf` and/or `max` is `+inf`, the interval has infinite
     *     size and is considered valid.
     *   - If `max < min` (regardless of whether one or both is infinite) or
     *     either is `NaN`, the interval is empty.
     *   - `min == max == +inf` or `min == max == -inf` is an error.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *     `min == max == inf` or `min == max == -inf`.
     */
    static IntervalD fromMinMax(Element min, Element max);

    /**
     *  Construct an interval from its lower bound and size.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] size      Size of interval.
     *
     *  If `size < 0` or any parameter is NaN, an empty interval is returned.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if `min`
     *       is infinite or if size is +infinity (-infinity yields an empty
     *       interval).
     */
    static IntervalD fromMinSize(Element min, Element size);

    /**
     *  Construct an interval from its upper bound and size.
     *
     *  @param[in] max       Maximum coordinate (inclusive).
     *  @param[in] size      Size of interval.
     *
     *  If `size < 0` or any parameter is NaN, an empty interval is returned.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if `min`
     *       is infinite or if size is +infinity (-infinity yields an empty
     *       interval).
     */
    static IntervalD fromMaxSize(Element max, Element size);

    /**
     *  Construct an interval centered on a particular point.
     *
     *  @param center The desired center of the interval.  May not be infinite.
     *  @param size   Number of pixels in interval.  May not be infinite.
     *
     *  If `size < 0` or any parameter is NaN, an empty interval is returned.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError Thrown if any
     *      parameter is infinite.  This includes the case where one parameter
     *      is NaN and the other is infinite.
     */
    static IntervalD fromCenterSize(double center, Element size);

    /**
     *  Construct a floating-point interval from an integer interval.
     *
     *  Integer to floating-point interval conversion is based on the concept
     *  that a pixel is not an infinitesimal point but rather a square of unit
     *  size centered on integer-valued coordinates.  While the output
     *  floating-point interval thus has the same size as the input integer
     *  interval, its minimum/maximum coordinates are 0.5 smaller/greater.
     */
    explicit IntervalD(IntervalI const& other) noexcept;

    /// Standard copy constructor.
    IntervalD(IntervalD const&) noexcept = default;

    /// Standard move constructor.
    IntervalD(IntervalD&&) noexcept = default;

    ~IntervalD() noexcept = default;

    void swap(IntervalD& other) noexcept {
        std::swap(_min, other._min);
        std::swap(_max, other._max);
    }

    /// Standard copy assignment operator.
    IntervalD& operator=(IntervalD const&) noexcept = default;

    /// Standard move assignment operator.
    IntervalD& operator=(IntervalD&&) noexcept = default;

    /**
     *  @name Min/Max Accessors
     *
     *  Return the minimum (inclusive) and maximum (exclusive) coordinates of
     *  the interval.
     */
    //@{
    Element getMin() const noexcept { return _min; }
    Element getMax() const noexcept { return _max; }
    //@}

    /**
     *  Return the size of the interval.
     *
     *  Empty intervals have zero size, but not all zero-size intervals are
     *  empty.  Intervals with an infinite bound have infinite size.
     */
    Element getSize() const noexcept;

    /**
     *  Return the center coordinate of the interval.
     *
     *  Returns NaN for empty intervals and infinite intervals.
     */
    Element getCenter() const noexcept;

    /// Return true if the interval contains no points.
    bool isEmpty() const noexcept { return std::isnan(_min); }

    /// Return true if the interval's size is finite.
    bool isFinite() const noexcept { return std::isfinite(getSize()); }

    /**
     *  Return true if the interval contains the point.
     *
     *   @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `point` is NaN.
     */
    bool contains(Element point) const;

    /**
     *  Return true if all points contained by other are also contained by
     *  this.
     *
     *  An empty interval is contained by every other interval, including
     *  other empty intervals.
     */
    bool contains(IntervalD const& other) const noexcept;

    //@{
    /**
     *  Return true if any points in other are also in this.
     *
     *  Any overlap operation involving an empty interval returns false.
     */
    bool overlaps(IntervalD const& other) const noexcept;
    bool intersects(IntervalD const& other) const noexcept { return overlaps(other); }
    //@}

    /**
     *  Return true if there are no points in both `this` and `other`.
     */
    bool isDisjointFrom(IntervalD const& other) const noexcept;

    /**
     *  Increase the size of the interval by the given amount in both
     *  directions (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to eroding by `-buffer`.
     *
     *  If the final size of the interval is less than zero, the new interval
     *  will be empty.
     *
     *  Empty intervals remain empty after dilation.  Infinite bounds are
     *  unaffected by dilation.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `buffer` is not finite.
     */
    IntervalD dilatedBy(Element buffer) const;

    /**
     *  Decrease the size of the interval by the given amount in both
     *  directions (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to dilating by `-buffer`.
     *
     *  If the final size of the interval is less than zero, the new interval
     *  will be empty.
     *
     *  Empty intervals remain empty after erosion.  Infinite bounds are
     *  unaffected by erosion.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `buffer` is not finite.
     */
    IntervalD erodedBy(Element buffer) const { return dilatedBy(-buffer); }

    /**
     *  Shift the position of the interval by the given offset (returning a
     *  new object).
     *
     *  Empty intervals remain empty when shifted.  Infinite bounds are
     *  unaffected by shifting.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `offset` is not finite.
     */
    IntervalD shiftedBy(Element offset) const;

    /**
     *  Reflect an interval about a point (returning a new object).
     *
     *  Empty intervals remain empty when reflected.  If an interval with an
     *  infinite bound is reflected, the opposite bound will become infinite
     *  with the opposite sign.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `point` is not finite.
     */
    IntervalD reflectedAbout(Element point) const;

    /**
     *  Expand an interval to ensure that `contains(other)` is true.
     *
     *  Expanding an empty interval with a single point yields an interval
     *  with `size == 0` at that point.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if `other`
     *       is not finite.
     */
    IntervalD expandedTo(Element other) const;

    /**
     *  Expand an interval to ensure that `contains(other)` is true.
     *
     *  Expanding an empty interval with a second interval is equivalent to
     *  assignment.
     *
     *  Expanding by an empty interval yields the original interval.
     */
    IntervalD expandedTo(IntervalD const& other) const noexcept;

    /**
     * Shrink an interval to ensure that it is contained by other (returning a
     * new object).
     *
     * In particular, if `other` and this interval do not overlap the new
     * interval will be empty.
     */
    IntervalD clippedTo(IntervalD const& other) const noexcept;

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator==(IntervalD const& other) const noexcept;

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator!=(IntervalD const& other) const noexcept;

    /// Return a hash of this object.
    std::size_t hash_value() const noexcept;

    std::string toString() const;

private:
    IntervalD(Element min, Element max);

    Element _min;
    Element _max;
};

std::ostream& operator<<(std::ostream& os, IntervalI const& interval);

std::ostream& operator<<(std::ostream& os, IntervalD const& interval);

inline void swap(IntervalI& a, IntervalI& b) noexcept { a.swap(b); }

inline void swap(IntervalD& a, IntervalD& b) noexcept { a.swap(b); }

}  // namespace geom
}  // namespace lsst

namespace std {

template <>
inline void swap<lsst::geom::IntervalI>(lsst::geom::IntervalI& a, lsst::geom::IntervalI& b) noexcept {
    a.swap(b);
}

template <>
inline void swap<lsst::geom::IntervalD>(lsst::geom::IntervalD& a, lsst::geom::IntervalD& b) noexcept {
    a.swap(b);
}

template <>
struct hash<lsst::geom::IntervalI> {
    using argument_type = lsst::geom::IntervalI;
    using result_type = size_t;
    size_t operator()(argument_type const& x) const noexcept { return x.hash_value(); }
};

template <>
struct hash<lsst::geom::IntervalD> {
    using argument_type = lsst::geom::IntervalD;
    using result_type = size_t;
    size_t operator()(argument_type const& x) const noexcept { return x.hash_value(); }
};

}  // namespace std

#endif
