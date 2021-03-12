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

#ifndef LSST_GEOM_BOX_H
#define LSST_GEOM_BOX_H

#include <vector>

#include "boost/format.hpp"
#include "ndarray.h"

#include "lsst/geom/Point.h"
#include "lsst/geom/Extent.h"
#include "lsst/geom/Interval.h"

namespace lsst {
namespace geom {

class Box2D;

/**
 *  An integer coordinate rectangle.
 *
 *  Box2I is an inclusive box that represents a rectangular region of pixels.  A box
 *  never has negative dimensions; the empty box is defined to have zero-size dimensions,
 *  and is treated as though it does not have a well-defined position (regardless of the
 *  return value of getMin() or getMax() for an empty box).
 *
 *  @internal
 *
 *  Box2I internally stores its minimum point and dimensions, because we expect
 *  these will be the most commonly accessed quantities.
 *
 *  Box2I sets the minimum point to the origin for an empty box, and returns -1 for both
 *  elements of the maximum point in that case.
 */
class Box2I final {
public:
    typedef Point2I Point;
    typedef Extent2I Extent;
    typedef int Element;

    using Interval = IntervalI;

    enum EdgeHandlingEnum { EXPAND, SHRINK };

    /// Construct an empty box.
    Box2I() noexcept : _minimum(0), _dimensions(0) {}

    /**
     *  Construct a box from its minimum and maximum points.
     *
     *  @param[in] minimum   Minimum (lower left) coordinate (inclusive).
     *  @param[in] maximum   Maximum (upper right) coordinate (inclusive).
     *  @param[in] invert    If true (default), swap the minimum and maximum coordinates if
     *                       minimum > maximum instead of creating an empty box.
     */
    Box2I(Point2I const& minimum, Point2I const& maximum, bool invert = true);

    /**
     *  Construct a box from one corner and dimensions.
     *
     *  @param[in] corner    Reference coordinate. This is the lower left corner if both
     *                       dimensions are positive, but a right corner or upper corner if
     *                       the corresponding dimension is negative and `invert` is set.
     *  @param[in] dimensions Box dimensions.  If either dimension coordinate is 0, the box will be empty.
     *  @param[in] invert     If true (default), invert any negative dimensions instead of creating
     *                        an empty box.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the maximum Point2I would overflow.
     */
    Box2I(Point2I const& corner, Extent2I const& dimensions, bool invert = true);

    /**
     *  Construct a box from a pair of intervals.
     *
     *  @param[in] x   Extent in x direction.
     *  @param[in] y   Extent in y direction.
     */
    Box2I(Interval const & x, Interval const & y) :
        Box2I(Point(x.getMin(), y.getMin()), Point(x.getMax(), y.getMax()), false)
    {}

    /**
     *  Construct an integer box from a floating-point box.
     *
     *  Floating-point to integer box conversion is based on the concept that a pixel
     *  is not an infinitesimal point but rather a square of unit size centered on
     *  integer-valued coordinates.  Converting a floating-point box to an integer box
     *  thus requires a choice on how to handle pixels which are only partially contained
     *  by the input floating-point box.
     *
     *  @param[in] other          A floating-point box to convert.
     *  @param[in] edgeHandling   If EXPAND, the integer box will contain any pixels that
     *                            overlap the floating-point box.  If SHRINK, the integer
     *                            box will contain only pixels completely contained by
     *                            the floating-point box.
     *
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if `other` is not finite.
     */
    explicit Box2I(Box2D const& other, EdgeHandlingEnum edgeHandling = EXPAND);

    /// Standard copy constructor.
    Box2I(Box2I const&) noexcept = default;
    Box2I(Box2I&&) noexcept = default;
    ~Box2I() noexcept = default;

    /**
     * Create a box centered as closely as possible on a particular point.
     *
     * @param center The desired center of the box.
     * @param size The desired width and height (in that order) of the box.
     *
     * @returns if `size` is positive, a box with size `size`; otherwise,
     *          an empty box. If the returned box is not empty, its center
     *          shall be within half a pixel of `center` in either dimension.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting box would overflow.
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if `center` is not finite.
     */
    static Box2I makeCenteredBox(Point2D const& center, Extent const& size);

    void swap(Box2I& other) noexcept {
        _minimum.swap(other._minimum);
        _dimensions.swap(other._dimensions);
    }

    /// Standard assignment operator.
    Box2I& operator=(Box2I const&) noexcept = default;
    Box2I& operator=(Box2I&&) noexcept = default;

    /**
     *  @name Min/Max Accessors
     *
     *  Return the minimum and maximum coordinates of the box (inclusive).
     */
    //@{
    Point2I const getMin() const noexcept { return _minimum; }
    int getMinX() const noexcept { return _minimum.getX(); }
    int getMinY() const noexcept { return _minimum.getY(); }

    Point2I const getMax() const noexcept { return _minimum + _dimensions - Extent2I(1); }
    int getMaxX() const noexcept { return _minimum.getX() + _dimensions.getX() - 1; }
    int getMaxY() const noexcept { return _minimum.getY() + _dimensions.getY() - 1; }
    //@}

    /**
     *  @name Begin/End Accessors
     *
     *  Return STL-style begin (inclusive) and end (exclusive) coordinates for the box.
     */
    //@{
    Point2I const getBegin() const noexcept { return _minimum; }
    int getBeginX() const noexcept { return _minimum.getX(); }
    int getBeginY() const noexcept { return _minimum.getY(); }

    Point2I const getEnd() const noexcept { return _minimum + _dimensions; }
    int getEndX() const noexcept { return _minimum.getX() + _dimensions.getX(); }
    int getEndY() const noexcept { return _minimum.getY() + _dimensions.getY(); }
    //@}

    /**
     *  @name Size Accessors
     *
     *  Return the size of the box in pixels.
     */
    //@{
    Extent2I const getDimensions() const noexcept { return _dimensions; }
    int getWidth() const noexcept { return _dimensions.getX(); }
    int getHeight() const noexcept { return _dimensions.getY(); }
    int getArea() const { return getWidth() * getHeight(); }
    //@}

    /**
     *  @name Center Accessors
     *
     *  Return the center coordinate of the box.
     */
    //@{
    Point2D const getCenter() const noexcept;
    double getCenterX() const noexcept { return this->getCenter().getX(); }
    double getCenterY() const noexcept { return this->getCenter().getY(); }
    //@}

    /// 1-d interval accessors
    //@{
    Interval getX() const { return Interval::fromMinSize(getMinX(), getWidth()); }
    Interval getY() const { return Interval::fromMinSize(getMinY(), getHeight()); }
    //@}

    /// Return slices to extract the box's region from an ndarray::Array.
    ndarray::View<boost::fusion::vector2<ndarray::index::Range, ndarray::index::Range> > getSlices() const;

    /// Return true if the box contains no points.
    bool isEmpty() const noexcept { return _dimensions.getX() == 0 && _dimensions.getY() == 0; }

    //@{
    /// Return true if the box contains the point.
    bool contains(Point2I const& point) const noexcept;
    bool contains(Element x, Element y) const noexcept { return contains(Point2I(x, y)); }
    //@}

    /**
     *  Return true if all points contained by other are also contained by this.
     *
     *  An empty box is contained by every other box, including other empty boxes.
     */
    bool contains(Box2I const& other) const noexcept;

    //@{
    /**
     *  Return true if any points in other are also in this.
     *
     *  Any overlap operation involving an empty box returns false.
     */
    bool overlaps(Box2I const& other) const noexcept;
    bool intersects(Box2I const& other) const noexcept { return overlaps(other); }
    //@}

    /**
     *  Return true if there are no points in both `this` and `other`.
     */
    bool isDisjointFrom(Box2I const & other) const noexcept;

    /**
     *  Increase the size of the box by the given buffer amount in all directions.
     *
     *  If a negative buffer is passed and the final size of the box is less than or
     *  equal to zero, the box will be made empty.
     */
    void grow(int buffer) { grow(Extent2I(buffer)); }

    /**
     *  Increase the size of the box by the given buffer amount in each direction.
     *
     *  If a negative buffer is passed and the final size of the box is less than or
     *  equal to zero, the box will be made empty.
     */
    void grow(Extent2I const& buffer);

    /// Shift the position of the box by the given offset.
    void shift(Extent2I const& offset);

    /// Flip a bounding box about the y-axis given a parent box of extent (xExtent).
    void flipLR(int xExtent);

    /// Flip a bounding box about the x-axis given a parent box of extent (yExtent).
    void flipTB(int yExtent);

    /// Expand this to ensure that this->contains(point).
    void include(Point2I const& point);

    /// Expand this to ensure that this->contains(other).
    void include(Box2I const& other);

    /** Shrink this to ensure that `other.contains(*this)`.
     *
     * In particular, if `other` and this box do not overlap this box will become empty.
     *
     * @param other the box that must contain this one
     */
    void clip(Box2I const& other) noexcept;

    /**
     *  Increase the size of the box by the given amount(s) on all sides
     *  (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to eroding by `-buffer`.
     *
     *  If a final dimension of the box is less than or equal to zero, the new
     *  box will be empty.
     *
     *  Empty boxes remain empty after dilation.
     */
    //@{
    Box2I dilatedBy(Extent const & buffer) const;
    Box2I dilatedBy(Element buffer) const {
        return dilatedBy(Extent(buffer, buffer));
    }
    //@}

    /**
     *  Decrease the size of the box by the given amount(s) on all sides
     *  (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to dilating by `-buffer`.
     *
     *  If a final dimension of the box is less than or equal to zero, the new
     *  box will be empty.
     *
     *  Empty boxes remain empty after erosion.
     */
    //@{
    Box2I erodedBy(Extent const & buffer) const { return dilatedBy(-buffer); }
    Box2I erodedBy(Element buffer) const { return dilatedBy(-buffer); }
    //@}

    /**
     * Shift the position of the box by the given offset (returning a new
     * object).
     *
     * Empty boxes remain empty when shifted.
     */
    Box2I shiftedBy(Extent const & offset) const;

    /**
     * Reflect the box about a vertical line (returning a new object).
     *
     * Empty boxes remain empty when reflected.
     */
    Box2I reflectedAboutX(Element x) const;

    /**
     * Reflect the box about a horizontal line (returning a new object).
     *
     * Empty boxes remain empty when reflected.
     */
    Box2I reflectedAboutY(Element y) const;

    //@{
    /**
     * Expand the box to ensure that `contains(other)` is true (returning a
     * new object).
     *
     * Expanding an empty box with a single point yields an box with size=1 at
     * that point; expanding an empty box with a second box is equivalent to
     * assignment.
     */
    Box2I expandedTo(Point const & other) const;
    Box2I expandedTo(Box2I const & other) const;
    //@}

    /**
     * Shrink an interval to ensure that it is contained by other (returning a
     * new object).
     *
     * In particular, if `other` and this interval do not overlap this
     * interval will become empty.
     */
    Box2I clippedTo(Box2I const & other) const noexcept;

    /**
     *  Compare two boxes for equality.
     *
     *  All empty boxes are equal.
     */
    bool operator==(Box2I const& other) const noexcept;

    /**
     *  Compare two boxes for equality.
     *
     *  All empty boxes are equal.
     */
    bool operator!=(Box2I const& other) const noexcept;

    /// Return a hash of this object.
    std::size_t hash_value() const noexcept;

    /**
     * Get the corner points
     *
     * The order is counterclockise, starting from the lower left corner, i.e.:
     *   (minX, minY), (maxX, maxY), (maxX, maxX), (minX, maxY)
     */
    std::vector<Point2I> getCorners() const;

    std::string toString() const {
        return (boost::format("Box2I(%s,%s)") % _minimum.toString() % _dimensions.toString()).str();
    }

private:
    Point2I _minimum;
    Extent2I _dimensions;
};

/**
 *  A floating-point coordinate rectangle geometry.
 *
 *  Box2D is a half-open (minimum is inclusive, maximum is exclusive) box. A box
 *  never has negative dimensions; the empty box is defined to zero-size dimensions
 *  and its minimum and maximum values set to NaN.  Only the empty box may have
 *  zero-size dimensions.
 *
 *  @internal
 *
 *  Box2D internally stores its minimum point and maximum point, instead of
 *  minimum point and dimensions, to ensure roundoff error does not affect
 *  whether points are contained by the box.
 *
 *  Despite some recommendations to the contrary, Box2D sets the minimum and maximum
 *  points to NaN for an empty box.  In almost every case, special checks for
 *  emptiness would have been necessary anyhow, so there was little to gain in
 *  using the minimum > maximum condition to denote an empty box, as was used in Box2I.
 */
class Box2D final {
public:
    typedef Point2D Point;
    typedef Extent2D Extent;
    typedef double Element;

    using Interval = IntervalD;

    /**
     *  Value the maximum coordinate is multiplied by to increase it by the smallest
     *  possible amount.
     */
    static double const EPSILON;

    /// Value used to specify undefined coordinate values.
    static double const INVALID;

    /// Construct an empty box.
    Box2D() noexcept;

    /**
     *  Construct a box from its minimum and maximum points.
     *
     *  If any(minimum == maximum), the box will always be empty (even if invert==true).
     *
     *  @param[in] minimum   Minimum (lower left) coordinate (inclusive).
     *  @param[in] maximum   Maximum (upper right) coordinate (exclusive).
     *  @param[in] invert    If true (default), swap the minimum and maximum coordinates if
     *                       minimum > maximum instead of creating an empty box.
     */
    Box2D(Point2D const& minimum, Point2D const& maximum, bool invert = true) noexcept;

    /**
     *  Construct a box from one corner and dimensions.
     *
     *  @param[in] corner    Reference coordinate (inclusive). This is the lower left corner if
     *                       both dimensions are positive, but a right corner or upper corner if
     *                       the corresponding dimension is negative and `invert` is set.
     *  @param[in] dimensions Box dimensions.  If either dimension coordinate is 0, the box will be empty.
     *  @param[in] invert     If true (default), invert any negative dimensions instead of creating
     *                        an empty box.
     */
    Box2D(Point2D const& corner, Extent2D const& dimensions, bool invert = true) noexcept;

    /**
     *  Construct a box from a pair of intervals.
     *
     *  @param[in] x   Extent in x direction.
     *  @param[in] y   Extent in y direction.
     */
    Box2D(Interval const & x, Interval const & y) :
        Box2D(Point(x.getMin(), y.getMin()), Point(x.getMax(), y.getMax()), false)
    {}

    /**
     *  Construct a floating-point box from an integer box.
     *
     *  Integer to floating-point box conversion is based on the concept that a pixel
     *  is not an infinitesimal point but rather a square of unit size centered on
     *  integer-valued coordinates.  While the output floating-point box thus has
     *  the same dimensions as the input integer box, its minimum/maximum coordinates
     *  are 0.5 smaller/greater.
     */
    explicit Box2D(Box2I const& other) noexcept;

    /// Standard copy constructor.
    Box2D(Box2D const&) noexcept = default;
    Box2D(Box2D&&) noexcept = default;

    ~Box2D() noexcept = default;

    /**
     * Create a box centered on a particular point.
     *
     * @param center The desired center of the box.
     * @param size The desired width and height (in that order) of the box.
     *
     * @returns if `size` is positive, a box with size `size`; otherwise,
     *          an empty box. If the returned box is not empty, it shall be
     *          centered on `center`. Behavior is undefined if either `center`
     *          or `size` is non-finite.
     */
    // It's hard to guarantee postconditions (especially size) with non-finite inputs
    static Box2D makeCenteredBox(Point2D const& center, Extent const& size) noexcept;

    void swap(Box2D& other) noexcept {
        _minimum.swap(other._minimum);
        _maximum.swap(other._maximum);
    }

    /// Standard assignment operator.
    Box2D& operator=(Box2D const&) noexcept = default;
    Box2D& operator=(Box2D&&) noexcept = default;

    /**
     *  @name Min/Max Accessors
     *
     *  Return the minimum (inclusive) and maximum (exclusive) coordinates of the box.
     */
    //@{
    Point2D const getMin() const noexcept { return _minimum; }
    double getMinX() const noexcept { return _minimum.getX(); }
    double getMinY() const noexcept { return _minimum.getY(); }

    Point2D const getMax() const noexcept { return _maximum; }
    double getMaxX() const noexcept { return _maximum.getX(); }
    double getMaxY() const noexcept { return _maximum.getY(); }
    //@}

    /**
     *  @name Size Accessors
     *
     *  Return the size of the box.
     */
    //@{
    Extent2D const getDimensions() const noexcept { return isEmpty() ? Extent2D(0.0) : _maximum - _minimum; }
    double getWidth() const noexcept { return isEmpty() ? 0 : _maximum.getX() - _minimum.getX(); }
    double getHeight() const noexcept { return isEmpty() ? 0 : _maximum.getY() - _minimum.getY(); }
    double getArea() const noexcept {
        Extent2D dim(getDimensions());
        return dim.getX() * dim.getY();
    }
    //@}

    //@{
    /// 1-d interval accessors
    Interval getX() const { return Interval::fromMinMax(getMinX(), getMaxX()); }
    Interval getY() const { return Interval::fromMinMax(getMinY(), getMaxY()); }
    //@}

    /**
     *  @name Center Accessors
     *
     *  Return the center coordinate of the box.
     */
    //@{
    Point2D const getCenter() const noexcept {
        return Point2D((_minimum.asEigen() + _maximum.asEigen()) * 0.5);
    }
    double getCenterX() const noexcept { return (_minimum.getX() + _maximum.getX()) * 0.5; }
    double getCenterY() const noexcept { return (_minimum.getY() + _maximum.getY()) * 0.5; }
    //@}

    /// Return true if the box contains no points.
    bool isEmpty() const noexcept { return _minimum.getX() != _minimum.getX(); }

    //@{
    /// Return true if the box contains the point.
    bool contains(Point2D const& point) const noexcept;
    bool contains(Element x, Element y) const noexcept { return contains(Point2D(x, y)); }
    //@}

    /**
     *  Return true if all points contained by other are also contained by this.
     *
     *  An empty box is contained by every other box, including other empty boxes.
     */
    bool contains(Box2D const& other) const;

    //@{
    /**
     *  Return true if any points in other are also in this.
     *
     *  Any overlap operation involving an empty box returns false.
     */
    bool overlaps(Box2D const& other) const noexcept;
    bool intersects(Box2D const& other) const noexcept { return overlaps(other); }
    //@}


    /**
     *  Return true if there are no points in both `this` and `other`.
     */
    bool isDisjointFrom(Box2D const & other) const noexcept;

    /**
     *  Increase the size of the box by the given buffer amount in all directions.
     *
     *  If a negative buffer is passed and the final size of the box is less than or
     *  equal to zero, the box will be made empty.
     */
    void grow(double buffer) { grow(Extent2D(buffer)); }

    /**
     *  Increase the size of the box by the given buffer amount in each direction.
     *
     *  If a negative buffer is passed and the final size of the box is less than or
     *  equal to zero, the box will be made empty.
     */
    void grow(Extent2D const& buffer);

    /// Shift the position of the box by the given offset.
    void shift(Extent2D const& offset);

    /// Flip a bounding box about the y-axis given a parent box of extent (xExtent).
    void flipLR(float xExtent);

    /// Flip a bounding box about the x-axis given a parent box of extent (yExtent).
    void flipTB(float yExtent);

    /**
     *  Expand this to ensure that this->contains(point).
     *
     *  If the point sets a new maximum value for the box, the maximum coordinate will
     *  be adjusted to ensure the point is actually contained
     *  by the box instead of sitting on its exclusive upper edge.
     */
    void include(Point2D const& point) noexcept;

    /// Expand this to ensure that this->contains(other).
    void include(Box2D const& other) noexcept;

    /** Shrink this to ensure that `other.contains(*this)`.
     *
     * In particular, if `other` and this box do not overlap this box will become empty.
     *
     * @param other the box that must contain this one
     */
    void clip(Box2D const& other) noexcept;

    //@{
    /**
     *  Increase the size of the box by the given amount(s) on all sides
     *  (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to eroding by `-buffer`.
     *
     *  If the final size of the box is less than zero in either dimension the
     *  new box will empty.
     *
     *  Empty boxes remain empty after dilation.  Infinite bounds are
     *  unaffected by dilation.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `buffer` is not finite.
     */
    Box2D dilatedBy(Extent const & buffer) const;
    Box2D dilatedBy(Element buffer) const {
        return dilatedBy(Extent(buffer, buffer));
    }
    //@}

    //@{
    /**
     *  Decrease the size of the box by the given amount(s) on all sides
     *  (returning a new object).
     *
     *  If `buffer` is negative, this is equivalent to dilating by `-buffer`.
     *
     *  If the final size of the box is less than zero in either dimension the
     *  new box will empty.
     *
     *  Empty boxes remain empty after erosion.  Infinite bounds are
     *  unaffected by erosion.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `buffer` is not finite.
     */
    Box2D erodedBy(Extent const & buffer) const { return dilatedBy(-buffer); }
    Box2D erodedBy(Element buffer) const { return dilatedBy(-buffer); }
    //@}

    /**
     *  Shift the position of the box by the given offset (returning a new
     *  object).
     *
     *  Empty boxes remain empty when shifted.  Infinite bounds are
     *  unaffected by shifting.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `offset` is not finite.
     */
    Box2D shiftedBy(Extent const & offset) const;

    /**
     * Reflect the box about a vertical line (returning a new object).
     *
     * Empty boxes remain empty when reflected.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `x` is not finite.
     */
    Box2D reflectedAboutX(Element x) const;

    /**
     * Reflect the box about a horizontal line (returning a new object).
     *
     * Empty boxes remain empty when reflected.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if
     *      `y` is not finite.
     */
    Box2D reflectedAboutY(Element y) const;

    /**
     *  Expand a box to ensure that `contains(other)` is true (returning a new
     *  object).
     *
     *  Expanding an empty box with a single point yields a box with
     *  `dimensions == (0, 0)` at that point.
     *
     *  @throws lsst::pex::exceptions::InvalidParameterError  Thrown if `other`
     *       is not finite.
     */
    Box2D expandedTo(Point const & other) const;

    /**
     *  Expand a box to ensure that `contains(other)` is true (returning a new
     *  object).
     *
     *  Expanding an empty box with a second box is equivalent to assignment.
     */
    Box2D expandedTo(Box2D const & other) const;

    /**
     * Shrink a box to ensure that it is contained by other (returning a new
     * object).
     *
     * In particular, if `other` and `this` do not overlap, the new box will
     * be empty.
     */
    Box2D clippedTo(Box2D const & other) const;

    /**
     *  Compare two boxes for equality.
     *
     *  All empty boxes are equal.
     */
    bool operator==(Box2D const& other) const noexcept;

    /**
     *  Compare two boxes for equality.
     *
     *  All empty boxes are equal.
     */
    bool operator!=(Box2D const& other) const noexcept;

    /// Return a hash of this object.
    std::size_t hash_value() const noexcept;

    /**
     * Get the corner points
     *
     * The order is counterclockise, starting from the lower left corner, i.e.:
     *   (minX, minY), (maxX, maxY), (maxX, maxX), (minX, maxY)
     */
    std::vector<Point2D> getCorners() const;

    std::string toString() const {
        return (boost::format("Box2D(%s,%s)") % _minimum.toString() % _maximum.toString()).str();
    }

private:
    void _tweakMax(int n) noexcept {
        if (_maximum[n] < 0.0) {
            _maximum[n] *= (1.0 - EPSILON);
        } else if (_maximum[n] > 0.0) {
            _maximum[n] *= (1.0 + EPSILON);
        } else {
            _maximum[n] = EPSILON;
        }
    }
    Point2D _minimum;
    Point2D _maximum;
};

typedef Box2D BoxD;
typedef Box2I BoxI;

std::ostream& operator<<(std::ostream& os, Box2I const& box);

std::ostream& operator<<(std::ostream& os, Box2D const& box);

}  // namespace geom
}  // namespace lsst

namespace std {
template <>
struct hash<lsst::geom::Box2I> {
    using argument_type = lsst::geom::Box2I;
    using result_type = size_t;
    size_t operator()(argument_type const& x) const noexcept { return x.hash_value(); }
};

template <>
struct hash<lsst::geom::Box2D> {
    using argument_type = lsst::geom::Box2D;
    using result_type = size_t;
    size_t operator()(argument_type const& x) const noexcept { return x.hash_value(); }
};
}  // namespace std

#endif
