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

#include <cmath>
#include <limits>
#include <type_traits>

#include "boost/format.hpp"

#include "lsst/utils/hashCombine.h"
#include "lsst/geom/Interval.h"

namespace lsst {
namespace geom {

namespace {

using BigElement = long long;

// Templated so we can check both floating point and BigElement values.
// Should never be invoked with IntervalI::Element values, as that won't
// accomplish anything.
template <typename T>
void checkForOverflow(T x, char const* where) {
    static_assert(sizeof(T) > sizeof(IntervalI::Element) || !std::is_integral<T>::value);
    if (x < std::numeric_limits<IntervalI::Element>::min() ||
        x > std::numeric_limits<IntervalI::Element>::max()) {
        throw LSST_EXCEPT(pex::exceptions::OverflowError,
                          (boost::format("Integer overflow (%d) in interval %s.") % x % where).str());
    }
}

}  // namespace

IntervalI IntervalI::fromMinMax(Element min, Element max) {
    return _fromMinMaxChecked(static_cast<BigElement>(min), static_cast<BigElement>(max));
}

IntervalI IntervalI::fromMinSize(Element min, Element size) {
    if (size <= 0) {
        return IntervalI();
    }
    BigElement max = static_cast<BigElement>(min) + static_cast<BigElement>(size) - 1;
    checkForOverflow(max, "maximum");
    return IntervalI(min, size);
}

IntervalI IntervalI::fromMaxSize(Element max, Element size) {
    if (size <= 0) {
        return IntervalI();
    }
    BigElement min = static_cast<BigElement>(max) - static_cast<BigElement>(size) + 1;
    checkForOverflow(min, "minimum");
    return IntervalI(static_cast<Element>(min), size);
}

IntervalI IntervalI::fromCenterSize(double center, Element size) {
    if (size <= 0) {
        return IntervalI();
    }
    double min = center - 0.5 * size;
    // compensate for IntervalI's coordinate conventions (where max = min + size - 1)
    min += 0.5;
    checkForOverflow(min, "minimum");
    checkForOverflow(min + size - 1, "maximum");
    return IntervalI(static_cast<BigElement>(min), size);
}

IntervalI::IntervalI(IntervalD const& other, EdgeHandlingEnum edgeHandling) : _min(), _size() {
    if (other.isEmpty()) {
        *this = IntervalI();
        return;
    }
    if (!std::isfinite(other.getMin()) || !std::isfinite(other.getMax())) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Cannot convert non-finite IntervalD to IntervalI");
    }
    BigElement min, max;
    switch (edgeHandling) {
        case EdgeHandlingEnum::EXPAND:
            min = static_cast<BigElement>(std::ceil(other.getMin() - 0.5));
            max = static_cast<BigElement>(std::floor(other.getMax() + 0.5));
            break;
        case EdgeHandlingEnum::SHRINK:
            min = static_cast<BigElement>(std::ceil(other.getMin() + 0.5));
            max = static_cast<BigElement>(std::floor(other.getMax() - 0.5));
            break;
        default:
            throw pex::exceptions::LogicError("Invalid enum value.");
    }
    *this = _fromMinMaxChecked(min, max);
}

ndarray::View<boost::fusion::vector1<ndarray::index::Range>> IntervalI::getSlice() const {
    return ndarray::view(getBegin(), getEnd());
}

bool IntervalI::contains(Element point) const noexcept {
    // The case where this->isEmpty() is handled implicitly by the invariant
    // that empty intervals have max < min.
    return point >= this->getMin() && point <= this->getMax();
}

bool IntervalI::contains(IntervalI const& other) const noexcept {
    // The case where this->isEmpty() is handled implicitly by the invariant
    // that empty intervals have max < min.
    return other.isEmpty() || (other.getMin() >= this->getMin() && other.getMax() <= this->getMax());
}

bool IntervalI::overlaps(IntervalI const& other) const noexcept { return !isDisjointFrom(other); }

bool IntervalI::isDisjointFrom(IntervalI const& other) const noexcept {
    if (isEmpty() || other.isEmpty()) {
        return true;
    }
    return getMin() > other.getMax() || getMax() < other.getMin();
}

IntervalI IntervalI::dilatedBy(Element buffer) const {
    if (isEmpty()) {
        return IntervalI();
    }
    BigElement min = static_cast<BigElement>(getMin()) - buffer;
    BigElement max = static_cast<BigElement>(getMax()) + buffer;
    return _fromMinMaxChecked(min, max);
}

IntervalI IntervalI::shiftedBy(Element offset) const {
    if (isEmpty()) {
        return IntervalI();
    }
    BigElement min = static_cast<BigElement>(getMin()) + offset;
    BigElement max = static_cast<BigElement>(getMax()) + offset;
    checkForOverflow(min, "minimum");
    checkForOverflow(max, "maximum");
    return IntervalI(static_cast<Element>(min), getSize());
}

IntervalI IntervalI::reflectedAbout(Element point) const {
    if (isEmpty()) {
        return IntervalI();
    }
    BigElement max = 2 * static_cast<BigElement>(point) - getMin();
    BigElement min = 2 * static_cast<BigElement>(point) - getMax();
    return _fromMinMaxChecked(min, max);
}

IntervalI IntervalI::expandedTo(Element point) const {
    if (isEmpty()) {
        return IntervalI(point, 1);
    }
    return _fromMinMaxChecked(std::min(static_cast<BigElement>(point), static_cast<BigElement>(getMin())),
                              std::max(static_cast<BigElement>(point), static_cast<BigElement>(getMax())));
}

IntervalI IntervalI::expandedTo(IntervalI const& other) const {
    if (other.isEmpty()) {
        return *this;
    }
    if (this->isEmpty()) {
        return other;
    }
    return _fromMinMaxChecked(
            std::min(static_cast<BigElement>(other.getMin()), static_cast<BigElement>(getMin())),
            std::max(static_cast<BigElement>(other.getMax()), static_cast<BigElement>(getMax())));
}

IntervalI IntervalI::clippedTo(IntervalI const& other) const noexcept {
    if (isEmpty() || other.isEmpty()) {
        return IntervalI();
    }
    return fromMinMax(std::max(getMin(), other.getMin()), std::min(getMax(), other.getMax()));
}

bool IntervalI::operator==(IntervalI const& other) const noexcept {
    return other._min == this->_min && other._size == this->_size;
}

bool IntervalI::operator!=(IntervalI const& other) const noexcept { return !(other == *this); }

std::size_t IntervalI::hash_value() const noexcept {
    // Completely arbitrary seed
    return utils::hashCombine(17, _min, _size);
}

std::string IntervalI::toString() const {
    return (boost::format("(min=%s, max=%s)") % getMin() % getMax()).str();
}

template <typename T>
IntervalI IntervalI::_fromMinMaxChecked(T min, T max) {
    if (max < min) {
        return IntervalI();
    }
    checkForOverflow(min, "minimum");
    checkForOverflow(max, "maximum");
    T size = 1 + max - min;
    checkForOverflow(size, "size");
    return IntervalI(static_cast<Element>(min), static_cast<Element>(size));
}

IntervalI::IntervalI(Element min, Element size) : _min(min), _size(size) {}

IntervalD::IntervalD() noexcept
        : _min(std::numeric_limits<IntervalD::Element>::quiet_NaN()),
          _max(std::numeric_limits<IntervalD::Element>::quiet_NaN()) {}

IntervalD IntervalD::fromMinMax(Element min, Element max) { return IntervalD(min, max); }

IntervalD IntervalD::fromMinSize(Element min, Element size) {
    if (std::isinf(min) || (std::isinf(size) && size > 0)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Ambiguously infinite interval parameters; use fromMinMax to "
                          "construct infinite intervals instead.");
    }
    return IntervalD(min, min + size);
}

IntervalD IntervalD::fromMaxSize(Element max, Element size) {
    if (std::isinf(max) || (std::isinf(size) && size > 0)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Ambiguously infinite interval parameters; use fromMinMax to "
                          "construct infinite intervals instead.");
    }
    return IntervalD(max - size, max);
}

IntervalD IntervalD::fromCenterSize(Element center, Element size) {
    if (std::isinf(center) || std::isinf(size)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Infinite values not supported.");
    }
    Element min = center - 0.5 * size;
    return lsst::geom::IntervalD::fromMinSize(min, size);
}

IntervalD::IntervalD(IntervalI const& other) noexcept
        : _min(other.getMin() - 0.5), _max(other.getMax() + 0.5) {
    if (other.isEmpty()) *this = IntervalD();
}

IntervalD::Element IntervalD::getSize() const noexcept { return isEmpty() ? 0.0 : (_max - _min); }

IntervalD::Element IntervalD::getCenter() const noexcept { return 0.5 * (_min + _max); }

bool IntervalD::contains(Element point) const {
    if (std::isnan(point)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Cannot test whether an interval contains NaN.");
    }
    return point >= this->getMin() && point <= this->getMax();
}

bool IntervalD::contains(IntervalD const& other) const noexcept {
    return other.isEmpty() || (other.getMin() >= this->getMin() && other.getMax() <= this->getMax());
}

bool IntervalD::overlaps(IntervalD const& other) const noexcept { return !isDisjointFrom(other); }

bool IntervalD::isDisjointFrom(IntervalD const& other) const noexcept {
    if (isEmpty() || other.isEmpty()) {
        return true;
    }
    return getMin() > other.getMax() || getMax() < other.getMin();
}

IntervalD IntervalD::dilatedBy(Element buffer) const {
    if (!std::isfinite(buffer)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Cannot dilate or erode with a non-finite buffer.");
    }
    return fromMinMax(_min - buffer, _max + buffer);
}

IntervalD IntervalD::shiftedBy(Element offset) const {
    if (!std::isfinite(offset)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Cannot shift with a non-finite offset.");
    }
    if (!isEmpty()) {
        return fromMinMax(_min + offset, _max + offset);
    } else {
        return IntervalD();
    }
}

IntervalD IntervalD::reflectedAbout(Element point) const {
    if (!std::isfinite(point)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Cannot reflect about a non-finite point.");
    }
    if (!isEmpty()) {
        return fromMinMax(2 * point - _max, 2 * point - _min);
    } else {
        return IntervalD();
    }
}

IntervalD IntervalD::expandedTo(Element point) const {
    if (!std::isfinite(point)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Cannot expand to a non-finite point.");
    }
    if (isEmpty()) {
        return fromMinMax(point, point);
    } else {
        return fromMinMax(std::min(point, _min), std::max(point, _max));
    }
}

IntervalD IntervalD::expandedTo(IntervalD const& other) const noexcept {
    if (other.isEmpty()) {
        return *this;
    } else if (this->isEmpty()) {
        return other;
    } else {
        return fromMinMax(std::min(this->getMin(), other.getMin()), std::max(this->getMax(), other.getMax()));
    }
}

IntervalD IntervalD::clippedTo(IntervalD const& other) const noexcept {
    if (this->isEmpty() || other.isEmpty()) {
        return IntervalD();
    } else {
        return fromMinMax(std::max(this->getMin(), other.getMin()), std::min(this->getMax(), other.getMax()));
    }
}

bool IntervalD::operator==(IntervalD const& other) const noexcept {
    return (other.isEmpty() && this->isEmpty()) || (other._min == this->_min && other._max == this->_max);
}

bool IntervalD::operator!=(IntervalD const& other) const noexcept { return !(other == *this); }

std::size_t IntervalD::hash_value() const noexcept {
    // Completely arbitrary seed
    return utils::hashCombine(17, _min, _max);
}

std::string IntervalD::toString() const {
    return (boost::format("(min=%s, max=%s)") % getMin() % getMax()).str();
}

IntervalD::IntervalD(Element min, Element max) : _min(min), _max(max) {
    if (max < min || std::isnan(min) || std::isnan(max)) {
        *this = IntervalD();
    } else if (std::isinf(min) && min > 0.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Cannot set interval minimum to +infinity.");
    } else if (std::isinf(max) && max < 0.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Cannot set interval maximum to -infinity.");
    }
}

std::ostream& operator<<(std::ostream& os, IntervalI const& interval) {
    if (interval.isEmpty()) return os << "IntervalI()";
    return os << "IntervalI" << interval.toString();
}

std::ostream& operator<<(std::ostream& os, IntervalD const& interval) {
    if (interval.isEmpty()) return os << "IntervalD()";
    return os << "IntervalD" << interval.toString();
}

}  // namespace geom
}  // namespace lsst
