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
#ifndef LSST_AFW_MATH_POLYNOMIALS_SafeSum_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_SafeSum_h_INCLUDED

#include <cmath>

namespace lsst { namespace geom { namespace polynomials {

/**
 *  Enum used to control how to sum polynomial terms.
 */
enum class SumMode {
    /// Summation using regular floating-point addition.
    FAST,
    /// Compensated summation using SafeSum.  Involves ~4x as many floating point operations.
    SAFE
};

/**
 *  A numerically stable summation algorithm for floating-point numbers.
 *
 *  SafeSum implements += and -= operators that can be used to accumulate
 *  floating point numbers with very different magnitudes, with accuracy
 *  limited only by the usual floating-point inaccuracy in represented
 *  the final sum.
 *
 *  SafeSum is *explicitly* convertible to and from its underlying
 *  floating-point type and only supports in-place addition and subtraction,
 *  in order to avoid cases where a combination of implicit conversion and
 *  multiple overloaded operators could lead to accidental use of regular
 *  floating-point operations.
 *
 *  SafeSum uses the Kahan-Neumaier algorithm (though this should be
 *  considered an implementation detail by callers), which maintains a
 *  lower-order-bit correction that compensates for the loss of precision in
 *  the main sum.  Particularly aggressive compiler optimizations that do not
 *  preserve IEEE floating point behavior (e.g. gcc's `-fassociative-math`)
 *  may optimize away the correction and reduce SafeSum's behavior to a
 *  standard unsafe sum.
 */
template <typename T>
class SafeSum {
public:

    explicit SafeSum(T initial=static_cast<T>(0)) noexcept :
        _sum(initial),
        _correction(static_cast<T>(0))
    {}

    SafeSum & operator=(T value) noexcept {
        _sum = value;
        _correction = static_cast<T>(0);
        return *this;
    }

    SafeSum & operator+=(T value) noexcept {
        T t = _sum + value;
        // update _correction to account for lost low-order bits of the
        // greater of _sum and value
        if (std::abs(_sum) >= std::abs(value)) {
            _correction += (_sum - t) + value;
        } else {
            _correction += (value - t) + _sum;
        }
        _sum = t;
        return *this;
    }

    SafeSum & operator-=(T value) noexcept {
        return operator+=(-value);
    }

    explicit operator T() const noexcept {
        return _sum + _correction;
    }

private:
    T _sum;
    T _correction;
};

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_SafeSum_h_INCLUDED
