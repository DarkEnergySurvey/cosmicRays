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
#ifndef LSST_AFW_MATH_POLYNOMIALS_RecurrenceBasis1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_RecurrenceBasis1d_h_INCLUDED

#include "lsst/geom/polynomials/SafeSum.h"
#include "lsst/geom/polynomials/ScaledBasis1d.h"

namespace lsst { namespace geom { namespace polynomials {

template <typename Basis>
class Function1d;

#ifdef DOXYGEN

/**
 *  A recurrence relation concept for RecurrenceBasis1d.
 *
 *  @note This class is only present in the documentation, as it represents an
 *        abstract interface for which C++ (prior to C++20, at least) has no
 *        language support.  It may be formalized into a true Concept when
 *        that language feature is available.
 */
struct Recurrence {

    /// Return the zeroth element of the basis, @f$B_0(x)@f$
    static double getB0(double x);

    /// Return the first element of the basis, @f$B_1(x)@f$
    static double getB1(double x);

    /**
     *  Return the next element in the recurrence.
     *
     *  @param[in]  x         The point at which the basis is being evaluated.
     *  @param[in]  n         The order of the current basis function.
     *  @param[in]  current   @f$B_{n}(x)@f$, the current basis function value at x.
     *  @param[in]  previous  @f$B_{n-1}(x)@f$, the previous basis function value at x.
     *
     *  @return @f$B_{n+1}(x)@f$, the next value of the basis function at x.
     */
    static double next(double x, std::size_t n, double current, double previous);

};

#endif  // DOXYGEN


/**
 *  A basis for 1-d series expansions defined by a recurrence relation.
 *
 *  The recurrence relations utilized by RecurrenceBasis1d must have the
 *  following form:
 *  @f[
 *      B_{n+1}(x) = R(x, n, B_n(x), B_{n-1}(x))
 *  @f]
 *  with explicit expressions for @f$B_0(x)@f$ and @f$B_1(x)@f$ also given.
 *  This includes all special polynomials (e.g. Chebyshev, Legendre, Hermite,
 *  Laguerre) and products of special polynomials with their natural weight
 *  functions (e.g. Gauss-Hermite functions).  The template parameter must be
 *  a model of the Recurrence concept.
 *
 *  RecurrenceBasis1d is a model of the Basis1d concept.
 */
template <typename Recurrence>
class RecurrenceBasis1d {
public:

    /// A Function1d object that uses this basis.
    using Function = Function1d<RecurrenceBasis1d>;

    /// The type returned by scale().
    using Scaled = ScaledBasis1d<RecurrenceBasis1d>;

    /// Construct a basis with the given order (inclusive).
    explicit RecurrenceBasis1d(std::size_t order) noexcept :
        _order(order)
    {}

    /// Default copy constructor.
    RecurrenceBasis1d(RecurrenceBasis1d const &) = default;

    /// Default move constructor.
    RecurrenceBasis1d(RecurrenceBasis1d &&) = default;

    /// Default copy assignment.
    RecurrenceBasis1d & operator=(RecurrenceBasis1d const &) = default;

    /// Default move assignment.
    RecurrenceBasis1d & operator=(RecurrenceBasis1d &&) = default;

    /// Return the order of the basis.
    std::size_t getOrder() const noexcept { return _order; }

    /// Return the number of elements in the basis.
    std::size_t size() const noexcept { return _order + 1; }

    /**
     *  Return a scaled basis with the same order and recurrence.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scaled(Scaling1d const & scaling) const noexcept {
        return Scaled(*this, scaling);
    }

    /**
     *  Evaluate a basis expansion with the given coefficients.
     *
     *  If the basis elements are @f$B_n(x)@f$ and the given coefficients are
     *  a vector @f$a_n@f$, this computes
     *  @f[
     *      \sum_{n = 0}^{n \le N} a_n B_n(x)
     *  @f]
     *
     *  @param[in] x             Point at which to evaluate the expansion.
     *  @param[in] coefficients  Coefficients vector.
     *                           See Basis1d::sumWith for more information.
     *  @param[in] mode          Enum indicating the tradeoff to make between
     *                           speed and numerical precision.
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              the same exception safety as it if it does.
     */
    template <typename Vector>
    double sumWith(double x, Vector const & coefficients, SumMode mode=SumMode::FAST) const {
        // This universal lambda lets us effectively template most of the
        // implementation of this function on double vs. SafeSum<double>
        // without having to define an external template.
        auto accumulate = [x, coefficients, this](auto & sum) {
            double previous = Recurrence::getB0(x);
            if (_order > 0u) {
                double current = Recurrence::getB1(x);
                sum += coefficients[1]*current;
                for (std::size_t n = 2; n <= _order; ++n) {
                    double next = Recurrence::next(x, n, current, previous);
                    sum += coefficients[n]*next;
                    previous = current;
                    current = next;
                }
            }
        };
        double result = 0.0;
        if (mode == SumMode::FAST) {
            double z = Recurrence::getB0(x)*coefficients[0];
            accumulate(z);
            result = z;
        } else {
            SafeSum<double> z(Recurrence::getB0(x)*coefficients[0]);
            accumulate(z);
            result = static_cast<double>(z);
        }
        return result;
    }

    /**
     *  Evaluate the basis at a given point.
     *
     *  @param[in] x      Point at which to evaluate the basis functions.
     *  @param[out] basis Output vector.
     *                    See Basis1d::fill more information.
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              basic exception safety if it does.
     */
    template <typename Vector>
    void fill(double x, Vector && basis) const {
        std::forward<Vector>(basis)[0] = Recurrence::getB0(x);
        if (_order > 0u) {
            std::forward<Vector>(basis)[1] = Recurrence::getB1(x);
            for (std::size_t n = 2; n <= _order; ++n) {
                std::forward<Vector>(basis)[n] = Recurrence::next(
                    x, n,
                    std::forward<Vector>(basis)[n - 1],
                    std::forward<Vector>(basis)[n - 2]
                );
            }
        }
    }

private:
    std::size_t _order;
};

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_RecurrenceBasis1d_h_INCLUDED
