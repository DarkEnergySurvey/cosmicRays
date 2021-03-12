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
#ifndef LSST_AFW_MATH_POLYNOMIALS_ScaledBasis1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_ScaledBasis1d_h_INCLUDED

#include "lsst/geom/polynomials/SafeSum.h"
#include "lsst/geom/polynomials/Scaling1d.h"

namespace lsst { namespace geom { namespace polynomials {

template <typename Basis>
class Function1d;

/**
 *  A 1-d basis that transforms all input points before evaluating nested basis.
 *
 *  If the nested basis is defined by basis functions @f$B_n(x)@f$, the scaled
 *  basis functions are @f$B_n(S(x))@f$, where @f$S(x)@f$ is the scaling
 *  transform.
 *
 *  Both the nested basis and ScaledBasis1d itself are models of the Basis1d
 *  concept.
 */
template <typename Nested>
class ScaledBasis1d {
public:

    /// A Function1d object that uses this basis.
    using Function = Function1d<ScaledBasis1d>;

    /// The type returned by scale().
    using Scaled = ScaledBasis1d<Nested>;

    /// Construct a scaled basis from a nested basis and a scaling transform.
    explicit ScaledBasis1d(Nested const & nested, Scaling1d const & scaling) :
        _nested(nested),
        _scaling(scaling)
    {}

    /**
     *  Construct a basis that remaps the given interval to [-1, 1] before
     *  evaluating the nested basis.
     *
     *  @param[in]  order    Maximum order of the basis (inclusive).
     *  @param[in]  min      Minimum point of the interval, mapped to -1.
     *  @param[in]  max      Maximum point of the interval, mapped to 1.
     *
     *  This constructor is particularly useful for Chebyshev polynomials, for
     *  which most of the special functions of the basis are only active when
     *  the domain is limited to [-1, 1].
     *
     *  This signature requires that Nested(order) be a valid constructor.
     */
    ScaledBasis1d(std::size_t order, double min, double max) :
        _nested(order),
        _scaling(makeUnitRangeScaling1d(min, max))
    {}

    /// Default copy constructor.
    ScaledBasis1d(ScaledBasis1d const &) = default;

    /// Default move constructor.
    ScaledBasis1d(ScaledBasis1d &&) = default;

    /// Default copy assignment.
    ScaledBasis1d & operator=(ScaledBasis1d const &) = default;

    /// Default move assignment.
    ScaledBasis1d & operator=(ScaledBasis1d &&) = default;

    /// Return the nested basis.
    Nested const & getNested() const noexcept { return _nested; }

    /// Return the scaling transform.
    Scaling1d const & getScaling() const noexcept { return _scaling; }

    /// Return the order of the basis.
    std::size_t getOrder() const { return getNested().getOrder(); }

    /// Return the number of elements in the basis.
    std::size_t size() const { return getNested().size(); }

    /**
     *  Return a further-scaled basis with the same order.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scaled(Scaling1d const & first) const {
        return getNested().scaled(first.then(getScaling()));
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
        return getNested().sumWith(getScaling().applyForward(x), coefficients, mode);
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
        return getNested().fill(getScaling().applyForward(x), std::forward<Vector>(basis));
    }

private:
    Nested _nested;
    Scaling1d _scaling;
};

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_ScaledBasis1d_h_INCLUDED
