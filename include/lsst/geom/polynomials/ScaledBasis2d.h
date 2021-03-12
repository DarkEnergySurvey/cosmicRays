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
#ifndef LSST_AFW_MATH_POLYNOMIALS_ScaledBasis2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_ScaledBasis2d_h_INCLUDED

#include "lsst/geom/polynomials/Scaling2d.h"

namespace lsst { namespace geom { namespace polynomials {

template <typename Basis>
class Function2d;

/**
 *  A 2-d basis that transforms all input points before evaluating nested basis.
 *
 *  If the nested basis is defined by basis functions @f$B_m(x)B_n(y)@f$, the
 *  scaled basis functions are @f$B_m(U(x))B_n(V(y))@f$, where @f$U(x)@f$ and
 *  @f$V(y)@f$ together represent the scaling transform.
 *
 *  Both the nested basis and ScaledBasis2d itself are models of the Basis2d
 *  concept.
 */
template <typename Nested>
class ScaledBasis2d {
public:

    /// A Function2d object that uses this basis.
    using Function = Function2d<ScaledBasis2d>;

    /// The type returned by scale().
    using Scaled = ScaledBasis2d<Nested>;

    /// The type returned by makeWorkspace().
    using Workspace = typename Nested::Workspace;

    /// The type returned by getIndices().
    using IndexRange = typename Nested::IndexRange;

    /// Construct a scaled basis from a nested basis and a scaling transform.
    explicit ScaledBasis2d(Nested const & nested, Scaling2d const & scaling) :
        _nested(nested),
        _scaling(scaling)
    {}

    /**
     *  Construct a basis that remaps the given box to [-1, 1]x[-1, 1] before
     *  evaluating the nested polynomials polynomials.
     *
     *  @param[in]  order    Maximum order of the basis (inclusive).
     *  @param[in]  box      Box to be mapped to [-1, 1]x[-1, 1].
     *
     *  This constructor is particularly useful for Chebyshev polynomials, for
     *  which most of the special functions of the basis are only active when
     *  the domain is limited to [-1, 1].
     *
     *  This signature requires that Nested(order) be a valid constructor.
     */
    ScaledBasis2d(std::size_t order, Box2D const & box) :
        _nested(order),
        _scaling(makeUnitRangeScaling2d(box))
    {}

    /// Default copy constructor.
    ScaledBasis2d(ScaledBasis2d const &) = default;

    /// Default move constructor.
    ScaledBasis2d(ScaledBasis2d &&) = default;

    /// Default copy assignment.
    ScaledBasis2d & operator=(ScaledBasis2d const &) = default;

    /// Default move assignment.
    ScaledBasis2d & operator=(ScaledBasis2d &&) = default;

    /// Return the nested basis.
    Nested const & getNested() const noexcept { return _nested; }

    /// Return the scaling transform.
    Scaling2d const & getScaling() const noexcept { return _scaling; }

    /// Return the order of the basis.
    std::size_t getOrder() const { return getNested().getOrder(); }

    /// Return the number of elements in the basis.
    std::size_t size() const { return getNested().size(); }

    /**
     *  Return a scaled basis that delegates to a copy of `this`.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scaled(Scaling2d const & first) const {
        return getNested().scaled(first.then(getScaling()));
    }

    /// Return the flattened index of the basis function with the given x and y orders.
    int index(int x, int y) const { return getNested().index(x, y); }

    /**
     *  Return a range of iterators that dereference to Index2d.
     *
     *  This is the recommended way to interpret the packed coefficients and basis functions
     *  utilized by PackedBasis2d; for example,
     *  ```
     *      // Evaluate basis functions at a point.
     *      geom::Point2D point(x, y);
     *      PolynomialBasis2d basis(order);
     *      std::vector<double> values(basis.size());
     *      basis.fill(point, values);
     *      // Iterate over tuples of flattened indices and x and y orders.
     *      for (auto const & index : basis.getIndices()) {
     *          double a = values[index.flat];
     *          // standard polynomial basis functions are just powers
     *          double b = std::pow(point.getX(), index.nx)*std::pow(point.getY(), index.ny)];
     *          assert(std::abs(a - b) < std::numeric_limits<double>::epsilon());
     *      }
     *  ```
     */
    IndexRange getIndices() const { return getNested().getIndices(); }

    /// Allocate a workspace that can be passed to sumWith() and fill() to avoid repeated memory allocations.
    Workspace makeWorkspace() const { return getNested().makeWorkspace();}

    /**
     *  Evaluate a basis expansion with the given coefficients.
     *
     *  If the 1-d basis elements are @f$B_n(x)@f$ and the given coefficients are
     *  a vector @f$a_{p, q}@f$, this computes
     *  @f[
     *      \sum_{p = 0, q = 0}^{p + q \le N} a_{p,q} B_{p}(x) B_{q}(y)
     *  @f]
     *
     *  @param[in] point         Point at which to evaluate the expansion.
     *  @param[in] coefficients  Flattened coefficients vector.
     *                           See Basis1d::sumWith for more information.
     *  @param[in] mode          Enum indicating the tradeoff to make between
     *                           speed and numerical precision.
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              the same exception safety as it if it does.
     */
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients,
                   SumMode mode=SumMode::FAST) const {
        return getNested().sumWith(getScaling().applyForward(point), coefficients, mode);
    }

    /// Evaluate a basis expansion with the given coefficients (external workspace version).
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients,
                   Workspace & workspace, SumMode mode=SumMode::FAST) const {
        return getNested().sumWith(getScaling().applyForward(point), coefficients, workspace, mode);
    }

    /**
     *  Evaluate the basis at a given point.
     *
     *  @param[in] point  Point at which to evaluate the basis functions.
     *  @param[out] basis Flattened output vector.
     *                    See Basis1d::fill more information.
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              basic exception safety if it does.
     */
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis) const {
        return getNested().fill(getScaling().applyForward(point),
                                std::forward<Vector>(basis));
    }

    /// Evaluate the basis at a given point (external workspace version).
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis, Workspace & workspace) const {
        return getNested().fill(getScaling().applyForward(point),
                                std::forward<Vector>(basis),
                                workspace);
    }

private:
    Nested _nested;
    Scaling2d _scaling;
};

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_ScaledBasis2d_h_INCLUDED
