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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Basis2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Basis2d_h_INCLUDED
#ifdef DOXYGEN

namespace lsst { namespace geom { namespace polynomials {

/**
 *  A basis interface for 2-d series expansions.
 *
 *  @note This class is only present in the documentation, as it represents an
 *        abstract interface for which C++ (prior to C++20, at least) has no
 *        language support.  It may be formalized into a true Concept when
 *        that language feature is available.
 */
template <typename Basis1d>
class Basis2d {
public:

    /// A Function2d object that uses this basis.
    using Function = ...;

    /// The type returned by scale().
    using Scaled = ...;

    /// The type returned by makeWorkspace().
    using Workspace = ...;

    /// Return the maximum order of the basis.
    std::size_t getOrder() const;

    /// Return the number of basis functions.
    std::size_t size() const;

    /**
     *  Return a scaled basis that delegates to a copy of `this`.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scaled(Scaling2d const & first) const;

    /// Allocate workspace that can be passed to sumWith() and fill() to avoid repeated memory allocations.
    Workspace makeWorkspace() const;

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
     */
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients) const;

    /// Evaluate a basis expansion with the given coefficients (external workspace version).
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients, Workspace & workspace) const;

    /**
     *  Evaluate the basis at a given point.
     *
     *  @param[in] point   Point at which to evaluate the basis functions.
     *  @param[out] basis  Flattened output vector.
     *                     See Basis1d::fill more information.
     */
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis) const;

    /// Evaluate the basis at a given point (external workspace version).
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis, Workspace & workspace) const;

private:
    Basis1d _basis1d;
};

}}} // namespace lsst::geom::polynomials

#endif // DOXYGEN
#endif // !LSST_AFW_MATH_POLYNOMIALS_Basis2d_h_INCLUDED
