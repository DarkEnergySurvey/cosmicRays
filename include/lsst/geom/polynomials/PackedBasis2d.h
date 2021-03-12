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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PackedBasis2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PackedBasis2d_h_INCLUDED

#include "lsst/geom/Point.h"
#include "lsst/geom/polynomials/SafeSum.h"
#include "lsst/geom/polynomials/PackedIndex.h"
#include "lsst/geom/polynomials/ScaledBasis2d.h"

namespace lsst { namespace geom { namespace polynomials {

template <typename Basis1d, PackingOrder packing>
class PackedBasis2d;


/**
 *  A workspace object that can be used to avoid extra memory allocations in
 *  repeated calls to PackedBasis2d methods.
 */
class PackedBasisWorkspace2d {
public:

    /// Construct workspace for a basis with the given order.
    explicit PackedBasisWorkspace2d(std::size_t order) : _x(order + 1), _y(order + 1) {}

    /// Return the maximum order this workspace supports.
    std::size_t getOrder() const { return _x.size() - 1; }

private:

    template <typename Recurrence, PackingOrder packing>
    friend class PackedBasis2d;

    Eigen::VectorXd _x;
    Eigen::VectorXd _y;
};

template <typename Basis>
class Function2d;

/**
 *  A Basis2d formed from the product of a Basis1d for each of x and y,
 *  truncated at the sum of their orders.
 *
 *  If @f$B_n(x)@f$ are the basis functions for the nested Basis1d, the basis
 *  functions of a PackedBasis2d with order @f$N@f$ are @f$B_m(x)B_n(y)@f$ for
 *  all combinations with @f$m + n \le N@f$.
 *
 *  The ordering of the products of 1-d basis functions in this 2-d basis is
 *  defined by the PackedIndexRange class, which is accessible as the
 *  PackedBasis2d::IndexRange. Note that while PackedBasis2d uses this ordering,
 *  Basis2d objects in general are not required to.
 */
template <typename Basis1d, PackingOrder packing>
class PackedBasis2d {
public:

    /// A Function2d object that uses this basis.
    using Function = Function2d<PackedBasis2d>;

    /// The type returned by scale().
    using Scaled = ScaledBasis2d<PackedBasis2d>;

    /// The type returned by makeWorkspace().
    using Workspace = PackedBasisWorkspace2d;

    /// The type returned by getIndices().
    using IndexRange = PackedIndexRange<packing>;

    /// Return the size of a PackedBasis with the given order.
    static constexpr std::size_t computeSize(std::size_t order) { return IndexRange::computeSize(order); }

    /// Construct from a 1-d basis that will be used for both x and y.
    explicit PackedBasis2d(Basis1d const & basis1d) : _basis1d(basis1d) {}

    /// Construct by forwarding all arguments to the 1-d basis constructor.
    template <typename ...Args>
    explicit PackedBasis2d(Args&& ...args) : _basis1d(std::forward<Args>(args)...) {}

    /// Default copy constructor.
    PackedBasis2d(PackedBasis2d const &) = default;

    /// Default move constructor.
    PackedBasis2d(PackedBasis2d &&) = default;

    /// Default copy assignment.
    PackedBasis2d & operator=(PackedBasis2d const &) = default;

    /// Default move assignment.
    PackedBasis2d & operator=(PackedBasis2d &&) = default;

    /// Return the maximum order of the basis.
    std::size_t getOrder() const noexcept { return _basis1d.getOrder(); }

    /// Return the number of basis functions.
    std::size_t size() const noexcept{ return IndexRange::computeSize(getOrder()); }

    /**
     *  Return a scaled basis that delegates to a copy of `this`.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scaled(Scaling2d const & first) const {
        return Scaled(*this, first);
    }

    /// Return the flattened index of the basis function with the given x and y orders.
    std::size_t index(std::size_t x, std::size_t y) const {
        return IndexRange::computeIndex(x, y);
    }

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
     *
     *  See PackedIndexIterator for documentation of the actual order.
     */
    IndexRange getIndices() const noexcept {
        return IndexRange(typename IndexRange::iterator(), IndexRange::iterator::makeEnd(getOrder()));
    }

    /// Allocate a workspace that can be passed to sumWith() and fill() to avoid repeated memory allocations.
    Workspace makeWorkspace() const { return Workspace(getOrder());}

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
     *  @param[in] workspace     Workspace object returned by makeWorkspace().
     *  @param[in] mode          Enum indicating the tradeoff to make between
     *                           speed and numerical precision.
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              the same exception safety as it if it does.
     */
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients,
                   Workspace & workspace, SumMode mode=SumMode::FAST) const {
        assert(workspace.getOrder() >= getOrder());
        _basis1d.fill(point.getX(), workspace._x);
        _basis1d.fill(point.getY(), workspace._y);
        // This universal lambda lets us effectively template most of the
        // implementation of this function on double vs. SafeSum<double>
        // without having to define an external template.
        auto accumulate = [coefficients, &workspace, this](auto & sum) {
            for (auto const & index : getIndices()) {
                sum += coefficients[index.flat]*workspace._x[index.nx]*workspace._y[index.ny];
            }
        };
        double result = 0.0;
        if (mode == SumMode::FAST) {
            double z = 0.0;
            accumulate(z);
            result = z;
        } else {
            SafeSum<double> z;
            accumulate(z);
            result = static_cast<double>(z);
        }
        return result;
    }

    /// Evaluate a basis expansion with the given coefficients (internal workspace version).
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients,
                   SumMode mode=SumMode::FAST) const {
        auto workspace = makeWorkspace();
        return sumWith(point, coefficients, workspace, mode);
    }

    /**
     *  Evaluate the basis at a given point.
     *
     *  @param[in] point       Point at which to evaluate the basis functions.
     *  @param[out] basis      Flattened output vector.
     *                         See Basis1d::fill more information.
     *  @param[in] workspace   Workspace object returned by makeWorkspace().
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              basic exception safety if it does.
     */
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis, Workspace & workspace) const {
        assert(workspace.getOrder() >= getOrder());
        _basis1d.fill(point.getX(), workspace._x);
        _basis1d.fill(point.getY(), workspace._y);
        for (auto const & index : getIndices()) {
            std::forward<Vector>(basis)[index.flat] = workspace._x[index.nx]*workspace._y[index.ny];
        }
    }

    /// Evaluate the basis at a given point (internal workspace version).
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis) const {
        auto workspace = makeWorkspace();
        fill(point, std::forward<Vector>(basis), workspace);
    }

private:
    Basis1d _basis1d;
};

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PackedBasis2d_h_INCLUDED
