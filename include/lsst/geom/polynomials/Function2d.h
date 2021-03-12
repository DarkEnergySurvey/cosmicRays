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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Function2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Function2d_h_INCLUDED

#include "Eigen/Core"
#include "lsst/geom/polynomials/SafeSum.h"
#include "lsst/geom/polynomials/ScaledBasis2d.h"

namespace lsst { namespace geom { namespace polynomials {

/**
 *  A 2-d function defined by a series expansion and its coefficients.
 *
 *  A Function1d combines a Basis2d that defines basis functions @f$B_{n}(x, y)@f$
 *  with a flattened vector of associated coefficients @f$a_{n}@f$.  The basis
 *  defines the order of coefficients.  Evaluating the function computes
 *  @f[
 *      \sum_{n}^{n \le N} a_{n} B_{n}(x, y)
 *  @f]
 */
template <typename Basis_>
class Function2d {
public:

    using iterator = double *;
    using const_iterator = double const *;

    /// The basis type used by this function.
    using Basis = Basis_;

    /// Type returned by makeWorkspace().
    using Workspace = typename Basis::Workspace;

    /// Construct with zero-valued coefficients.
    explicit Function2d(Basis const & basis) :
        _basis(basis),
        _coefficients(Eigen::VectorXd::Zero(basis.size()))
    {}

    /// Construct with coefficients from an Eigen object.
    explicit Function2d(Basis const & basis, Eigen::VectorXd const & coefficients) :
        _basis(basis),
        _coefficients(coefficients)
    {
        assert(basis.size() == static_cast<std::size_t>(_coefficients.size()));
    }

    /// Construct by copying coefficients from an STL iterator range.
    template <typename Iterator>
    explicit Function2d(Basis const & basis, Iterator first, Iterator last) :
        _basis(basis),
        _coefficients(basis.size())
    {
        assert(std::distance(first, last) == static_cast<std::ptrdiff_t>(basis.size()));
        std::copy(first, last, &_coefficients[0]);
    }

    /// Default copy constructor.
    Function2d(Function2d const &) = default;

    /// Default move constructor.
    Function2d(Function2d &&) = default;

    /// Default copy assignment.
    Function2d & operator=(Function2d const &) = default;

    /// Default move assignment.
    Function2d & operator=(Function2d &&) = default;

    //@{
    /// Iterators over coefficients.
    iterator begin() { return _coefficients.data(); }
    iterator end() { return begin() + size(); }
    const_iterator cbegin() const { return _coefficients.data(); }
    const_iterator cend() const { return begin() + size(); }
    const_iterator begin() const { return _coefficients.data(); }
    const_iterator end() const { return begin() + size(); }
    //@}

    /// Return the associated Basis2d object.
    Basis const & getBasis() const { return _basis; }

    /// Return the number of coefficients.
    std::size_t size() const { return _basis.size(); }

    /// Allocate workspace that can be passed to operator() to avoid repeated memory allocations.
    Workspace makeWorkspace() const { return _basis.makeWorkspace(); }

    /// Evaluate the function at the given point.
    double operator()(geom::Point2D const & point, SumMode mode=SumMode::FAST) const {
        return _basis.sumWith(point, _coefficients, mode);
    }

    /// Evaluate the function at the given point.
    double operator()(geom::Point2D const & point, Workspace & workspace, SumMode mode=SumMode::FAST) const {
        return _basis.sumWith(point, _coefficients, workspace, mode);
    }

    //@{
    /**
     *  Return the coefficient associated with the nth basis function.
     *
     *  Caller is responsible for ensuring that the given index is valid.
     */
    double & operator[](std::size_t n) { return  begin()[n]; }
    double const & operator[](std::size_t n) const { return begin()[n]; }
    //@}

    //@{
    /**
     *  Return the coefficient vector as an Eigen matrix-like object.
     *
     *  The exact type of the returned object is unspecified, but it is
     *  guaranteed to be a view.
     */
    auto getCoefficients() {
        // Return a block view to ensure the caller only modify the values, not the size.
        return _coefficients.head(size());
    }
    auto getCoefficients() const { return _coefficients.head(size()); }
    //@}

    /// Return a new function that applies the given scaling to all points before evaluation.
    Function2d<typename Basis::Scaled> scaled(Scaling2d const & scaling) const {
        return Function2d<typename Basis::Scaled>(getBasis().scaled(scaling), _coefficients);
    }

private:
    Basis _basis;
    Eigen::VectorXd _coefficients;
};

/// Create a Function2d of the appropriate type from a Basis2d and an Eigen object containing coefficients.
template <typename Basis>
Function2d<Basis> makeFunction2d(Basis const & basis, Eigen::VectorXd const & coefficients) {
    return Function2d<Basis>(basis, coefficients);
}

/// Create a Function2d of the appropriate type from a Basis2d and an iterator range to copy coefficients from.
template <typename Basis, typename Iterator>
Function2d<Basis> makeFunction2d(Basis const & basis, Iterator first, Iterator last) {
    return Function2d<Basis>(basis, first, last);
}

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Function2d_h_INCLUDED
