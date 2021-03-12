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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Function1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Function1d_h_INCLUDED

#include "Eigen/Core"
#include "lsst/geom/polynomials/SafeSum.h"
#include "lsst/geom/polynomials/ScaledBasis1d.h"

namespace lsst { namespace geom { namespace polynomials {

/**
 *  A 1-d function defined by a series expansion and its coefficients.
 *
 *  A Function1d combines a Basis1d that defines basis functions @f$B_n(x)@f$
 *  with a vector of associated coefficients @f$a_n@f$.  Evaluating the
 *  function computes
 *  @f[
 *      \sum_{n=0}^{n \le N} a_n B_n(x)
 *  @f]
 */
template <typename Basis_>
class Function1d {
public:

    using iterator = double *;
    using const_iterator = double const *;

    /// The basis type used by this function.
    using Basis = Basis_;

    /// Construct with zero-valued coefficients.
    explicit Function1d(Basis const & basis) :
        _basis(basis),
        _coefficients(Eigen::VectorXd::Zero(basis.size()))
    {}

    /// Construct with coefficients from an Eigen object.
    Function1d(Basis const & basis, Eigen::VectorXd const & coefficients) :
        _basis(basis),
        _coefficients(coefficients)
    {
        assert(basis.size() == static_cast<std::size_t>(_coefficients.size()));
    }

    /// Construct by copying coefficients from an STL iterator range.
    template <typename Iterator>
    Function1d(Basis const & basis, Iterator first, Iterator last) :
        _basis(basis),
        _coefficients(basis.size())
    {
        assert(std::distance(first, last) == static_cast<std::ptrdiff_t>(basis.size()));
        std::copy(first, last, &_coefficients[0]);
    }

    /// Default copy constructor.
    Function1d(Function1d const &) = default;

    /// Default move constructor.
    Function1d(Function1d &&) = default;

    /// Default copy assignment.
    Function1d & operator=(Function1d const &) = default;

    /// Default move assignment.
    Function1d & operator=(Function1d &&) = default;

    //@{
    /// Iterators over coefficients
    iterator begin() { return _coefficients.data(); }
    iterator end() { return begin() + size(); }
    const_iterator cbegin() const { return _coefficients.data(); }
    const_iterator cend() const { return begin() + size(); }
    const_iterator begin() const { return _coefficients.data(); }
    const_iterator end() const { return begin() + size(); }
    //@}

    /// Return the associated Basis1d object.
    Basis const & getBasis() const { return _basis; }

    /// Return the number of coefficients.
    std::size_t size() const { return _basis.size(); }

    /// Evaluate the function at the given point.
    double operator()(double x, SumMode mode=SumMode::FAST) const {
        return _basis.sumWith(x, _coefficients, mode);
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
    Function1d<typename Basis::Scaled> scaled(Scaling1d const & scaling) const {
        return Function1d<typename Basis::Scaled>(getBasis().scaled(scaling), _coefficients);
    }

private:
    Basis _basis;
    Eigen::VectorXd _coefficients;
};

/// Create a Function1d of the appropriate type from a Basis1d and an Eigen object containing coefficients.
template <typename Basis>
Function1d<Basis> makeFunction1d(Basis const & basis, Eigen::VectorXd const & coefficients) {
    return Function1d<Basis>(basis, coefficients);
}

/// Create a Function1d of the appropriate type from a Basis1d and an iterator range to copy coefficients from.
template <typename Basis, typename Iterator>
Function1d<Basis> makeFunction1d(Basis const & basis, Iterator first, Iterator last) {
    return Function1d<Basis>(basis, first, last);
}


}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Function1d_h_INCLUDED
