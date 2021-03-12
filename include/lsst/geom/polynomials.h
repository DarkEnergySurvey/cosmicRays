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
#ifndef LSST_AFW_MATH_polynomials_h_INCLUDED
#define LSST_AFW_MATH_polynomials_h_INCLUDED

#include "lsst/geom/polynomials/BinomialMatrix.h"
#include "lsst/geom/polynomials/PolynomialFunction1d.h"
#include "lsst/geom/polynomials/Chebyshev1Function1d.h"
#include "lsst/geom/polynomials/PolynomialFunction2d.h"
#include "lsst/geom/polynomials/Chebyshev1Function2d.h"

namespace lsst { namespace geom {

/**
 *  @namespace lsst::geom::polynomials Low-level polynomials (including special polynomials) in C++.
 *
 *  The geom::polynomials library provides low-level classes for
 *  efficiently evaluating polynomial basis functions and expansions in C++.
 *  The classes here:
 *   - are not available in Python;
 *   - are not polymorphic (no virtual functions);
 *   - provide workspace objects to minimize memory allocations when
 *     appropriate;
 *   - do not throw exceptions (users are responsible for providing valid
 *     inputs);
 *   - have almost no outside dependencies (just Eigen);
 *   - use templates to allow them to work with any array/vector objects (not
 *     just Eigen).
 *
 *  They are intended to be used as the building blocks of higher-level
 *  objects that are visible to Python users and persisted with our data
 *  products, such as afw::math::ChebyshevBoundedField and the
 *  afw::math::Function hierarchy.
 *
 *  At present, the library only includes support for 1-d and 2-d standard
 *  polynomials and Chebyshev polynomials of the first kind, but adding
 *  support for any other function defined by a recurrence relation (i.e. any
 *  other special polynomial) should be extremely easy (see
 *  RecurrenceBasis1d), and need not be done within the polynomials library
 *  itself.
 *
 *  For both 1-d and 2-d, the library contains the following kinds of objects:
 *
 *   - Basis1d and Basis2d: objects that evaluate basis functions at
 *     individual points.  The only concrete Basis1d implementation provided
 *     at present is RecurrenceBasis1d template class, which can be used (with
 *     the recurrence relation as a template parameter) to implement most
 *     special polynomials.  Recurrence relations for standard polynomials
 *     (PolynomialRecurrence) and Chebyshev polynomials of the first kind
 *     (Chebyshev1Recurrence) are provided here.   The PackedBasis2d template
 *     provides an implementation of Basis2d that evaluates a Basis1d over
 *     each dimension.
 *     @ref PolynomialBasis1d, @ref Chebyshev1Basis1d, @ref PolynomialBasis2d,
 *     and @ref Chebyshev1Basis2d provide typedefs to instantiations of these
 *     that should generally be preferred to explicit use of these templates.
 *
 *   - Function1d and Function2d: templates that combine a Basis1d or Basis2d
 *     with a vector of coefficients.
 *
 *   - Scaling1d and Scaling2d: transformation objects that, with the
 *     ScaledBasis1d and ScaledBasis2d templates, allow a basis or function to
 *     be constructed that remaps points as part of evaluating the basis
 *     functions.  For example, a @ref ScaledChebyshev1Basis1d combines both
 *     the a Chebyshev polynomial basis and the scaling from some domain to
 *     [-1, 1] that enables most special properties of Chebyshevs.  When
 *     necessary, the simplified() functions can be used to convert a
 *     ScaledPolynomialFunction1d or ScaledPolynomialFunction2d into an
 *     equivalent PolynomialFunction1d or PolynomialFunction2d by folding the
 *     scaling into the coefficients.
 *
 *  The library also includes a few utility classes and functions:
 *
 *   - BinomialMatrix provides an efficient and stable way to get
 *     binomial coefficients.
 *
 *   - PackedIndexIterator and PackedIndexRange implement PackingOrders that
 *     allow a conceptual 2-d triangular matrix to be implemented on top of
 *     a 1-d array.
 */
namespace polynomials {

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_polynomials_h_INCLUDED
