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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction2d_h_INCLUDED

#include "lsst/geom/polynomials/Function2d.h"
#include "lsst/geom/polynomials/PolynomialBasis2d.h"

namespace lsst { namespace geom { namespace polynomials {

/// A Function2d for standard polynomials.
template <PackingOrder packing>
using PolynomialFunction2d = Function2d<PolynomialBasis2d<packing>>;

/// A Function2d for scaled standard polynomials.
template <PackingOrder packing>
using ScaledPolynomialFunction2d = Function2d<ScaledPolynomialBasis2d<packing>>;

/// A Function2d for standard polynomials, ordered via PackingOrder::XY.
using PolynomialFunction2dXY = PolynomialFunction2d<PackingOrder::XY>;

/// A Function2d for standard polynomials, ordered via PackingOrder::YX.
using PolynomialFunction2dYX = PolynomialFunction2d<PackingOrder::YX>;

/// A Function2d for scaled standard polynomials, ordered via PackingOrder::XY.
using ScaledPolynomialFunction2dXY = ScaledPolynomialFunction2d<PackingOrder::XY>;

/// A Function2d for scaled standard polynomials, ordered via PackingOrder::YX.
using ScaledPolynomialFunction2dYX = ScaledPolynomialFunction2d<PackingOrder::YX>;

/**
 *  Calculate the standard polynomial function that is equivalent to a scaled
 *  standard polynomial function.
 *
 *  The coefficients of the returned polynomial will be different from those
 *  of the input in order to fold in the scaling.  This is primarily useful
 *  in contexts where external code does not support the (more numerically
 *  stable) scaled representation, such as the FITS WCS SIP convention.
 */
template <PackingOrder packing>
PolynomialFunction2d<packing> simplified(ScaledPolynomialFunction2d<packing> const & f);

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction2d_h_INCLUDED
