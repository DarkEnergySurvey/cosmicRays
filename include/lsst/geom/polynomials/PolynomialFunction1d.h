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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction1d_h_INCLUDED

#include "lsst/geom/polynomials/Function1d.h"
#include "lsst/geom/polynomials/PolynomialBasis1d.h"

namespace lsst { namespace geom { namespace polynomials {

/// A Function1d for standard polynomials.
using PolynomialFunction1d = Function1d<PolynomialBasis1d>;

/// A Function1d for scaled standard polynomials.
using ScaledPolynomialFunction1d = Function1d<ScaledPolynomialBasis1d>;

/**
 *  Calculate the standard polynomial function that is equivalent to a scaled
 *  standard polynomial function.
 *
 *  The returned polynomial will course have different coefficients than the
 *  input one, as these need to account for the scaling without it being
 *  explicitly applied.
 */
PolynomialFunction1d simplified(ScaledPolynomialFunction1d const & f);

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction1d_h_INCLUDED
