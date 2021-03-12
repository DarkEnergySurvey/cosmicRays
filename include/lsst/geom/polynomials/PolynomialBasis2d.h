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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis2d_h_INCLUDED

#include "lsst/geom/polynomials/PolynomialBasis1d.h"
#include "lsst/geom/polynomials/PackedBasis2d.h"
#include "lsst/geom/polynomials/ScaledBasis2d.h"
#include "lsst/geom/polynomials/Scaling2d.h"

namespace lsst { namespace geom { namespace polynomials {

/// A Basis2d for standard polynomials, templated on packing order.
template <PackingOrder packing>
using PolynomialBasis2d = PackedBasis2d<PolynomialBasis1d, packing>;

/// A Basis2d for scaled standard polynomials, templated on packing order.
template <PackingOrder packing>
using ScaledPolynomialBasis2d = ScaledBasis2d<PolynomialBasis2d<packing>>;

/// A Basis2d for standard polynomials, ordered via PackingOrder::XY.
using PolynomialBasis2dXY = PolynomialBasis2d<PackingOrder::XY>;

/// A Basis2d for standard polynomials, ordered via PackingOrder::YX.
using PolynomialBasis2dYX = PolynomialBasis2d<PackingOrder::YX>;

/// A Basis2d for scaled standard polynomials, ordered via PackingOrder::XY.
using ScaledPolynomialBasis2dXY = ScaledPolynomialBasis2d<PackingOrder::XY>;

/// A Basis2d for scaled standard polynomials, ordered via PackingOrder::YX.
using ScaledPolynomialBasis2dYX = ScaledPolynomialBasis2d<PackingOrder::YX>;

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis2d_h_INCLUDED
