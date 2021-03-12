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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis2d_h_INCLUDED

#include "lsst/geom/polynomials/Chebyshev1Basis1d.h"
#include "lsst/geom/polynomials/PackedBasis2d.h"
#include "lsst/geom/polynomials/ScaledBasis2d.h"
#include "lsst/geom/polynomials/Scaling2d.h"

namespace lsst { namespace geom { namespace polynomials {

/// A Basis2d for Chebyshev polynomials of the first kind, templated on packing order.
template <PackingOrder packing>
using Chebyshev1Basis2d = PackedBasis2d<Chebyshev1Basis1d, packing>;

/// A Basis2d for scaled Chebyshev polynomials of the first kind, templated on packing order.
template <PackingOrder packing>
using ScaledChebyshev1Basis2d = ScaledBasis2d<Chebyshev1Basis2d<packing>>;

/// A Basis2d for Chebyshev polynomials of the first kind, ordered via PackingOrder::XY.
using Chebyshev1Basis2dXY = Chebyshev1Basis2d<PackingOrder::XY>;

/// A Basis2d for Chebyshev polynomials of the first kind, ordered via PackingOrder::YX.
using Chebyshev1Basis2dYX = Chebyshev1Basis2d<PackingOrder::YX>;

/// A Basis2d for scaled Chebyshev polynomials of the first kind, ordered via PackingOrder::XY.
using ScaledChebyshev1Basis2dXY = ScaledChebyshev1Basis2d<PackingOrder::XY>;

/// A Basis2d for scaled Chebyshev polynomials of the first kind, ordered via PackingOrder::YX.
using ScaledChebyshev1Basis2dYX = ScaledChebyshev1Basis2d<PackingOrder::YX>;

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis2d_h_INCLUDED
