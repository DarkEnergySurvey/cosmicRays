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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis1d_h_INCLUDED

#include "lsst/geom/polynomials/Scaling1d.h"
#include "lsst/geom/polynomials/RecurrenceBasis1d.h"
#include "lsst/geom/polynomials/ScaledBasis1d.h"

namespace lsst { namespace geom { namespace polynomials {

/**
 *  A Recurrence for standard polynomials.
 *
 *  @see RecurrenceBasis1d.
 */
class PolynomialRecurrence {
public:

    static double getB0(double x) noexcept { return 1; }

    static double getB1(double x) noexcept { return x; }

    static double next(double x, std::size_t n, double current, double previous) noexcept {
        return current*x;
    }

};

/// A Basis1d for standard polynomials.
using PolynomialBasis1d = RecurrenceBasis1d<PolynomialRecurrence>;

/// A ScaledBasis1d for standard polynomials.
using ScaledPolynomialBasis1d = ScaledBasis1d<PolynomialBasis1d>;

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis1d_h_INCLUDED
