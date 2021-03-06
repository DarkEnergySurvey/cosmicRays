/*
 * LSST Data Management System
 * Copyright 2014-2015 AURA/LSST.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
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
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

/// \file
/// \brief This file contains the Angle class implementation.

#include "lsst/sphgeom/Angle.h"

#include <cstdio>
#include <ostream>

namespace lsst {
namespace sphgeom {

std::ostream & operator<<(std::ostream & os, Angle const & a) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.17g", a.asRadians());
    return os << buf;
}

double Angle::asDegrees() const { return _rad * DEG_PER_RAD; }

double Angle::asRadians() const { return _rad; }

bool Angle::isNormalized() const { return _rad >= 0.0 && _rad <= 2.0 * PI; }

bool Angle::isNan() const { return std::isnan(_rad); }

}} // namespace lsst::sphgeom
