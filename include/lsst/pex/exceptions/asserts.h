/*
 * This file is part of pex_exceptions.
 *
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

#include "boost/format.hpp"

#include "lsst/pex/exceptions/Exception.h"

/**
 * Check whether the given values are equal, and throw an LSST Exception if they are not.
 *
 * The given message must include two Boost.Format placeholders for the two numbers.
 *
 * For example:
 *
 *     LSST_THROW_IF_NE(3, 4, LengthError, "size of foo (%d) is not equal to size of bar (%d)");
 *
 */
#define LSST_THROW_IF_NE(N1, N2, EXC_CLASS, MSG) \
    if ((N1) != (N2)) throw LSST_EXCEPT(EXC_CLASS, (boost::format(MSG) % (N1) % (N2)).str())
