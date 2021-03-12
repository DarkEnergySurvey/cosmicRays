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

#include "lsst/geom/CoordinateBase.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/Extent.h"

namespace lsst {
namespace geom {

template <typename Derived, typename T, int N>
bool allclose(CoordinateBase<Derived, T, N> const &a, CoordinateBase<Derived, T, N> const &b, T rtol,
              T atol) noexcept(std::is_nothrow_copy_constructible<T>::value
                                       &&std::is_nothrow_copy_assignable<T>::value) {
    Eigen::Array<T, N, 1> diff = (a.asEigen().array() - b.asEigen().array()).abs();
    Eigen::Array<T, N, 1> rhs = (0.5 * (a.asEigen().array() + b.asEigen().array())).abs();
    rhs *= rtol;
    rhs += atol;
    return (diff <= rhs).all();
}

template bool allclose<Point2D, double, 2>(CoordinateBase<Point2D, double, 2> const &,
                                           CoordinateBase<Point2D, double, 2> const &, double, double);
template bool allclose<Point3D, double, 3>(CoordinateBase<Point3D, double, 3> const &,
                                           CoordinateBase<Point3D, double, 3> const &, double, double);
template bool allclose<Extent2D, double, 2>(CoordinateBase<Extent2D, double, 2> const &,
                                            CoordinateBase<Extent2D, double, 2> const &, double, double);
template bool allclose<Extent3D, double, 3>(CoordinateBase<Extent3D, double, 3> const &,
                                            CoordinateBase<Extent3D, double, 3> const &, double, double);

}  // namespace geom
}  // namespace lsst
