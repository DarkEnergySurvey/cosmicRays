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

/*
 * A boolean pair class used to express the output of spatial predicates on Point and Extent.
 */
#ifndef LSST_GEOM_COORDINATEEXPR_H
#define LSST_GEOM_COORDINATEEXPR_H

#include "lsst/geom/CoordinateBase.h"

namespace lsst {
namespace geom {

/**
 *  A boolean coordinate.
 *
 *  CoordinateExpr is intended to be used as a temporary in coordinate comparisons:
 *
 *      Point2D a(3.5,1.2);
 *      Point2D b(-1.5,4.3);
 *      std::cout << all(a.lt(b)) << std::endl;  // false
 *      std::cout << any(a.lt(b)) << std::endl;  // true
 *
 *  CoordinateExpr is not a true lazy-evaluation expression template, as that seems unnecessary when
 *  the object is typically only two bools large (smaller than the raw pointers necessary to implement
 *  a lazy solution).  The consequence is that there's no short-circuiting of logical operators, but I don't
 *  think that will even remotely matter for most use cases.  The any() and all() functions do support
 *  short-circuiting.
 */
template <int N>
class CoordinateExpr : public CoordinateBase<CoordinateExpr<N>, bool, N> {
    typedef CoordinateBase<CoordinateExpr<N>, bool, N> Super;

public:
    /// Construct a CoordinateExpr with all elements set to the same scalar value.
    explicit CoordinateExpr(bool val = false) noexcept : Super(val) {}

    /// Construct a CoordinateExpr from an Eigen vector.
    template <typename Vector>
    explicit CoordinateExpr(Eigen::MatrixBase<Vector> const& vector) : Super(vector) {}

    CoordinateExpr(CoordinateExpr const&) noexcept = default;
    CoordinateExpr(CoordinateExpr&&) noexcept = default;
    CoordinateExpr& operator=(CoordinateExpr const&) noexcept = default;
    CoordinateExpr& operator=(CoordinateExpr&&) noexcept = default;
    ~CoordinateExpr() noexcept = default;

    /**
     *  @name Logical operators
     *
     *  These operators do not provide interoperability with scalars.
     */
    //@{
    CoordinateExpr and_(CoordinateExpr const& rhs) const noexcept;
    CoordinateExpr or_(CoordinateExpr const& rhs) const noexcept;
    CoordinateExpr not_() const noexcept;
    //@}
};

/// Return true if all elements are true.
template <int N>
inline bool all(CoordinateExpr<N> const& expr) noexcept {
    for (int n = 0; n < N; ++n)
        if (!expr[n]) return false;
    return true;
}

/// Return true if any elements are true.
template <int N>
inline bool any(CoordinateExpr<N> const& expr) noexcept {
    for (int n = 0; n < N; ++n)
        if (expr[n]) return true;
    return false;
}

typedef CoordinateExpr<2> CoordinateExpr2;
typedef CoordinateExpr<3> CoordinateExpr3;

}  // namespace geom
}  // namespace lsst

#endif
