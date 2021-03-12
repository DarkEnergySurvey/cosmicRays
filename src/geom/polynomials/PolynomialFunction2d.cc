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

#include <vector>

#include "lsst/geom/polynomials/PolynomialFunction2d.h"
#include "lsst/geom/polynomials/BinomialMatrix.h"
#include "lsst/geom/polynomials/SafeSum.h"


namespace lsst { namespace geom { namespace polynomials {

namespace {

Eigen::VectorXd computePowers(double x, int n) {
    Eigen::VectorXd r(n + 1);
    r[0] = 1.0;
    for (int i = 1; i <= n; ++i) {
        r[i] = r[i - 1]*x;
    }
    return r;
}

} // anonymous


template <PackingOrder packing>
PolynomialFunction2d<packing> simplified(ScaledPolynomialFunction2d<packing> const & f) {
    auto const & basis = f.getBasis();
    std::vector<SafeSum<double>> sums(basis.size());
    std::size_t const n = basis.getOrder();
    auto rPow = computePowers(basis.getScaling().getX().getScale(), n);
    auto sPow = computePowers(basis.getScaling().getY().getScale(), n);
    auto uPow = computePowers(basis.getScaling().getX().getShift(), n);
    auto vPow = computePowers(basis.getScaling().getY().getShift(), n);
    BinomialMatrix binomial(basis.getNested().getOrder());
    for (auto const & i : basis.getIndices()) {
        for (std::size_t j = 0; j <= i.nx; ++j) {
            double tmp = binomial(i.nx, j)*uPow[j] *
                f[i.flat]*rPow[i.nx]*sPow[i.ny];
            for (std::size_t k = 0; k <= i.ny; ++k) {
                sums[basis.index(i.nx - j, i.ny - k)] +=
                    binomial(i.ny, k)*vPow[k]*tmp;
            }
        }
    }
    Eigen::VectorXd result = Eigen::VectorXd::Zero(basis.size());
    for (std::size_t i = 0; i < basis.size(); ++i) {
        result[i] = static_cast<double>(sums[i]);
    }
    return makeFunction2d(basis.getNested(), result);
}

template PolynomialFunction2d<PackingOrder::XY> simplified(
    ScaledPolynomialFunction2d<PackingOrder::XY> const &
);
template PolynomialFunction2d<PackingOrder::YX> simplified(
    ScaledPolynomialFunction2d<PackingOrder::YX> const &
);

}}} // namespace lsst::geom::polynomials
