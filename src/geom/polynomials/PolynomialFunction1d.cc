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

#include "lsst/geom/polynomials/PolynomialFunction1d.h"
#include "lsst/geom/polynomials/BinomialMatrix.h"
#include "lsst/geom/polynomials/SafeSum.h"


namespace lsst { namespace geom { namespace polynomials {

PolynomialFunction1d simplified(ScaledPolynomialFunction1d const & f) {
    auto const & basis = f.getBasis();
    std::vector<SafeSum<double>> sums(basis.size());
    double const s = basis.getScaling().getScale();
    double const v = basis.getScaling().getShift();
    double sn = 1; // s^n
    BinomialMatrix binomial(basis.getNested().getOrder());
    for (std::size_t n = 0; n < basis.size(); ++n, sn *= s) {
        double vk = 1; // v^k
        for (std::size_t k = 0; k <= n; ++k, vk *= v) {
            sums[n - k] += sn*binomial(n, k)*f[n]*vk;
        }
    }
    Eigen::VectorXd result = Eigen::VectorXd::Zero(basis.size());
    for (std::size_t n = 0; n < basis.size(); ++n) {
        result[n] = static_cast<double>(sums[n]);
    }
    return makeFunction1d(basis.getNested(), result);
}

}}} // namespace lsst::geom::polynomials
