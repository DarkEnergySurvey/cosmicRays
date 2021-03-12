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

#include "lsst/geom/polynomials/BinomialMatrix.h"

namespace lsst { namespace geom { namespace polynomials {

BinomialMatrix::BinomialMatrix(int nMax) : _matrix(Eigen::MatrixXd::Zero(nMax + 1, nMax + 1)) {
    _matrix.col(0).setConstant(1);
    _matrix.diagonal().setConstant(1);
    for (int i = 2; i <= nMax; ++i) {
        _matrix(i, 0) = 1.0;
        _matrix(i, i) = 1.0;
        for (int j = 1; j < i; ++j) {
            _matrix(i, j) = _matrix(i - 1, j - 1) *
                (static_cast<double>(i) / static_cast<double>(j));
        }
    }
}

}}} // namespace lsst::geom::polynomials
