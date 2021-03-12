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
#ifndef LSST_AFW_MATH_POLYNOMIALS_BinomialMatrix_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_BinomialMatrix_h_INCLUDED

#include "Eigen/Core"

namespace lsst { namespace geom { namespace polynomials {

/**
 *  A class that computes binomial coefficients up to a certain power.
 *
 *  The binomial coefficient is defined as:
 *  @f[
 *     \left(\begin{array}{ c }
 *       n \\
 *       k
 *     \end{array}\right)
 *     = \frac{n!}{k!(n-k)!}
 *  @f]
 *  with both @f$n@f$ and @f$k@f$ nonnegative integers and @f$k \le n@f$
 *
 *  This class uses recurrence relations to avoid computing factorials directly,
 *  making it both more efficient and numerically stable.
 */
class BinomialMatrix {
public:

    /**
     *  Construct an object that can compute binomial coefficients with @f$n@f$
     *  up to and including the given value.
     */
    explicit BinomialMatrix(int nMax);

    /**
     *  Return the binomial coefficient.
     *
     *  No error checking is performed; the behavior of this method is
     *  undefined if the given values do not satisfy
     *  @code
     *  n <= nMax && k <= n && n >=0 && k >= 0
     *  @endcode
     */
    double operator()(int n, int k) const {
        return _matrix(n, k);
    }

private:
    Eigen::MatrixXd _matrix;
};

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_BinomialMatrix_h_INCLUDED
