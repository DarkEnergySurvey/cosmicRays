/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

/*
 * Support routines for
 */
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/table/io/Persistable.hpp"

namespace lsst {
namespace afw {
namespace math {

template <typename ReturnT>
std::vector<double> PolynomialFunction2<ReturnT>::getDFuncDParameters(double x, double y) const {
    std::vector<double> coeffs(this->getNParameters());

    //
    // Go through params order by order, evaluating x^r y^s;  we do this by first evaluating
    // y^s for a complete order, then going through again multiplying by x^r
    //
    int i0 = 0;  // starting index for this order's coefficients
    for (int order = 0; order <= this->_order; ++order) {
        coeffs[i0] = 1;
        double zn = y;  // y^s
        for (int i = 1; i <= order; ++i) {
            coeffs[i0 + i] = zn;
            zn *= y;
        }

        zn = x;  // x^r
        for (int i = order - 1; i >= 0; --i) {
            coeffs[i0 + i] *= zn;
            zn *= x;
        }

        i0 += order + 1;
    }

    assert(i0 == static_cast<int>(coeffs.size()));

    return coeffs;
}

}  // namespace math

/// @cond
#define INSTANTIATE(TYPE)                                                                                 \
    template std::shared_ptr<math::Function<TYPE>>                                                        \
    table::io::PersistableFacade<math::Function<TYPE>>::dynamicCast(                                      \
            std::shared_ptr<table::io::Persistable> const&);                                              \
    template std::shared_ptr<math::Function1<TYPE>>                                                       \
    table::io::PersistableFacade<math::Function1<TYPE>>::dynamicCast(                                     \
            std::shared_ptr<table::io::Persistable> const&);                                              \
    template std::shared_ptr<math::Function2<TYPE>>                                                       \
    table::io::PersistableFacade<math::Function2<TYPE>>::dynamicCast(                                     \
            std::shared_ptr<table::io::Persistable> const&);                                              \
    template std::vector<double> math::PolynomialFunction2<TYPE>::getDFuncDParameters(double x, double y) \
            const

INSTANTIATE(double);
INSTANTIATE(float);
/// @endcond
}  // namespace afw
}  // namespace lsst
