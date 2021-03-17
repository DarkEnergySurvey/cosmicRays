#include "pybind/ip_bind.h"
#include "lsst/ip/isr/isr.h"
#include "lsst/ip/isr/applyLookupTable.h"

namespace lsst {
namespace ip {
namespace isr {
namespace {
template <typename PixelT>
static void declareApplyLookupTable(py::module& mod) {
    mod.def("applyLookupTable", &applyLookupTable<PixelT>, "image"_a, "table"_a, "indOffset"_a);
}

template <typename PixelT>
static void declareCountMaskedPixels(py::module& mod, std::string const& suffix) {
    py::class_<CountMaskedPixels<PixelT>, std::shared_ptr<CountMaskedPixels<PixelT>>> cls(
            mod, ("CountMaskedPixels" + suffix).c_str());

    cls.def("reset", &CountMaskedPixels<PixelT>::reset);
    cls.def("apply", &CountMaskedPixels<PixelT>::apply, "image"_a, "bitmask"_a);
    cls.def("getCount", &CountMaskedPixels<PixelT>::getCount);
}

/**
 * Wrap all code in Isr.h for a given template parameter
 *
 * @tparam PixelT  Pixel type; typically `float` or `double` (potentially could also be
 *                  and integer class, but so far we have not needed those)
 * @param mod  pybind11 module to which to add the wrappers.
 * @param[in] suffix  Class name suffix associated with `PixelT`, e.g. "F" for `float` and "D" for `double`
 *
 * Note that the second (function type) template parameter of `fitOverscanImage` is always `double`.
 */
template <typename PixelT>
static void declareAll(py::module& mod, std::string const& suffix) {
    declareCountMaskedPixels<PixelT>(mod, suffix);

    mod.def("maskNans", &maskNans<PixelT>, "maskedImage"_a, "maskVal"_a, "allow"_a = 0);
    mod.def("fitOverscanImage", &fitOverscanImage<PixelT, double>, "overscanFunction"_a, "overscan"_a,
            "stepSize"_a = 1.1, "sigma"_a = 1);
}

WRAP(ApplyLookupTable) {
    declareApplyLookupTable<float>(mod);
    declareApplyLookupTable<double>(mod);
}

WRAP(Isr) {
    declareAll<float>(mod, "F");
    declareAll<double>(mod, "D");
}

}

WRAP(IP) {
    wrapIsr(mod);
    wrapApplyLookupTable(mod);
}

}
}
}
