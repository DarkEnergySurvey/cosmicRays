
#include "lsst/utils/python.h"
#include "lsst/utils/Backtrace.h"
#include "lsst/utils/packaging.h"
#include "lsst/utils/Demangle.h"

#include "pybind/utils_bind.h"
#include "pybind/base_bind.h"
#include "pybind/ip_bind.h"
#include "pybind/geom_bind.h"
#include "pybind/pex_bind.h"
#include "pybind/daf_bind.h"
#include "pybind/sphgeom_bind.h"
#include "pybind/afw_bind.h"
#include "pybind/meas_bind.h"
#include <typeinfo>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "pybind11/eigen.h"
#include "pybind11/numpy.h"


namespace py = pybind11;
using namespace pybind11::literals;


PYBIND11_MODULE(libcosmicRays, mod) {

    auto utilbase = mod.def_submodule("utils");
    lsst::utils::wrapUtils(mod);

    auto base = mod.def_submodule("base");
    lsst::base::wrapBase(base);

    auto pexbase = mod.def_submodule("pex");
    lsst::pex::wrapPex(pexbase);


    auto dafbase = mod.def_submodule("daf");
    lsst::daf::wrapDaf(dafbase);

    auto sphgeombase = mod.def_submodule("sphgeom");
    lsst::sphgeom::wrapSphgeom(sphgeombase);

    auto geombase = mod.def_submodule("geom");
    lsst::geom::wrapGeom(geombase);

    auto afwbase = mod.def_submodule("afw");
    lsst::afw::wrapAfw(afwbase);

    auto ipbase = mod.def_submodule("ip");
    lsst::ip::isr::wrapIP(ipbase);

    auto measbase = mod.def_submodule("meas");
    lsst::meas::wrapMeas(measbase);
}

