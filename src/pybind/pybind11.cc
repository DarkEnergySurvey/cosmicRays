
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
    py::print("CALL MOD1");

    auto utilbase = mod.def_submodule("utils");
    lsst::utils::wrapUtils(mod);
    py::print("CALL MOD2");
    auto base = mod.def_submodule("base");
    lsst::base::wrapBase(base);
    py::print("CALL MOD3");
    auto pexbase = mod.def_submodule("pex");
    lsst::pex::wrapPex(pexbase);
    py::print("CALL MOD4");

    auto dafbase = mod.def_submodule("daf");
    lsst::daf::wrapDaf(dafbase);
    py::print("CALL MOD5");
    auto sphgeombase = mod.def_submodule("sphgeom");
    lsst::sphgeom::wrapSphgeom(sphgeombase);
    py::print("CALL MOD6");
    auto geombase = mod.def_submodule("geom");
    lsst::geom::wrapGeom(geombase);
    py::print("CALL MOD7");
    auto afwbase = mod.def_submodule("afw");
    lsst::afw::wrapAfw(afwbase);
    py::print("CALL MOD8");
    auto ipbase = mod.def_submodule("ip");
    lsst::ip::isr::wrapIP(ipbase);
    py::print("CALL MOD9");
    auto measbase = mod.def_submodule("meas");
    lsst::meas::wrapMeas(measbase);
}

