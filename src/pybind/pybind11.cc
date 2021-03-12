
#include "lsst/utils/python.h"
#include "lsst/utils/Backtrace.h"
#include "lsst/utils/packaging.h"
#include "lsst/utils/Demangle.h"

#include "pybind/geom_bind.h"
#include "pybind/pex_bind.h"
#include "pybind/daf_bind.h"
#include "pybind/sphgeom_bind.h"
#include "pybind/afw_bind.h"
#include <typeinfo>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "pybind11/eigen.h"
#include "pybind11/numpy.h"


namespace py = pybind11;
using namespace pybind11::literals;
namespace utils = lsst::utils;
namespace python = lsst::utils::python;

/*void wrapBacktrace(utils::python::WrapperCollection & wrappers) {
    wrappers.wrap(
        [](auto & mod) {
            utils::Backtrace &backtrace = utils::Backtrace::get();
            // Trick to tell the compiler backtrace is used and should not be
            // optimized away, as well as convenient way to check if backtrace
            // is enabled.
            mod.def("isEnabled", [&backtrace]() -> bool { return backtrace.isEnabled(); });
        }
    );
}

void wrapPackaging(utils::python::WrapperCollection & wrappers) {
    wrappers.wrap(
        [](auto & mod) {
            mod.def("getPackageDir", utils::getPackageDir);
        }
    );
}

void wrapDemangle(utils::python::WrapperCollection & wrappers) {
    wrappers.wrap(
        [](auto & mod) {
            mod.def("demangleType", utils::demangleType);
        }
    );
}
*/






PYBIND11_MODULE(libcosmicRays, mod) {
    py::print("CALL MOD");

    auto pexbase = mod.def_submodule("pex");
    lsst::pex::wrapPex(pexbase);


    /*utils::python::WrapperCollection wrappers(mod, "_utils");
    {
        auto backtraceWrappers = wrappers.makeSubmodule("backtrace");
        wrapBacktrace(backtraceWrappers);
        wrappers.collectSubmodule(std::move(backtraceWrappers));
    }
    wrapPackaging(wrappers);
    wrapDemangle(wrappers);
    wrappers.finish();
*/


    auto dafbase = mod.def_submodule("daf");
    lsst::daf::wrapDaf(dafbase);
    auto sphgeombase = mod.def_submodule("sphgeom");
    lsst::sphgeom::wrapSphgeom(sphgeombase);
    auto geombase = mod.def_submodule("geom");
    lsst::geom::wrapGeom(geombase);
    auto afwbase = mod.def_submodule("afw");
    lsst::afw::wrapAfw(afwbase);
}

