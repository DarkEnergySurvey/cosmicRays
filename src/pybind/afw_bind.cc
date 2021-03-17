#include "pybind/afw_bind.h"

namespace lsst {
namespace afw {

void wrapAfw(pybind11::module_ &mod) {
    auto fmod = mod.def_submodule("formatters");
    formatters::wrapFormatters(fmod);
    auto tmod = mod.def_submodule("table");
    table::wrapTable(tmod);
    auto tymod = mod.def_submodule("typehandling");
    typehandling::wrapTypehandling(tymod);
    auto gmod = mod.def_submodule("geom");
    geom::wrapGeom(gmod);
    auto fimod = mod.def_submodule("fits");
    fits::wrapFits(fimod);

    auto imod = mod.def_submodule("image");
    image::wrapInitial(imod);
    auto camod = mod.def_submodule("cameraGeom");
    cameraGeom::wrapCameraGeom(camod);
    auto comod = mod.def_submodule("coord");
    coord::wrapCoord(comod);
    auto dmod = mod.def_submodule("detection");
    detection::wrapDetection(dmod);

    image::wrapImage(imod);

    auto mmod = mod.def_submodule("math");
    math::wrapMath(mmod);

    auto dimod = mod.def_submodule("display");
    display::wrapDisplay(dimod);
}

}
}
