#include "pybind/afw_bind.h"

namespace lsst {
namespace afw {

void wrapAfw(pybind11::module_ &mod) {
    formatters::wrapFormatters(mod);
    typehandling::wrapTypehandling(mod);
    geom::wrapGeom(mod);
    table::wrapTable(mod);
    fits::wrapFits(mod);

    cameraGeom::wrapCameraGeom(mod);
    detection::wrapDetection(mod);
    coord::wrapCoord(mod);

    image::wrapImage(mod);
    math::wrapMath(mod);

    display::wrapDisplay(mod);
}

}
}
