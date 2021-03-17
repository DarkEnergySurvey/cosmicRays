#pragma once

#include "pybind/bind.h"
#include "lsst/utils/python.h"

namespace lsst {
namespace afw {

WRAP(Afw);

namespace typehandling {
WRAP(Typehandling);
}

namespace table {
WRAP(Table);

}

namespace math {
WRAP(Math);
}

namespace image {
WRAP(Initial);
WRAP(Image);
}

namespace geom {
WRAP(Geom);
}
namespace formatters {
WRAP(Formatters);
}

namespace fits {
WRAP(Fits);
}

namespace display {
WRAP(Display);
}

namespace detection {
WRAP(Detection);
}

namespace coord {
WRAP(Coord);
}

namespace cameraGeom {
WRAP(CameraGeom);
}
}
}
