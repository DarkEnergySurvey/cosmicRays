#pragma once
#include "pybind/bind.h"
#include "lsst/sphgeom/Angle.h"
#include "lsst/sphgeom/AngleInterval.h"
#include "lsst/sphgeom/Box.h"
#include "lsst/sphgeom/Box3d.h"
#include "lsst/sphgeom/Chunker.h"
#include "lsst/sphgeom/Circle.h"
#include "lsst/sphgeom/ConvexPolygon.h"
#include "lsst/sphgeom/Ellipse.h"
#include "lsst/sphgeom/HtmPixelization.h"
#include "lsst/sphgeom/Interval1d.h"
#include "lsst/sphgeom/LonLat.h"
#include "lsst/sphgeom/Matrix3d.h"
#include "lsst/sphgeom/Mq3cPixelization.h"
#include "lsst/sphgeom/NormalizedAngle.h"
#include "lsst/sphgeom/NormalizedAngleInterval.h"
#include "lsst/sphgeom/Pixelization.h"
#include "lsst/sphgeom/Q3cPixelization.h"
#include "lsst/sphgeom/RangeSet.h"
#include "lsst/sphgeom/Region.h"
#include "lsst/sphgeom/UnitVector3d.h"
#include "lsst/sphgeom/Vector3d.h"
#include "lsst/sphgeom/curve.h"
#include "lsst/sphgeom/orientation.h"
#include "lsst/sphgeom/python.h"
#include "lsst/sphgeom/python/interval.h"
#include "lsst/sphgeom/python/relationship.h"
#include "lsst/sphgeom/python/utils.h"
#include "lsst/sphgeom/utils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace sphgeom {

WRAP(Sphgeom);
}
}
