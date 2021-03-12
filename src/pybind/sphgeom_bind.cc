#include "pybind/sphgeom_bind.h"

namespace lsst {

namespace sphgeom {
namespace {
py::bytes encode(sphgeom::Region const &self) {
    std::vector<uint8_t> bytes = self.encode();
    return py::bytes(reinterpret_cast<char const *>(bytes.data()),
                     bytes.size());
}
uint64_t _uint64(py::handle const &obj) {
    try {
        return obj.cast<uint64_t>();
    } catch (py::cast_error const &) {
        throw py::value_error(
                    "RangeSet elements and range beginning and "
                    "end points must be non-negative integers "
                    "less than 2**64");
    }
}

/// Make a RangeSet from an iterable. Each item must be an integer that fits
/// in a uint64_t, or a sequence of two such integers.
sphgeom::RangeSet makeRangeSet(py::iterable iterable) {
    sphgeom::RangeSet rs;
    for (py::handle item : iterable) {
        PyObject *o = item.ptr();
        if (PySequence_Check(o) && PySequence_Size(o) == 2) {
            uint64_t first = _uint64(py::reinterpret_steal<py::object>(
                    PySequence_GetItem(o, 0)));
            uint64_t last = _uint64(py::reinterpret_steal<py::object>(
                    PySequence_GetItem(o, 1)));
            rs.insert(first, last);
        } else {
            rs.insert(_uint64(item));
        }
    }
    return rs;
}

std::unique_ptr<sphgeom::Ellipse> decode(py::bytes bytes) {
    uint8_t const *buffer = reinterpret_cast<uint8_t const *>(
            PYBIND11_BYTES_AS_STRING(bytes.ptr()));
    size_t n = static_cast<size_t>(PYBIND11_BYTES_SIZE(bytes.ptr()));
    return sphgeom::Ellipse::decode(buffer, n);
}

/// Make a python list of the ranges in the given RangeSet.
py::list ranges(sphgeom::RangeSet const &self) {
    py::list list;
    for (auto t : self) {
        list.append(py::make_tuple(py::int_(std::get<0>(t)),
                                   py::int_(std::get<1>(t))));
    }
    return list;
}

sphgeom::Vector3d getRow(sphgeom::Matrix3d const &self, py::int_ row) {
    return self.getRow(static_cast<int>(sphgeom::python::convertIndex(3, row)));
}

std::unique_ptr<sphgeom::ConvexPolygon> condecode(py::bytes bytes) {
    uint8_t const *buffer = reinterpret_cast<uint8_t const *>(
            PYBIND11_BYTES_AS_STRING(bytes.ptr()));
    size_t n = static_cast<size_t>(PYBIND11_BYTES_SIZE(bytes.ptr()));
    return sphgeom::ConvexPolygon::decode(buffer, n);
}

std::unique_ptr<sphgeom::Box> boxdecode(py::bytes bytes) {
    uint8_t const *buffer = reinterpret_cast<uint8_t const *>(
            PYBIND11_BYTES_AS_STRING(bytes.ptr()));
    size_t n = static_cast<size_t>(PYBIND11_BYTES_SIZE(bytes.ptr()));
    return sphgeom::Box::decode(buffer, n);
}

std::unique_ptr<sphgeom::Circle> circdecode(py::bytes bytes) {
    uint8_t const *buffer = reinterpret_cast<uint8_t const *>(
            PYBIND11_BYTES_AS_STRING(bytes.ptr()));
    size_t n = static_cast<size_t>(PYBIND11_BYTES_SIZE(bytes.ptr()));
    return sphgeom::Circle::decode(buffer, n);
}

py::str toString(sphgeom::Chunker const &self) {
    return py::str("Chunker({!s}, {!s})")
            .format(self.getNumStripes(), self.getNumSubStripesPerStripe());
}

}

template <>
void defineClass(py::class_<Chunker, std::shared_ptr<Chunker>> &cls) {
    cls.def(py::init<int32_t, int32_t>(), "numStripes"_a,
            "numSubStripesPerStripe"_a);

    cls.def("__eq__", &Chunker::operator==, py::is_operator());
    cls.def("__ne__", &Chunker::operator!=, py::is_operator());

    cls.def_property_readonly("numStripes", &Chunker::getNumStripes);
    cls.def_property_readonly("numSubStripesPerStripe",
                              &Chunker::getNumSubStripesPerStripe);

    cls.def("getChunksIntersecting", &Chunker::getChunksIntersecting,
            "region"_a);
    cls.def("getSubChunksIntersecting",
            [](Chunker const &self, Region const &region) {
                py::list results;
                for (auto const &sc : self.getSubChunksIntersecting(region)) {
                    results.append(py::make_tuple(sc.chunkId, sc.subChunkIds));
                }
                return results;
            },
            "region"_a);
    cls.def("getAllChunks", &Chunker::getAllChunks);
    cls.def("getAllSubChunks", &Chunker::getAllSubChunks, "chunkId"_a);

    cls.def("getChunkBoundingBox", &Chunker::getChunkBoundingBox, "stripe"_a, "chunk"_a);
    cls.def("getSubChunkBoundingBox", &Chunker::getSubChunkBoundingBox, "subStripe"_a, "subChunk"_a);

    cls.def("getStripe", &Chunker::getStripe, "chunkId"_a);
    cls.def("getChunk", &Chunker::getChunk, "chunkId"_a, "stripe"_a);


    cls.def("__str__", &toString);
    cls.def("__repr__", &toString);

    cls.def("__reduce__", [cls](Chunker const &self) {
        return py::make_tuple(cls,
                              py::make_tuple(self.getNumStripes(),
                                             self.getNumSubStripesPerStripe()));
    });
}

template <>
void defineClass(py::class_<Circle, std::unique_ptr<Circle>, Region> &cls) {
    cls.attr("TYPE_CODE") = py::int_(Circle::TYPE_CODE);

    cls.def_static("empty", &Circle::empty);
    cls.def_static("full", &Circle::full);
    cls.def_static("squaredChordLengthFor", &Circle::squaredChordLengthFor,
                   "openingAngle"_a);
    cls.def_static("openingAngleFor", &Circle::openingAngleFor,
                   "squaredChordLength"_a);

    cls.def(py::init<>());
    cls.def(py::init<UnitVector3d const &>(), "center"_a);
    cls.def(py::init<UnitVector3d const &, Angle>(), "center"_a, "angle"_a);
    cls.def(py::init<UnitVector3d const &, double>(), "center"_a,
            "squaredChordLength"_a);
    cls.def(py::init<Circle const &>(), "circle"_a);

    cls.def("__eq__", &Circle::operator==, py::is_operator());
    cls.def("__ne__", &Circle::operator!=, py::is_operator());
    cls.def("__contains__",
            (bool (Circle::*)(Circle const &) const) & Circle::contains,
            py::is_operator());
    // Rewrap this base class method since there are overloads in this subclass
    cls.def("__contains__",
            (bool (Circle::*)(UnitVector3d const &) const) & Circle::contains,
            py::is_operator());

    cls.def("isEmpty", &Circle::isEmpty);
    cls.def("isFull", &Circle::isFull);
    cls.def("getCenter", &Circle::getCenter);
    cls.def("getSquaredChordLength", &Circle::getSquaredChordLength);
    cls.def("getOpeningAngle", &Circle::getOpeningAngle);
    cls.def("contains",
            (bool (Circle::*)(Circle const &) const) & Circle::contains);
    // Rewrap this base class method since there are overloads in this subclass
    cls.def("contains",
            (bool (Circle::*)(UnitVector3d const &) const) & Circle::contains);

    cls.def("isDisjointFrom",
            (bool (Circle::*)(UnitVector3d const &) const) &
                    Circle::isDisjointFrom);
    cls.def("isDisjointFrom",
            (bool (Circle::*)(Circle const &) const) & Circle::isDisjointFrom);
    cls.def("intersects",
            (bool (Circle::*)(UnitVector3d const &) const) &
                    Circle::intersects);
    cls.def("intersects",
            (bool (Circle::*)(Circle const &) const) & Circle::intersects);
    cls.def("isWithin",
            (bool (Circle::*)(UnitVector3d const &) const) & Circle::isWithin);
    cls.def("isWithin",
            (bool (Circle::*)(Circle const &) const) & Circle::isWithin);
    cls.def("clipTo",
            (Circle & (Circle::*)(UnitVector3d const &)) & Circle::clipTo);
    cls.def("clipTo", (Circle & (Circle::*)(Circle const &)) & Circle::clipTo);
    cls.def("clippedTo",
            (Circle(Circle::*)(UnitVector3d const &) const) &
                    Circle::clippedTo);
    cls.def("clippedTo",
            (Circle(Circle::*)(Circle const &) const) & Circle::clippedTo);
    cls.def("expandTo",
            (Circle & (Circle::*)(UnitVector3d const &)) & Circle::expandTo);
    cls.def("expandTo",
            (Circle & (Circle::*)(Circle const &)) & Circle::expandTo);
    cls.def("expandedTo",
            (Circle(Circle::*)(UnitVector3d const &) const) &
                    Circle::expandedTo);
    cls.def("expandedTo",
            (Circle(Circle::*)(Circle const &) const) & Circle::expandedTo);
    cls.def("dilateBy", &Circle::dilateBy, "radius"_a);
    cls.def("dilatedBy", &Circle::dilatedBy, "radius"_a);
    cls.def("erodeBy", &Circle::erodeBy, "radius"_a);
    cls.def("erodedBy", &Circle::erodedBy, "radius"_a);
    cls.def("getArea", &Circle::getArea);
    cls.def("complement", &Circle::complement);
    cls.def("complemented", &Circle::complemented);

    // Note that the Region interface has already been wrapped.

    // The lambda is necessary for now; returning the unique pointer
    // directly leads to incorrect results and crashes.
    cls.def_static("decode",
                   [](py::bytes bytes) { return circdecode(bytes).release(); },
                   "bytes"_a);

    cls.def("__str__", [](Circle const &self) {
        return py::str("Circle({!s}, {!s})")
                .format(self.getCenter(), self.getOpeningAngle());
    });
    cls.def("__repr__", [](Circle const &self) {
        return py::str("Circle({!r}, {!r})")
                .format(self.getCenter(), self.getOpeningAngle());
    });
    cls.def(py::pickle(
            [](const Circle &self) { return python::encode(self); },
            [](py::bytes bytes) { return circdecode(bytes).release(); }));
}

template <>
void defineClass(py::class_<Angle> &cls) {
    cls.def_static("nan", &Angle::nan);
    cls.def_static("fromDegrees", &Angle::fromDegrees);
    cls.def_static("fromRadians", &Angle::fromRadians);

    cls.def(py::init<>());
    cls.def(py::init<double>(), "radians"_a);
    cls.def(py::init<Angle>(), "angle"_a);
    // Construct an Angle from a NormalizedAngle, enabling implicit
    // conversion from NormalizedAngle to Angle in python via
    // py::implicitly_convertible
    cls.def(py::init(
            [](NormalizedAngle &a) {
                return new Angle(a.asRadians());
            }),
            "normalizedAngle"_a);

    cls.def("__eq__", &Angle::operator==, py::is_operator());
    cls.def("__ne__", &Angle::operator!=, py::is_operator());
    cls.def("__lt__", &Angle::operator<, py::is_operator());
    cls.def("__gt__", &Angle::operator>, py::is_operator());
    cls.def("__le__", &Angle::operator<=, py::is_operator());
    cls.def("__ge__", &Angle::operator>=, py::is_operator());

    cls.def("__neg__", (Angle(Angle::*)() const) & Angle::operator-);
    cls.def("__add__", &Angle::operator+, py::is_operator());
    cls.def("__sub__",
            (Angle(Angle::*)(Angle const &) const) & Angle::operator-,
            py::is_operator());
    cls.def("__mul__", &Angle::operator*, py::is_operator());
    cls.def("__rmul__", &Angle::operator*, py::is_operator());
    cls.def("__truediv__", (Angle(Angle::*)(double) const) & Angle::operator/,
            py::is_operator());
    cls.def("__truediv__",
            (double (Angle::*)(Angle const &) const) & Angle::operator/,
            py::is_operator());

    cls.def("__iadd__", &Angle::operator+=);
    cls.def("__isub__", &Angle::operator-=);
    cls.def("__imul__", &Angle::operator*=);
    cls.def("__itruediv__", &Angle::operator/=);

    cls.def("asDegrees", &Angle::asDegrees);
    cls.def("asRadians", &Angle::asRadians);
    cls.def("isNormalized", &Angle::isNormalized);
    cls.def("isNan", &Angle::isNan);

    cls.def("__str__", [](Angle const &self) {
        return py::str("{!s}").format(self.asRadians());
    });
    cls.def("__repr__", [](Angle const &self) {
        return py::str("Angle({!r})").format(self.asRadians());
    });

    cls.def("__reduce__", [cls](Angle const &self) {
        return py::make_tuple(cls, py::make_tuple(self.asRadians()));
    });
}

template <>
void defineClass(
        py::class_<AngleInterval, std::shared_ptr<AngleInterval>> &cls) {
    python::defineInterval<decltype(cls), AngleInterval, Angle>(cls);

    cls.def_static("fromDegrees", &AngleInterval::fromDegrees, "x"_a, "y"_a);
    cls.def_static("fromRadians", &AngleInterval::fromRadians, "x"_a, "y"_a);
    cls.def_static("empty", &AngleInterval::empty);
    cls.def_static("full", &AngleInterval::full);

    cls.def(py::init<>());
    cls.def(py::init<Angle>(), "x"_a);
    cls.def(py::init<Angle, Angle>(), "x"_a, "y"_a);
    cls.def(py::init<AngleInterval const &>(), "interval"_a);

    cls.def("__str__", [](AngleInterval const &self) {
        return py::str("[{!s}, {!s}]")
                .format(self.getA().asRadians(), self.getB().asRadians());
    });
    cls.def("__repr__", [](AngleInterval const &self) {
        return py::str("AngleInterval.fromRadians({!r}, {!r})")
                .format(self.getA().asRadians(), self.getB().asRadians());
    });
}


template <>
void defineClass(py::class_<Box3d, std::shared_ptr<Box3d>> &cls) {
    cls.def_static("empty", &Box3d::empty);
    cls.def_static("full", &Box3d::full);
    cls.def_static("aroundUnitSphere", &Box3d::aroundUnitSphere);

    cls.def(py::init<>());
    cls.def(py::init<Vector3d const &>(), "vector"_a);
    cls.def(py::init<Vector3d const &, Vector3d const &>(), "vector1"_a,
            "vector2"_a);
    cls.def(py::init<Vector3d const &, double, double, double>(), "center"_a,
            "halfWidth"_a, "halfHeight"_a, "halfDepth"_a);
    cls.def(py::init<Interval1d const &, Interval1d const &,
                     Interval1d const &>(),
            "x"_a, "y"_a, "z"_a);
    cls.def(py::init<Box3d const &>(), "box3d"_a);

    cls.def("__eq__",
            (bool (Box3d::*)(Box3d const &) const) & Box3d::operator==,
            py::is_operator());
    cls.def("__eq__",
            (bool (Box3d::*)(Vector3d const &) const) & Box3d::operator==,
            py::is_operator());
    cls.def("__ne__",
            (bool (Box3d::*)(Box3d const &) const) & Box3d::operator!=,
            py::is_operator());
    cls.def("__ne__",
            (bool (Box3d::*)(Vector3d const &) const) & Box3d::operator!=,
            py::is_operator());
    cls.def("__contains__",
            (bool (Box3d::*)(Vector3d const &) const) & Box3d::contains,
            py::is_operator());
    cls.def("__contains__",
            (bool (Box3d::*)(Box3d const &) const) & Box3d::contains,
            py::is_operator());
    cls.def("__len__", [](Box3d const &self) { return py::int_(3); });
    cls.def("__getitem__", [](Box3d const &self, py::int_ row) {
        return self(static_cast<int>(python::convertIndex(3, row)));
    });

    cls.def("x", &Box3d::x);
    cls.def("y", &Box3d::y);
    cls.def("z", &Box3d::z);
    cls.def("isEmpty", &Box3d::isEmpty);
    cls.def("isFull", &Box3d::isFull);
    cls.def("getCenter", &Box3d::getCenter);
    cls.def("getWidth", &Box3d::getWidth);
    cls.def("getHeight", &Box3d::getHeight);
    cls.def("getDepth", &Box3d::getDepth);

    cls.def("contains",
            (bool (Box3d::*)(Vector3d const &) const) & Box3d::contains);
    cls.def("contains",
            (bool (Box3d::*)(Box3d const &) const) & Box3d::contains);
    cls.def("isDisjointFrom",
            (bool (Box3d::*)(Vector3d const &) const) & Box3d::isDisjointFrom);
    cls.def("isDisjointFrom",
            (bool (Box3d::*)(Box3d const &) const) & Box3d::isDisjointFrom);
    cls.def("intersects",
            (bool (Box3d::*)(Vector3d const &) const) & Box3d::intersects);
    cls.def("intersects",
            (bool (Box3d::*)(Box3d const &) const) & Box3d::intersects);
    cls.def("isWithin",
            (bool (Box3d::*)(Vector3d const &) const) & Box3d::isWithin);
    cls.def("isWithin",
            (bool (Box3d::*)(Box3d const &) const) & Box3d::isWithin);

    cls.def("clipTo", (Box3d & (Box3d::*)(Vector3d const &)) & Box3d::clipTo);
    cls.def("clipTo", (Box3d & (Box3d::*)(Box3d const &)) & Box3d::clipTo);
    cls.def("clippedTo",
            (Box3d(Box3d::*)(Vector3d const &) const) & Box3d::clippedTo);
    cls.def("clippedTo",
            (Box3d(Box3d::*)(Box3d const &) const) & Box3d::clippedTo);
    cls.def("expandTo",
            (Box3d & (Box3d::*)(Vector3d const &)) & Box3d::expandTo);
    cls.def("expandTo", (Box3d & (Box3d::*)(Box3d const &)) & Box3d::expandTo);
    cls.def("expandedTo",
            (Box3d(Box3d::*)(Vector3d const &) const) & Box3d::expandedTo);
    cls.def("expandedTo",
            (Box3d(Box3d::*)(Box3d const &) const) & Box3d::expandedTo);

    cls.def("dilateBy", (Box3d & (Box3d::*)(double)) & Box3d::dilateBy,
            "radius"_a);
    cls.def("dilateBy",
            (Box3d & (Box3d::*)(double, double, double)) & Box3d::dilateBy,
            "width"_a, "height"_a, "depth"_a);
    cls.def("dilatedBy", (Box3d(Box3d::*)(double) const) & Box3d::dilatedBy,
            "radius"_a);
    cls.def("dilatedBy",
            (Box3d(Box3d::*)(double, double, double) const) & Box3d::dilatedBy,
            "width"_a, "height"_a, "depth"_a);
    cls.def("erodeBy", (Box3d & (Box3d::*)(double)) & Box3d::erodeBy,
            "radius"_a);
    cls.def("erodeBy",
            (Box3d & (Box3d::*)(double, double, double)) & Box3d::erodeBy,
            "width"_a, "height"_a, "depth"_a);
    cls.def("erodedBy", (Box3d(Box3d::*)(double) const) & Box3d::erodedBy,
            "radius"_a);
    cls.def("erodedBy",
            (Box3d(Box3d::*)(double, double, double) const) & Box3d::erodedBy,
            "width"_a, "height"_a, "depth"_a);

    cls.def("relate",
            (Relationship(Box3d::*)(Vector3d const &) const) & Box3d::relate);
    cls.def("relate",
            (Relationship(Box3d::*)(Box3d const &) const) & Box3d::relate);

    cls.def("__str__", [](Box3d const &self) {
        return py::str("[{!s},\n"
                       " {!s},\n"
                       " {!s}]")
                .format(self.x(), self.y(), self.z());
    });
    cls.def("__repr__", [](Box3d const &self) {
        return py::str("Box3d({!r},\n"
                       "      {!r},\n"
                       "      {!r})")
                .format(self.x(), self.y(), self.z());
    });
    cls.def("__reduce__", [cls](Box3d const &self) {
        return py::make_tuple(cls,
                              py::make_tuple(self.x(), self.y(), self.z()));
    });
}

template <>
void defineClass(py::class_<Box, std::unique_ptr<Box>, Region> &cls) {
    cls.attr("TYPE_CODE") = py::int_(Box::TYPE_CODE);

    cls.def_static("fromDegrees", &Box::fromDegrees, "lon1"_a, "lat1"_a,
                   "lon2"_a, "lat2"_a);
    cls.def_static("fromRadians", &Box::fromRadians, "lon1"_a, "lat1"_a,
                   "lon2"_a, "lat2"_a);
    cls.def_static("empty", &Box::empty);
    cls.def_static("full", &Box::full);
    cls.def_static("halfWidthForCircle", &Box::halfWidthForCircle, "radius"_a,
                   "lat"_a);
    cls.def_static("allLongitudes", &Box::allLongitudes);
    cls.def_static("allLatitudes", &Box::allLatitudes);

    cls.def(py::init<>());
    cls.def(py::init<LonLat const &>(), "point"_a);
    cls.def(py::init<LonLat const &, LonLat const &>(), "point1"_a, "point2"_a);
    cls.def(py::init<LonLat const &, Angle, Angle>(), "center"_a, "width"_a,
            "height"_a);
    cls.def(py::init<NormalizedAngleInterval const &, AngleInterval const &>(),
            "lon"_a, "lat"_a);
    cls.def(py::init<Box const &>(), "box"_a);

    cls.def("__eq__", (bool (Box::*)(Box const &) const) & Box::operator==,
            py::is_operator());
    cls.def("__eq__", (bool (Box::*)(LonLat const &) const) & Box::operator==,
            py::is_operator());
    cls.def("__ne__", (bool (Box::*)(Box const &) const) & Box::operator!=,
            py::is_operator());
    cls.def("__ne__", (bool (Box::*)(LonLat const &) const) & Box::operator!=,
            py::is_operator());
    cls.def("__contains__",
            (bool (Box::*)(LonLat const &) const) & Box::contains,
            py::is_operator());
    cls.def("__contains__", (bool (Box::*)(Box const &) const) & Box::contains,
            py::is_operator());
    // Rewrap this base class method since there are overloads in this subclass
    cls.def("__contains__",
            (bool (Box::*)(UnitVector3d const &) const) & Box::contains,
            py::is_operator());

    cls.def("getLon", &Box::getLon);
    cls.def("getLat", &Box::getLat);
    cls.def("isEmpty", &Box::isEmpty);
    cls.def("isFull", &Box::isFull);
    cls.def("getCenter", &Box::getCenter);
    cls.def("getWidth", &Box::getWidth);
    cls.def("getHeight", &Box::getHeight);
    cls.def("contains", (bool (Box::*)(LonLat const &) const) & Box::contains);
    cls.def("contains", (bool (Box::*)(Box const &) const) & Box::contains);
    // Rewrap this base class method since there are overloads in this subclass
    cls.def("contains",
            (bool (Box::*)(UnitVector3d const &) const) & Box::contains);
    cls.def("isDisjointFrom",
            (bool (Box::*)(LonLat const &) const) & Box::isDisjointFrom);
    cls.def("isDisjointFrom",
            (bool (Box::*)(Box const &) const) & Box::isDisjointFrom);
    cls.def("intersects",
            (bool (Box::*)(LonLat const &) const) & Box::intersects);
    cls.def("intersects", (bool (Box::*)(Box const &) const) & Box::intersects);
    cls.def("isWithin", (bool (Box::*)(LonLat const &) const) & Box::isWithin);
    cls.def("isWithin", (bool (Box::*)(Box const &) const) & Box::isWithin);
    cls.def("clipTo", (Box & (Box::*)(LonLat const &)) & Box::clipTo);
    cls.def("clipTo", (Box & (Box::*)(Box const &)) & Box::clipTo);
    cls.def("clippedTo", (Box(Box::*)(LonLat const &) const) & Box::clippedTo);
    cls.def("clippedTo", (Box(Box::*)(Box const &) const) & Box::clippedTo);
    cls.def("expandTo", (Box & (Box::*)(LonLat const &)) & Box::expandTo);
    cls.def("expandTo", (Box & (Box::*)(Box const &)) & Box::expandTo);
    cls.def("expandedTo",
            (Box(Box::*)(LonLat const &) const) & Box::expandedTo);
    cls.def("expandedTo", (Box(Box::*)(Box const &) const) & Box::expandedTo);
    cls.def("dilateBy", (Box & (Box::*)(Angle)) & Box::dilateBy, "angle"_a);
    cls.def("dilateBy", (Box & (Box::*)(Angle, Angle)) & Box::dilateBy,
            "width"_a, "height"_a);
    cls.def("dilatedBy", (Box(Box::*)(Angle) const) & Box::dilatedBy,
            "angle"_a);
    cls.def("dilatedBy", (Box(Box::*)(Angle, Angle) const) & Box::dilatedBy,
            "width"_a, "height"_a);
    cls.def("erodeBy", (Box & (Box::*)(Angle)) & Box::erodeBy, "angle"_a);
    cls.def("erodeBy", (Box & (Box::*)(Angle, Angle)) & Box::erodeBy, "width"_a,
            "height"_a);
    cls.def("erodedBy", (Box(Box::*)(Angle) const) & Box::erodedBy, "angle"_a);
    cls.def("erodedBy", (Box(Box::*)(Angle, Angle) const) & Box::erodedBy,
            "width"_a, "height"_a);
    cls.def("getArea", &Box::getArea);
    cls.def("relate",
            (Relationship(Box::*)(LonLat const &) const) & Box::relate,
            "point"_a);
    // Rewrap this base class method since there are overloads in this subclass
    cls.def("relate",
            (Relationship(Box::*)(Region const &) const) & Box::relate,
            "region"_a);

    // Note that the Region interface has already been wrapped.

    // The lambda is necessary for now; returning the unique pointer
    // directly leads to incorrect results and crashes.
    cls.def_static("decode",
                   [](py::bytes bytes) { return boxdecode(bytes).release(); },
                   "bytes"_a);

    cls.def("__str__", [](Box const &self) {
        return py::str("Box({!s}, {!s})").format(self.getLon(), self.getLat());
    });
    cls.def("__repr__", [](Box const &self) {
        return py::str("Box({!r}, {!r})").format(self.getLon(), self.getLat());
    });
    cls.def(py::pickle(
            [](const Box &self) { return python::encode(self); },
            [](py::bytes bytes) { return boxdecode(bytes).release(); }));
}

template <>
void defineClass(py::class_<ConvexPolygon, std::unique_ptr<ConvexPolygon>,
                            Region> &cls) {
    cls.attr("TYPE_CODE") = py::int_(ConvexPolygon::TYPE_CODE);

    cls.def_static("convexHull", &ConvexPolygon::convexHull, "points"_a);

    cls.def(py::init<std::vector<UnitVector3d> const &>(), "points"_a);
    // Do not wrap the two unsafe (3 and 4 vertex) constructors
    cls.def(py::init<ConvexPolygon const &>(), "convexPolygon"_a);

    cls.def("__eq__", &ConvexPolygon::operator==, py::is_operator());
    cls.def("__ne__", &ConvexPolygon::operator!=, py::is_operator());

    cls.def("getVertices", &ConvexPolygon::getVertices);
    cls.def("getCentroid", &ConvexPolygon::getCentroid);

    // Note that much of the Region interface has already been wrapped. Here are bits that have not:
    cls.def("contains", py::overload_cast<UnitVector3d const &>(&ConvexPolygon::contains, py::const_));
    cls.def("contains", py::overload_cast<Region const &>(&ConvexPolygon::contains, py::const_));
    cls.def("isDisjointFrom", &ConvexPolygon::isDisjointFrom);
    cls.def("intersects", &ConvexPolygon::intersects);
    cls.def("isWithin", &ConvexPolygon::isWithin);

    // The lambda is necessary for now; returning the unique pointer
    // directly leads to incorrect results and crashes.
    cls.def_static("decode",
                   [](py::bytes bytes) { return condecode(bytes).release(); },
                   "bytes"_a);

    cls.def("__repr__", [](ConvexPolygon const &self) {
        return py::str("ConvexPolygon({!r})").format(self.getVertices());
    });
    cls.def(py::pickle(
            [](const ConvexPolygon &self) { return python::encode(self); },
            [](py::bytes bytes) { return condecode(bytes).release(); }));
}

void defineCurve(py::module &mod) {
    mod.def("log2", (uint8_t(*)(uint64_t)) & sphgeom::log2);
    mod.def("mortonIndex", (uint64_t(*)(uint32_t, uint32_t)) & sphgeom::mortonIndex,
            "x"_a, "y"_a);
    mod.def("mortonIndexInverse",
            (std::tuple<uint32_t, uint32_t>(*)(uint64_t)) & sphgeom::mortonIndexInverse,
            "z"_a);
    mod.def("mortonToHilbert", &sphgeom::mortonToHilbert, "z"_a, "m"_a);
    mod.def("hilbertToMorton", &sphgeom::hilbertToMorton, "h"_a, "m"_a);
    mod.def("hilbertIndex",
            (uint64_t(*)(uint32_t, uint32_t, int)) &sphgeom::hilbertIndex, "x"_a, "y"_a,
            "m"_a);
    mod.def("hilbertIndexInverse",
            (std::tuple<uint32_t, uint32_t>(*)(uint64_t, int))
                    &sphgeom::hilbertIndexInverse,
            "h"_a, "m"_a);
}

template <>
void defineClass(py::class_<Ellipse, std::unique_ptr<Ellipse>, Region> &cls) {
    cls.attr("TYPE_CODE") = py::int_(Ellipse::TYPE_CODE);

    cls.def_static("empty", &Ellipse::empty);
    cls.def_static("full", &Ellipse::full);

    cls.def(py::init<>());
    cls.def(py::init<Circle const &>(), "circle"_a);
    cls.def(py::init<UnitVector3d const &, Angle>(), "center"_a,
            "angle"_a = Angle(0.0));
    cls.def(py::init<UnitVector3d const &, UnitVector3d const &, Angle>(),
            "focus1"_a, "focus2"_a, "alpha"_a);
    cls.def(py::init<UnitVector3d const &, Angle, Angle, Angle>(), "center"_a,
            "alpha"_a, "beta"_a, "orientation"_a);
    cls.def(py::init<Ellipse const &>(), "ellipse"_a);

    cls.def("__eq__", &Ellipse::operator==, py::is_operator());
    cls.def("__ne__", &Ellipse::operator!=, py::is_operator());

    cls.def("isEmpty", &Ellipse::isEmpty);
    cls.def("isFull", &Ellipse::isFull);
    cls.def("isGreatCircle", &Ellipse::isGreatCircle);
    cls.def("isCircle", &Ellipse::isCircle);
    cls.def("getTransformMatrix", &Ellipse::getTransformMatrix);
    cls.def("getCenter", &Ellipse::getCenter);
    cls.def("getF1", &Ellipse::getF1);
    cls.def("getF2", &Ellipse::getF2);
    cls.def("getAlpha", &Ellipse::getAlpha);
    cls.def("getBeta", &Ellipse::getBeta);
    cls.def("getGamma", &Ellipse::getGamma);
    cls.def("complement", &Ellipse::complement);
    cls.def("complemented", &Ellipse::complemented);

    // Note that the Region interface has already been wrapped.

    // The lambda is necessary for now; returning the unique pointer
    // directly leads to incorrect results and crashes.
    cls.def_static("decode",
                   [](py::bytes bytes) { return decode(bytes).release(); },
                   "bytes"_a);

    cls.def("__str__", [](Ellipse const &self) {
        return py::str("Ellipse({!s}, {!s}, {!s})")
                .format(self.getF1(), self.getF2(), self.getAlpha());
    });
    cls.def("__repr__", [](Ellipse const &self) {
        return py::str("Ellipse({!r}, {!r}, {!r})")
                .format(self.getF1(), self.getF2(), self.getAlpha());
    });
    cls.def(py::pickle(
            [](const Ellipse &self) { return python::encode(self); },
            [](py::bytes bytes) { return decode(bytes).release(); }));
}

template <>
void defineClass(py::class_<LonLat, std::shared_ptr<LonLat>> &cls) {
    cls.def_static("fromDegrees", &LonLat::fromDegrees);
    cls.def_static("fromRadians", &LonLat::fromRadians);
    cls.def_static("latitudeOf", &LonLat::latitudeOf);
    cls.def_static("longitudeOf", &LonLat::longitudeOf);

    cls.def(py::init<>());
    cls.def(py::init<LonLat const &>());
    cls.def(py::init<NormalizedAngle, Angle>(), "lon"_a, "lat"_a);
    cls.def(py::init<Vector3d const &>(), "vector"_a);

    cls.def("__eq__", &LonLat::operator==, py::is_operator());
    cls.def("__nq__", &LonLat::operator!=, py::is_operator());

    cls.def("getLon", &LonLat::getLon);
    cls.def("getLat", &LonLat::getLat);

    cls.def("__len__", [](LonLat const &self) { return py::int_(2); });
    cls.def("__getitem__", [](LonLat const &self, py::object key) {
        auto t = py::make_tuple(self.getLon(), self.getLat());
        return t.attr("__getitem__")(key);
    });
    cls.def("__iter__", [](LonLat const &self) {
        auto t = py::make_tuple(self.getLon(), self.getLat());
        return t.attr("__iter__")();
    });

    cls.def("__str__", [](LonLat const &self) {
        return py::str("[{!s}, {!s}]")
                .format(self.getLon().asRadians(), self.getLat().asRadians());
    });
    cls.def("__repr__", [](LonLat const &self) {
        return py::str("LonLat.fromRadians({!r}, {!r})")
                .format(self.getLon().asRadians(), self.getLat().asRadians());
    });
    cls.def("__reduce__", [cls](LonLat const &self) {
        return py::make_tuple(cls,
                              py::make_tuple(self.getLon(), self.getLat()));
    });
}

template <>
void defineClass(py::class_<Interval1d, std::shared_ptr<Interval1d>> &cls) {
    python::defineInterval<decltype(cls), Interval1d, double>(cls);

    cls.def_static("empty", &Interval1d::empty);
    cls.def_static("full", &Interval1d::full);

    cls.def(py::init<>());
    cls.def(py::init<double>(), "x"_a);
    cls.def(py::init<double, double>(), "x"_a, "y"_a);
    cls.def(py::init<Interval1d const &>(), "interval"_a);

    cls.def("isFull", &Interval1d::isFull);

    cls.def("__str__", [](Interval1d const &self) {
        return py::str("[{!s}, {!s}]").format(self.getA(), self.getB());
    });
    cls.def("__repr__", [](Interval1d const &self) {
        return py::str("Interval1d({!r}, {!r})")
                .format(self.getA(), self.getB());
    });
}

template <>
void defineClass(py::class_<HtmPixelization, Pixelization> &cls) {
    cls.attr("MAX_LEVEL") = py::int_(HtmPixelization::MAX_LEVEL);

    cls.def_static("level", &HtmPixelization::level, "i"_a);
    cls.def_static("triangle", &HtmPixelization::triangle, "i"_a);
    cls.def_static("asString", &HtmPixelization::asString, "i"_a);

    cls.def(py::init<int>(), "level"_a);
    cls.def(py::init<HtmPixelization const &>(), "htmPixelization"_a);

    cls.def("getLevel", &HtmPixelization::getLevel);

    cls.def("__eq__",
            [](HtmPixelization const &self, HtmPixelization const &other) {
                return self.getLevel() == other.getLevel();
            });
    cls.def("__ne__",
            [](HtmPixelization const &self, HtmPixelization const &other) {
                return self.getLevel() != other.getLevel();
            });
    cls.def("__repr__", [](HtmPixelization const &self) {
        return py::str("HtmPixelization({!s})").format(self.getLevel());
    });
    cls.def("__reduce__", [cls](HtmPixelization const &self) {
        return py::make_tuple(cls, py::make_tuple(self.getLevel()));
    });
}

template <>
void defineClass(py::class_<Matrix3d, std::shared_ptr<Matrix3d>> &cls) {
    cls.def(py::init<>());
    cls.def(py::init<double, double, double, double, double, double, double,
                     double, double>(),
            "m00"_a, "m01"_a, "m02"_a, "m10"_a, "m11"_a, "m12"_a, "m20"_a,
            "m21"_a, "m22"_a);
    cls.def(py::init<Vector3d const &>(), "diagonal"_a);
    cls.def(py::init<double>(), "scale"_a);
    cls.def(py::init<Matrix3d const &>(), "matrix"_a);

    cls.def("__eq__", &Matrix3d::operator==, py::is_operator());
    cls.def("__ne__", &Matrix3d::operator!=, py::is_operator());

    // Add bounds checking to getRow and getColumn
    cls.def("getRow", &getRow, "row"_a);
    cls.def("getColumn",
            [](Matrix3d const &self, py::int_ col) {
                return self.getColumn(
                        static_cast<int>(python::convertIndex(3, col)));
            },
            "col"_a);

    cls.def("__len__", [](Matrix3d const &self) { return py::int_(3); });
    cls.def("__getitem__", &getRow, py::is_operator());
    cls.def("__getitem__",
            [](Matrix3d const &self, py::tuple t) {
                if (t.size() > 2) {
                    throw py::index_error("Too many indexes for Matrix3d");
                } else if (t.size() == 0) {
                    return py::cast(self);
                } else if (t.size() == 1) {
                    return py::cast(getRow(self, t[0].cast<py::int_>()));
                }
                return py::cast(
                        self(python::convertIndex(3, t[0].cast<py::int_>()),
                             python::convertIndex(3, t[1].cast<py::int_>())));
            },
            py::is_operator());

    cls.def("inner", &Matrix3d::inner, "matrix"_a);
    cls.def("getSquaredNorm", &Matrix3d::getSquaredNorm);
    cls.def("getNorm", &Matrix3d::getNorm);

    cls.def("__mul__",
            (Vector3d(Matrix3d::*)(Vector3d const &) const) &
                    Matrix3d::operator*,
            "vector"_a, py::is_operator());
    cls.def("__mul__",
            (Matrix3d(Matrix3d::*)(Matrix3d const &) const) &
                    Matrix3d::operator*,
            "matrix"_a, py::is_operator());
    cls.def("__add__", &Matrix3d::operator+, py::is_operator());
    cls.def("__sub__", &Matrix3d::operator-, py::is_operator());

    cls.def("cwiseProduct", &Matrix3d::cwiseProduct);
    cls.def("transpose", &Matrix3d::transpose);
    cls.def("inverse", &Matrix3d::inverse);

    cls.def("__str__", [](Matrix3d const &self) {
        return py::str("[[{!s}, {!s}, {!s}],\n"
                       " [{!s}, {!s}, {!s}],\n"
                       " [{!s}, {!s}, {!s}]]")
                .format(self(0, 0), self(0, 1), self(0, 2), self(1, 0),
                        self(1, 1), self(1, 2), self(2, 0), self(2, 1),
                        self(2, 2));
    });
    cls.def("__repr__", [](Matrix3d const &self) {
        return py::str("Matrix3d({!r}, {!r}, {!r},\n"
                       "         {!r}, {!r}, {!r},\n"
                       "         {!r}, {!r}, {!r})")
                .format(self(0, 0), self(0, 1), self(0, 2), self(1, 0),
                        self(1, 1), self(1, 2), self(2, 0), self(2, 1),
                        self(2, 2));
    });
    cls.def("__reduce__", [cls](Matrix3d const &self) {
        auto args = py::make_tuple(self(0, 0), self(0, 1), self(0, 2),
                                   self(1, 0), self(1, 1), self(1, 2),
                                   self(2, 0), self(2, 1), self(2, 2));
        return py::make_tuple(cls, args);
    });
}

template <>
void defineClass(py::class_<Pixelization> &cls) {
    cls.def("universe", &Pixelization::universe);
    cls.def("pixel", &Pixelization::pixel, "i"_a);
    cls.def("index", &Pixelization::index, "i"_a);
    cls.def("toString", &Pixelization::toString, "i"_a);
    cls.def("envelope", &Pixelization::envelope, "region"_a, "maxRanges"_a = 0);
    cls.def("interior", &Pixelization::interior, "region"_a, "maxRanges"_a = 0);
}

void defineOrientation(py::module &mod) {
    mod.def("orientationExact", &sphgeom::orientationExact, "a"_a, "b"_a, "c"_a);
    mod.def("orientation", &sphgeom::orientation, "a"_a, "b"_a, "c"_a);
    mod.def("orientationX", &sphgeom::orientationX, "b"_a, "c"_a);
    mod.def("orientationY", &sphgeom::orientationY, "b"_a, "c"_a);
    mod.def("orientationZ", &sphgeom::orientationZ, "b"_a, "c"_a);
}

template <>
void defineClass(py::class_<Mq3cPixelization, Pixelization> &cls) {
    cls.attr("MAX_LEVEL") = py::int_(Mq3cPixelization::MAX_LEVEL);

    cls.def_static("level", &Mq3cPixelization::level);
    cls.def_static("quad", &Mq3cPixelization::quad);
    cls.def_static("neighborhood", &Mq3cPixelization::neighborhood);
    cls.def_static("asString", &Mq3cPixelization::asString);

    cls.def(py::init<int>(), "level"_a);
    cls.def(py::init<Mq3cPixelization const &>(), "mq3cPixelization"_a);

    cls.def("getLevel", &Mq3cPixelization::getLevel);

    cls.def("__eq__",
            [](Mq3cPixelization const &self, Mq3cPixelization const &other) {
                return self.getLevel() == other.getLevel();
            });
    cls.def("__ne__",
            [](Mq3cPixelization const &self, Mq3cPixelization const &other) {
                return self.getLevel() != other.getLevel();
            });
    cls.def("__repr__", [](Mq3cPixelization const &self) {
        return py::str("Mq3cPixelization({!s})").format(self.getLevel());
    });
    cls.def("__reduce__", [cls](Mq3cPixelization const &self) {
        return py::make_tuple(cls, py::make_tuple(self.getLevel()));
    });
}

template <>
void defineClass(py::class_<NormalizedAngleInterval,
                            std::shared_ptr<NormalizedAngleInterval>> &cls) {
    python::defineInterval<decltype(cls), NormalizedAngleInterval,
                           NormalizedAngle>(cls);

    cls.def_static("fromDegrees", &NormalizedAngleInterval::fromDegrees, "x"_a,
                   "y"_a);
    cls.def_static("fromRadians", &NormalizedAngleInterval::fromRadians, "x"_a,
                   "y"_a);
    cls.def_static("empty", &NormalizedAngleInterval::empty);
    cls.def_static("full", &NormalizedAngleInterval::full);

    cls.def(py::init<>());
    cls.def(py::init<Angle>(), "x"_a);
    cls.def(py::init<NormalizedAngle>(), "x"_a);
    cls.def(py::init<Angle, Angle>(), "x"_a, "y"_a);
    cls.def(py::init<NormalizedAngle, NormalizedAngle>(), "x"_a, "y"_a);
    cls.def(py::init<NormalizedAngleInterval const &>(), "angleInterval"_a);

    cls.def("isEmpty", &NormalizedAngleInterval::isEmpty);
    cls.def("isFull", &NormalizedAngleInterval::isFull);
    cls.def("wraps", &NormalizedAngleInterval::wraps);

    cls.def("__str__", [](NormalizedAngleInterval const &self) {
        return py::str("[{!s}, {!s}]")
                .format(self.getA().asRadians(), self.getB().asRadians());
    });
    cls.def("__repr__", [](NormalizedAngleInterval const &self) {
        return py::str("NormalizedAngleInterval.fromRadians({!r},"
                       " {!r})")
                .format(self.getA().asRadians(), self.getB().asRadians());
    });
}

template <>
void defineClass(py::class_<NormalizedAngle> &cls) {
    // Provide the equivalent of the NormalizedAngle to Angle C++ cast
    // operator in Python
    py::implicitly_convertible<NormalizedAngle, Angle>();

    cls.def_static("nan", &NormalizedAngle::nan);
    cls.def_static("fromDegrees", &NormalizedAngle::fromDegrees);
    cls.def_static("fromRadians", &NormalizedAngle::fromRadians);
    cls.def_static("between", &NormalizedAngle::between, "a"_a, "b"_a);
    cls.def_static("center", &NormalizedAngle::center, "a"_a, "b"_a);

    cls.def(py::init<>());
    cls.def(py::init<NormalizedAngle const &>());
    cls.def(py::init<Angle const &>());
    cls.def(py::init<double>(), "radians"_a);
    cls.def(py::init<LonLat const &, LonLat const &>(), "a"_a, "b"_a);
    cls.def(py::init<Vector3d const &, Vector3d const &>(), "a"_a, "b"_a);

    cls.def("__eq__", &NormalizedAngle::operator==, py::is_operator());
    cls.def("__ne__", &NormalizedAngle::operator!=, py::is_operator());
    cls.def("__lt__", &NormalizedAngle::operator<, py::is_operator());
    cls.def("__gt__", &NormalizedAngle::operator>, py::is_operator());
    cls.def("__le__", &NormalizedAngle::operator<=, py::is_operator());
    cls.def("__ge__", &NormalizedAngle::operator>=, py::is_operator());

    cls.def("__neg__",
            (Angle(NormalizedAngle::*)() const) & NormalizedAngle::operator-);
    cls.def("__add__", &NormalizedAngle::operator+, py::is_operator());
    cls.def("__sub__",
            (Angle(NormalizedAngle::*)(Angle const &) const) &
                    NormalizedAngle::operator-,
            py::is_operator());
    cls.def("__mul__", &NormalizedAngle::operator*, py::is_operator());
    cls.def("__rmul__", &NormalizedAngle::operator*, py::is_operator());
    cls.def("__truediv__",
            (Angle(NormalizedAngle::*)(double) const) &
                    NormalizedAngle::operator/,
            py::is_operator());
    cls.def("__truediv__",
            (double (NormalizedAngle::*)(Angle const &) const) &
                    NormalizedAngle::operator/,
            py::is_operator());

    cls.def("asDegrees", &NormalizedAngle::asDegrees);
    cls.def("asRadians", &NormalizedAngle::asRadians);
    cls.def("isNan", &NormalizedAngle::isNan);
    cls.def("getAngleTo", &NormalizedAngle::getAngleTo);

    cls.def("__str__", [](NormalizedAngle const &self) {
        return py::str("{!s}").format(self.asRadians());
    });
    cls.def("__repr__", [](NormalizedAngle const &self) {
        return py::str("NormalizedAngle({!r})").format(self.asRadians());
    });

    cls.def("__reduce__", [cls](NormalizedAngle const &self) {
        return py::make_tuple(cls, py::make_tuple(self.asRadians()));
    });
}

template <>
void defineClass(py::class_<sphgeom::Vector3d, std::shared_ptr<sphgeom::Vector3d>> &cls) {
    cls.def(py::init<>());
    cls.def(py::init<double, double, double>(), "x"_a, "y"_a, "z"_a);
    cls.def(py::init<sphgeom::Vector3d const &>(), "vector"_a);
    // Construct a Vector3d from a UnitVector3d, enabling implicit
    // conversion from UnitVector3d to Vector3d in python via
    // py::implicitly_convertible
    cls.def(py::init([](sphgeom::UnitVector3d const &u) {
        return new sphgeom::Vector3d(u.x(), u.y(), u.z());
    }));

    cls.def("__eq__", &sphgeom::Vector3d::operator==, py::is_operator());
    cls.def("__ne__", &sphgeom::Vector3d::operator!=, py::is_operator());
    cls.def("__neg__", (sphgeom::Vector3d(sphgeom::Vector3d::*)() const) & sphgeom::Vector3d::operator-);
    cls.def("__add__", &sphgeom::Vector3d::operator+, py::is_operator());
    cls.def("__sub__",
            (sphgeom::Vector3d(sphgeom::Vector3d::*)(sphgeom::Vector3d const &) const) &
                    sphgeom::Vector3d::operator-,
            py::is_operator());
    cls.def("__mul__", &sphgeom::Vector3d::operator*, py::is_operator());
    cls.def("__truediv__", &sphgeom::Vector3d::operator/, py::is_operator());

    cls.def("__iadd__", &sphgeom::Vector3d::operator+=);
    cls.def("__isub__", &sphgeom::Vector3d::operator-=);
    cls.def("__imul__", &sphgeom::Vector3d::operator*=);
    cls.def("__itruediv__", &sphgeom::Vector3d::operator/=);

    cls.def("x", &sphgeom::Vector3d::x);
    cls.def("y", &sphgeom::Vector3d::y);
    cls.def("z", &sphgeom::Vector3d::z);
    cls.def("dot", &sphgeom::Vector3d::dot);
    cls.def("getSquaredNorm", &sphgeom::Vector3d::getSquaredNorm);
    cls.def("getNorm", &sphgeom::Vector3d::getNorm);
    cls.def("isZero", &sphgeom::Vector3d::isZero);
    cls.def("normalize", &sphgeom::Vector3d::normalize);
    cls.def("isNormalized", &sphgeom::Vector3d::isNormalized);
    cls.def("cross", &sphgeom::Vector3d::cross);
    cls.def("cwiseProduct", &sphgeom::Vector3d::cwiseProduct);
    cls.def("rotatedAround", &sphgeom::Vector3d::rotatedAround, "axis"_a, "angle"_a);

    cls.def("__len__", [](sphgeom::Vector3d const &self) { return py::int_(3); });
    cls.def("__getitem__", [](sphgeom::Vector3d const &self, py::int_ i) {
        return self(python::convertIndex(3, i));
    });

    cls.def("__str__", [](sphgeom::Vector3d const &self) {
        return py::str("[{!s}, {!s}, {!s}]")
                .format(self.x(), self.y(), self.z());
    });
    cls.def("__repr__", [](sphgeom::Vector3d const &self) {
        return py::str("Vector3d({!r}, {!r}, {!r})")
                .format(self.x(), self.y(), self.z());
    });

    cls.def("__reduce__", [cls](sphgeom::Vector3d const &self) {
        return py::make_tuple(cls,
                              py::make_tuple(self.x(), self.y(), self.z()));
    });
}

void defineUtils(py::module &mod) {
    mod.def("getMinSquaredChordLength", &sphgeom::getMinSquaredChordLength, "v"_a, "a"_a,
            "b"_a, "n"_a);
    mod.def("getMaxSquaredChordLength", &sphgeom::getMaxSquaredChordLength, "v"_a, "a"_a,
            "b"_a, "n"_a);
    mod.def("getMinAngleToCircle", &sphgeom::getMinAngleToCircle, "x"_a, "c"_a);
    mod.def("getMaxAngleToCircle", &sphgeom::getMaxAngleToCircle, "x"_a, "c"_a);
    mod.def("getWeightedCentroid", &sphgeom::getWeightedCentroid, "vector0"_a,
            "vector1"_a, "vector2"_a);
}

template <>
void defineClass(py::class_<Region, std::unique_ptr<Region>> &cls) {
    cls.def("clone", [](Region const &self) { return self.clone().release(); });
    cls.def("getBoundingBox", &Region::getBoundingBox);
    cls.def("getBoundingBox3d", &Region::getBoundingBox3d);
    cls.def("getBoundingCircle", &Region::getBoundingCircle);
    cls.def("contains", &Region::contains, "unitVector"_a);
    cls.def("__contains__", &Region::contains, "unitVector"_a,
            py::is_operator());
    // The per-subclass relate() overloads are used to implement
    // double-dispatch in C++, and are not needed in Python.
    cls.def("relate",
            (Relationship(Region::*)(Region const &) const) & Region::relate,
            "region"_a);
    cls.def("encode", &python::encode);

    cls.def_static(
            "decode",
            [](py::bytes bytes) {
                uint8_t const *buffer = reinterpret_cast<uint8_t const *>(
                        PYBIND11_BYTES_AS_STRING(bytes.ptr()));
                size_t n =
                        static_cast<size_t>(PYBIND11_BYTES_SIZE(bytes.ptr()));
                return Region::decode(buffer, n).release();
            },
            "bytes"_a);
}

template <>
void defineClass(py::class_<UnitVector3d, std::shared_ptr<UnitVector3d>> &cls) {
    // Provide the equivalent of the UnitVector3d to Vector3d C++ cast
    // operator in Python
    py::implicitly_convertible<UnitVector3d, Vector3d>();

    cls.def_static(
            "orthogonalTo",
            (UnitVector3d(*)(Vector3d const &)) & UnitVector3d::orthogonalTo,
            "vector"_a);
    cls.def_static("orthogonalTo",
                   (UnitVector3d(*)(Vector3d const &, Vector3d const &)) &
                           UnitVector3d::orthogonalTo,
                   "vector1"_a, "vector2"_a);
    cls.def_static("orthogonalTo",
                   (UnitVector3d(*)(NormalizedAngle const &)) &
                           UnitVector3d::orthogonalTo,
                   "meridian"_a);
    cls.def_static("northFrom", &UnitVector3d::northFrom, "vector"_a);
    cls.def_static("X", &UnitVector3d::X);
    cls.def_static("Y", &UnitVector3d::Y);
    cls.def_static("Z", &UnitVector3d::Z);
    // The fromNormalized static factory functions are not exposed to
    // Python, as they are easy to misuse and intended only for performance
    // critical code (i.e. not Python).

    cls.def(py::init<>());
    cls.def(py::init<UnitVector3d const &>(), "unitVector"_a);
    cls.def(py::init<Vector3d const &>(), "vector"_a);
    cls.def(py::init<double, double, double>(), "x"_a, "y"_a, "z"_a);
    cls.def(py::init<LonLat const &>(), "lonLat"_a);
    cls.def(py::init<Angle, Angle>(), "lon"_a, "lat"_a);

    cls.def("__eq__", &UnitVector3d::operator==, py::is_operator());
    cls.def("__ne__", &UnitVector3d::operator!=, py::is_operator());
    cls.def("__neg__",
            (UnitVector3d(UnitVector3d::*)() const) & UnitVector3d::operator-);
    cls.def("__add__", &UnitVector3d::operator+, py::is_operator());
    cls.def("__sub__",
            (Vector3d(UnitVector3d::*)(Vector3d const &) const) &
                    UnitVector3d::operator-,
            py::is_operator());
    cls.def("__mul__", &UnitVector3d::operator*, py::is_operator());
    cls.def("__truediv__", &UnitVector3d::operator/, py::is_operator());

    cls.def("x", &UnitVector3d::x);
    cls.def("y", &UnitVector3d::y);
    cls.def("z", &UnitVector3d::z);
    cls.def("x", &UnitVector3d::dot);
    cls.def("dot", &UnitVector3d::dot);
    cls.def("cross", &UnitVector3d::cross);
    cls.def("robustCross", &UnitVector3d::robustCross);
    cls.def("cwiseProduct", &UnitVector3d::cwiseProduct);
    cls.def("rotatedAround", &UnitVector3d::rotatedAround, "axis"_a, "angle"_a);

    cls.def("__len__", [](UnitVector3d const &self) { return py::int_(3); });
    cls.def("__getitem__", [](UnitVector3d const &self, py::int_ i) {
        return self(python::convertIndex(3, i));
    });

    cls.def("__str__", [](UnitVector3d const &self) {
        return py::str("[{!s}, {!s}, {!s}]")
                .format(self.x(), self.y(), self.z());
    });
    cls.def("__repr__", [](UnitVector3d const &self) {
        return py::str("UnitVector3d({!r}, {!r}, {!r})")
                .format(self.x(), self.y(), self.z());
    });

    // Do not implement __reduce__ for pickling. Why? Given:
    //
    //    u = UnitVector3d(x, y, z)
    //    v = UnitVector3d(u.x(), u.y(), u.z())
    //
    // u may not be identical to v, since the constructor normalizes its input
    // components. Furthermore, UnitVector3d::fromNormalized is not visible to
    // Python, and even if it were, pybind11 is currently incapable of returning
    // a picklable reference to it.
    cls.def(py::pickle([](UnitVector3d const &self) { return py::make_tuple(self.x(), self.y(), self.z()); },
                       [](py::tuple t) {
                           if (t.size() != 3) {
                               throw std::runtime_error("Tuple size = " + std::to_string(t.size()) +
                                                        "; must be 3 for a UnitVector3d");
                           }
                           return new UnitVector3d(UnitVector3d::fromNormalized(
                                   t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>()));
                       }));
}

template <>
void defineClass(py::class_<RangeSet, std::shared_ptr<RangeSet>> &cls) {
    cls.def(py::init<>());
    cls.def(py::init<uint64_t>(), "integer"_a);
    cls.def(py::init([](uint64_t a, uint64_t b) {
                return new RangeSet(a, b);
            }),
            "first"_a, "last"_a);
    cls.def(py::init<RangeSet const &>(), "rangeSet"_a);
    cls.def(py::init(
            [](py::iterable iterable) {
                return new RangeSet(makeRangeSet(iterable));
            }),
            "iterable"_a);
    cls.def("__eq__", &RangeSet::operator==, py::is_operator());
    cls.def("__ne__", &RangeSet::operator!=, py::is_operator());

    cls.def("insert", (void (RangeSet::*)(uint64_t)) & RangeSet::insert,
            "integer"_a);
    cls.def("insert",
            (void (RangeSet::*)(uint64_t, uint64_t)) & RangeSet::insert,
            "first"_a, "last"_a);
    cls.def("erase", (void (RangeSet::*)(uint64_t)) & RangeSet::erase,
            "integer"_a);
    cls.def("erase", (void (RangeSet::*)(uint64_t, uint64_t)) & RangeSet::erase,
            "first"_a, "last"_a);

    cls.def("complement", &RangeSet::complement);
    cls.def("complemented", &RangeSet::complemented);
    cls.def("intersection", &RangeSet::intersection, "rangeSet"_a);
    // In C++, the set union function is named join because union is a keyword.
    // Python does not suffer from the same restriction.
    cls.def("union", &RangeSet::join, "rangeSet"_a);
    cls.def("difference", &RangeSet::difference, "rangeSet"_a);
    cls.def("symmetricDifference", &RangeSet::symmetricDifference,
            "rangeSet"_a);
    cls.def("__invert__", &RangeSet::operator~, py::is_operator());
    cls.def("__and__", &RangeSet::operator&, py::is_operator());
    cls.def("__or__", &RangeSet::operator|, py::is_operator());
    cls.def("__sub__", &RangeSet::operator-, py::is_operator());
    cls.def("__xor__", &RangeSet::operator^, py::is_operator());
    cls.def("__iand__", &RangeSet::operator&=);
    cls.def("__ior__", &RangeSet::operator|=);
    cls.def("__isub__", &RangeSet::operator-=);
    cls.def("__ixor__", &RangeSet::operator^=);

    cls.def("__len__", &RangeSet::size);
    cls.def("__getitem__", [](RangeSet const &self, py::int_ i) {
        auto j = python::convertIndex(static_cast<ptrdiff_t>(self.size()), i);
        return py::cast(self.begin()[j]);
    });

    cls.def("intersects",
            (bool (RangeSet::*)(uint64_t) const) & RangeSet::intersects,
            "integer"_a);
    cls.def("intersects",
            (bool (RangeSet::*)(uint64_t, uint64_t) const) &
                    RangeSet::intersects,
            "first"_a, "last"_a);
    cls.def("intersects",
            (bool (RangeSet::*)(RangeSet const &) const) & RangeSet::intersects,
            "rangeSet"_a);

    cls.def("contains",
            (bool (RangeSet::*)(uint64_t) const) & RangeSet::contains,
            "integer"_a);
    cls.def("contains",
            (bool (RangeSet::*)(uint64_t, uint64_t) const) & RangeSet::contains,
            "first"_a, "last"_a);
    cls.def("contains",
            (bool (RangeSet::*)(RangeSet const &) const) & RangeSet::contains,
            "rangeSet"_a);
    cls.def("__contains__",
            (bool (RangeSet::*)(uint64_t) const) & RangeSet::contains,
            "integer"_a, py::is_operator());
    cls.def("__contains__",
            (bool (RangeSet::*)(uint64_t, uint64_t) const) & RangeSet::contains,
            "first"_a, "last"_a, py::is_operator());
    cls.def("__contains__",
            (bool (RangeSet::*)(RangeSet const &) const) & RangeSet::contains,
            "rangeSet"_a, py::is_operator());

    cls.def("isWithin",
            (bool (RangeSet::*)(uint64_t) const) & RangeSet::isWithin,
            "integer"_a);
    cls.def("isWithin",
            (bool (RangeSet::*)(uint64_t, uint64_t) const) & RangeSet::isWithin,
            "first"_a, "last"_a);
    cls.def("isWithin",
            (bool (RangeSet::*)(RangeSet const &) const) & RangeSet::isWithin,
            "rangeSet"_a);

    cls.def("isDisjointFrom",
            (bool (RangeSet::*)(uint64_t) const) & RangeSet::isDisjointFrom,
            "integer"_a);
    cls.def("isDisjointFrom",
            (bool (RangeSet::*)(uint64_t, uint64_t) const) &
                    RangeSet::isDisjointFrom,
            "first"_a, "last"_a);
    cls.def("isDisjointFrom",
            (bool (RangeSet::*)(RangeSet const &) const) &
                    RangeSet::isDisjointFrom,
            "rangeSet"_a);

    cls.def("simplify", &RangeSet::simplify, "n"_a);
    cls.def("simplified", &RangeSet::simplified, "n"_a);
    cls.def("scale", &RangeSet::scale, "factor"_a);
    cls.def("scaled", &RangeSet::scaled, "factor"_a);
    cls.def("fill", &RangeSet::fill);
    cls.def("clear", &RangeSet::clear);
    cls.def("empty", &RangeSet::empty);
    cls.def("full", &RangeSet::full);
    cls.def("size", &RangeSet::size);
    cls.def("cardinality", &RangeSet::cardinality);
    // max_size() and swap() are omitted. The former is a C++ container
    // requirement, and the latter doesn't seem relevant to Python.
    cls.def("isValid", &RangeSet::cardinality);
    cls.def("ranges", &ranges);

    cls.def("__str__",
            [](RangeSet const &self) { return py::str(ranges(self)); });
    cls.def("__repr__", [](RangeSet const &self) {
        return py::str("RangeSet({!s})").format(ranges(self));
    });

    cls.def("__reduce__", [cls](RangeSet const &self) {
        return py::make_tuple(cls, py::make_tuple(ranges(self)));
    });
}

template <>
void defineClass(py::class_<Q3cPixelization, Pixelization> &cls) {
    cls.attr("MAX_LEVEL") = py::int_(Q3cPixelization::MAX_LEVEL);

    cls.def(py::init<int>(), "level"_a);
    cls.def(py::init<Q3cPixelization const &>(), "q3cPixelization"_a);

    cls.def("getLevel", &Q3cPixelization::getLevel);
    cls.def("quad", &Q3cPixelization::quad);
    cls.def("neighborhood", &Q3cPixelization::neighborhood);

    cls.def("__eq__",
            [](Q3cPixelization const &self, Q3cPixelization const &other) {
                return self.getLevel() == other.getLevel();
            });
    cls.def("__ne__",
            [](Q3cPixelization const &self, Q3cPixelization const &other) {
                return self.getLevel() != other.getLevel();
            });
    cls.def("__repr__", [](Q3cPixelization const &self) {
        return py::str("Q3cPixelization({!s})").format(self.getLevel());
    });
    cls.def("__reduce__", [cls](Q3cPixelization const &self) {
        return py::make_tuple(cls, py::make_tuple(self.getLevel()));
    });
}

void defineRelationship(py::module &mod) {
    mod.attr("DISJOINT") = py::cast(sphgeom::DISJOINT.to_ulong());
    mod.attr("INTERSECTS") = py::cast(sphgeom::INTERSECTS.to_ulong());
    mod.attr("CONTAINS") = py::cast(sphgeom::CONTAINS.to_ulong());
    mod.attr("WITHIN") = py::cast(sphgeom::WITHIN.to_ulong());

    mod.def("invert", &sphgeom::invert, "relationship"_a);
}


WRAP(Sphgeom) {
    py::class_<sphgeom::Angle> angle(mod, "Angle");
    py::class_<sphgeom::NormalizedAngle> normalizedAngle(mod, "NormalizedAngle");
    py::class_<sphgeom::LonLat, std::shared_ptr<sphgeom::LonLat>> lonLat(mod, "LonLat");
    py::class_<sphgeom::Vector3d, std::shared_ptr<sphgeom::Vector3d>> vector3d(mod, "Vector3d");
    py::class_<sphgeom::UnitVector3d, std::shared_ptr<sphgeom::UnitVector3d>> unitVector3d(
            mod, "UnitVector3d");
    py::class_<sphgeom::Matrix3d, std::shared_ptr<sphgeom::Matrix3d>> matrix3d(mod, "Matrix3d");

    py::class_<sphgeom::AngleInterval, std::shared_ptr<sphgeom::AngleInterval>> angleInterval(
            mod, "AngleInterval");
    py::class_<sphgeom::NormalizedAngleInterval,
               std::shared_ptr<sphgeom::NormalizedAngleInterval>>
            normalizedAngleInterval(mod, "NormalizedAngleInterval");
    py::class_<sphgeom::Interval1d, std::shared_ptr<sphgeom::Interval1d>> interval1d(
            mod, "Interval1d");

    py::class_<sphgeom::Box3d, std::shared_ptr<sphgeom::Box3d>> box3d(mod, "Box3d");

    py::class_<sphgeom::Region, std::unique_ptr<sphgeom::Region>> region(mod, "Region");
    py::class_<sphgeom::Box, std::unique_ptr<sphgeom::Box>, sphgeom::Region> box(mod, "Box");
    py::class_<sphgeom::Circle, std::unique_ptr<sphgeom::Circle>, sphgeom::Region> circle(mod, "Circle");
    py::class_<sphgeom::ConvexPolygon, std::unique_ptr<sphgeom::ConvexPolygon>, sphgeom::Region>
            convexPolygon(mod, "ConvexPolygon");
    py::class_<sphgeom::Ellipse, std::unique_ptr<sphgeom::Ellipse>, sphgeom::Region> ellipse(mod,
                                                                  "Ellipse");

    py::class_<sphgeom::RangeSet, std::shared_ptr<sphgeom::RangeSet>> rangeSet(mod, "RangeSet");

    py::class_<sphgeom::Pixelization> pixelization(mod, "Pixelization");
    py::class_<sphgeom::HtmPixelization, sphgeom::Pixelization> htmPixelization(
            mod, "HtmPixelization");
    py::class_<sphgeom::Mq3cPixelization, sphgeom::Pixelization> mq3cPixelization(
            mod, "Mq3cPixelization");
    py::class_<sphgeom::Q3cPixelization, sphgeom::Pixelization> q3cPixelization(
            mod, "Q3cPixelization");

    py::class_<sphgeom::Chunker, std::shared_ptr<sphgeom::Chunker>> chunker(mod, "Chunker");

    defineClass(angle);
    defineClass(normalizedAngle);
    defineClass(lonLat);
    defineClass(vector3d);
    defineClass(unitVector3d);
    defineClass(matrix3d);

    defineClass(angleInterval);
    defineClass(normalizedAngleInterval);
    defineClass(interval1d);

    defineClass(box3d);

    defineClass(region);
    defineClass(box);
    defineClass(circle);
    defineClass(convexPolygon);
    defineClass(ellipse);

    defineClass(rangeSet);

    defineClass(pixelization);
    defineClass(htmPixelization);
    defineClass(mq3cPixelization);
    defineClass(q3cPixelization);

    defineClass(chunker);

    // Define C++ functions.

    defineCurve(mod);
    defineOrientation(mod);
    defineRelationship(mod);
    defineUtils(mod);

}

}
}
