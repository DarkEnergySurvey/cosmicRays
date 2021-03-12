#include "pybind/afw_bind.h"
#include <vector>

#include "lsst/utils/python.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/afw/cameraGeom/Camera.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/cameraGeom/Amplifier.h"
#include "lsst/afw/cameraGeom/CameraSys.h"
#include "lsst/afw/cameraGeom/Orientation.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/afw/cameraGeom/CameraSys.h"
#include "lsst/afw/cameraGeom/TransformMap.h"
#include "lsst/afw/cameraGeom/DetectorCollection.h"
#include "lsst/geom.h"
#include "lsst/afw/typehandling/Storable.h"
#include "lsst/afw/cameraGeom/Detector.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace afw {
namespace cameraGeom {
using PyCamera = py::class_<Camera, DetectorCollection, std::shared_ptr<Camera>>;
using PyCameraBuilder = py::class_<Camera::Builder, DetectorCollectionBase<Detector::InCameraBuilder>,
                                   std::shared_ptr<Camera::Builder>>;
using PyAmplifier = py::class_<Amplifier, std::shared_ptr<Amplifier>>;
using PyAmplifierBuilder = py::class_<Amplifier::Builder, Amplifier, std::shared_ptr<Amplifier::Builder>>;
using PyTransformMap = py::class_<TransformMap, std::shared_ptr<TransformMap>>;
using PyTransformMapConnection = py::class_<TransformMap::Connection,
                                            std::shared_ptr<TransformMap::Connection>>;
using PyDetectorBase = py::class_<DetectorBase, std::shared_ptr<DetectorBase>>;
using PyDetector = py::class_<Detector, DetectorBase, std::shared_ptr<Detector>, typehandling::Storable>;
using PyDetectorBuilder = py::class_<Detector::Builder, DetectorBase, std::shared_ptr<Detector::Builder>>;
using PyDetectorPartialRebuilder = py::class_<Detector::PartialRebuilder, Detector::Builder,
                                              std::shared_ptr<Detector::PartialRebuilder>>;
using PyDetectorInCameraBuilder = py::class_<Detector::InCameraBuilder, Detector::Builder,
                                             std::shared_ptr<Detector::InCameraBuilder>>;

namespace {
// Bindings here are ordered to match the order of the declarations in
// Camera.h to the greatest extent possible; modifications to this file should
// attempt to preserve this.

void declareCameraBuilder(PyCamera & parent);

void declareCamera(py::module & mod) {
    PyCamera cls(mod, "Camera");
    declareCameraBuilder(cls);
    cls.def("rebuild", &Camera::rebuild);
    cls.def("getName", &Camera::getName);
    cls.def("getPupilFactoryName", &Camera::getPupilFactoryName);
    cls.def("findDetectors", &Camera::findDetectors, "point"_a, "cameraSys"_a);
    cls.def("findDetectorsList", &Camera::findDetectorsList, "pointList"_a, "cameraSys"_a);
    // transform methods are wrapped with lambdas that translate exceptions for backwards compatibility
    cls.def(
        "getTransform",
        [](Camera const & self, CameraSys const & fromSys, CameraSys const & toSys) {
            try {
                return self.getTransform(fromSys, toSys);
            } catch (pex::exceptions::NotFoundError & err) {
                PyErr_SetString(PyExc_KeyError, err.what());
                throw py::error_already_set();
            }
        },
        "fromSys"_a, "toSys"_a
    );
    cls.def("getTransformMap", &Camera::getTransformMap);
    cls.def(
        "transform",
        [](
            Camera const & self,
            lsst::geom::Point2D const & point,
            CameraSys const & fromSys,
            CameraSys const & toSys
        ) {
            try {
                return self.transform(point, fromSys, toSys);
            } catch (pex::exceptions::NotFoundError & err) {
                PyErr_SetString(PyExc_KeyError, err.what());
                throw py::error_already_set();
            }
        },
        "point"_a, "fromSys"_a, "toSys"_a
    );
    cls.def(
        "transform",
        [](
            Camera const & self,
            std::vector<lsst::geom::Point2D> const & points,
            CameraSys const & fromSys,
            CameraSys const & toSys
        ) {
            try {
                return self.transform(points, fromSys, toSys);
            } catch (pex::exceptions::NotFoundError & err) {
                PyErr_SetString(PyExc_KeyError, err.what());
                throw py::error_already_set();
            }
        },
        "points"_a, "fromSys"_a, "toSys"_a
    );
    table::io::python::addPersistableMethods(cls);
}

void declareCameraBuilder(PyCamera & parent) {
    PyCameraBuilder cls(parent, "Builder");
    cls.def(py::init<std::string const &>(), "name"_a);
    cls.def(py::init<Camera const &>(), "camera"_a);
    cls.def("finish", &Camera::Builder::finish);
    cls.def("getName", &Camera::Builder::getName);
    cls.def("setName", &Camera::Builder::setName);
    cls.def("getPupilFactoryName", &Camera::Builder::getPupilFactoryName);
    cls.def("setPupilFactoryName", &Camera::Builder::setPupilFactoryName);
    cls.def("setPupilFactoryClass",
            [](Camera::Builder & self, py::object pupilFactoryClass) {
                std::string pupilFactoryName = "lsst.afw.cameraGeom.pupil.PupilFactory";
                if (!pupilFactoryClass.is(py::none())) {
                    pupilFactoryName = py::str("{}.{}").format(
                        pupilFactoryClass.attr("__module__"),
                        pupilFactoryClass.attr("__name__")
                    );
                }
                self.setPupilFactoryName(pupilFactoryName);
            });
    cls.def("setTransformFromFocalPlaneTo", &Camera::Builder::setTransformFromFocalPlaneTo,
            "toSys"_a, "transform"_a);
    cls.def("discardTransformFromFocalPlaneTo",&Camera::Builder::discardTransformFromFocalPlaneTo);
    cls.def("add", &Camera::Builder::add);
    cls.def("__delitem__", py::overload_cast<int>(&Camera::Builder::remove));
    cls.def("__delitem__", py::overload_cast<std::string const &>(&Camera::Builder::remove));
}

WRAP(Camera) {

    declareCamera(mod);
}

PyAmplifier declarePyAmplifier(py::module & mod) {
    py::enum_<ReadoutCorner>(mod, "ReadoutCorner")
        .value("LL", ReadoutCorner::LL)
        .value("LR", ReadoutCorner::LR)
        .value("UR", ReadoutCorner::UR)
        .value("UL", ReadoutCorner::UL);
    py::enum_<AssemblyState>(mod, "AssemblyState")
        .value("RAW", AssemblyState::RAW)
        .value("SCIENCE", AssemblyState::SCIENCE);
    PyAmplifier cls(mod, "Amplifier");
    cls.def_static("getRecordSchema", &Amplifier::getRecordSchema);
    cls.def("toRecord", &Amplifier::toRecord);
    cls.def("rebuild", &Amplifier::rebuild);
    cls.def("getName", &Amplifier::getName);
    cls.def("getBBox", &Amplifier::getBBox);
    cls.def("getGain", &Amplifier::getGain);
    cls.def("getReadNoise", &Amplifier::getReadNoise);
    cls.def("getSaturation", &Amplifier::getSaturation);
    cls.def("getSuspectLevel", &Amplifier::getSuspectLevel);
    cls.def("getReadoutCorner", &Amplifier::getReadoutCorner);
    cls.def("getLinearityCoeffs", &Amplifier::getLinearityCoeffs);
    cls.def("getLinearityType", &Amplifier::getLinearityType);
    cls.def("getLinearityThreshold", &Amplifier::getLinearityThreshold);
    cls.def("getLinearityMaximum", &Amplifier::getLinearityMaximum);
    cls.def("getLinearityUnits", &Amplifier::getLinearityUnits);
    cls.def("getRawBBox", &Amplifier::getRawBBox);
    cls.def("getRawDataBBox", &Amplifier::getRawDataBBox);
    cls.def("getRawFlipX", &Amplifier::getRawFlipX);
    cls.def("getRawFlipY", &Amplifier::getRawFlipY);
    cls.def("getRawXYOffset", &Amplifier::getRawXYOffset);
    cls.def("getRawHorizontalOverscanBBox", &Amplifier::getRawHorizontalOverscanBBox);
    cls.def("getRawVerticalOverscanBBox", &Amplifier::getRawVerticalOverscanBBox);
    cls.def("getRawPrescanBBox", &Amplifier::getRawPrescanBBox);
    cls.def("getRawSerialOverscanBBox", &Amplifier::getRawSerialOverscanBBox);
    cls.def("getRawParallelOverscanBBox", &Amplifier::getRawParallelOverscanBBox);
    cls.def("getRawSerialPrescanBBox", &Amplifier::getRawSerialPrescanBBox);
    cls.def("getRawHorizontalPrescanBBox", &Amplifier::getRawHorizontalPrescanBBox);
    return cls;
}

void declarePyAmplifierBuilder(PyAmplifier & parent) {
    PyAmplifierBuilder cls(parent, "Builder");
    cls.def_static("fromRecord", &Amplifier::Builder::fromRecord);
    cls.def(py::init());
    cls.def("finish", &Amplifier::Builder::finish);
    cls.def("assign", [](Amplifier::Builder & self, Amplifier const & other) { self = other; });
    cls.def("setName", &Amplifier::Builder::setName, "name"_a);
    cls.def("setBBox", &Amplifier::Builder::setBBox, "bbox"_a);
    cls.def("setGain", &Amplifier::Builder::setGain, "gain"_a);
    cls.def("setReadNoise", &Amplifier::Builder::setReadNoise, "readNoise"_a);
    cls.def("setSaturation", &Amplifier::Builder::setSaturation, "saturation"_a);
    cls.def("setSuspectLevel", &Amplifier::Builder::setSuspectLevel, "suspectLevel"_a);
    cls.def("setReadoutCorner", &Amplifier::Builder::setReadoutCorner, "corner"_a);
    cls.def("setLinearityCoeffs", &Amplifier::Builder::setLinearityCoeffs, "coeffs"_a);
    // Backwards compatibility: accept std::vector (list in Python) in
    // addition to ndarray::Array (np.ndarray)
    cls.def("setLinearityCoeffs",
            [](Amplifier::Builder & self, std::vector<double> const & coeffs) {
                ndarray::Array<double, 1, 1> array = ndarray::allocate(coeffs.size());
                std::copy(coeffs.begin(), coeffs.end(), array.begin());
                self.setLinearityCoeffs(array);
            });
    cls.def("setLinearityType", &Amplifier::Builder::setLinearityType, "type"_a);
    cls.def("setLinearityThreshold", &Amplifier::Builder::setLinearityThreshold, "threshold"_a);
    cls.def("setLinearityMaximum", &Amplifier::Builder::setLinearityMaximum, "maximum"_a);
    cls.def("setLinearityUnits", &Amplifier::Builder::setLinearityUnits, "units"_a);
    cls.def("setRawBBox", &Amplifier::Builder::setRawBBox, "bbox"_a);
    cls.def("setRawDataBBox", &Amplifier::Builder::setRawDataBBox, "bbox"_a);
    cls.def("setRawFlipX", &Amplifier::Builder::setRawFlipX, "rawFlipX"_a);
    cls.def("setRawFlipY", &Amplifier::Builder::setRawFlipY, "rawFlipY"_a);
    cls.def("setRawXYOffset", &Amplifier::Builder::setRawXYOffset, "offset"_a);
    cls.def("setRawHorizontalOverscanBBox", &Amplifier::Builder::setRawHorizontalOverscanBBox, "bbox"_a);
    cls.def("setRawVerticalOverscanBBox", &Amplifier::Builder::setRawVerticalOverscanBBox, "bbox"_a);
    cls.def("setRawPrescanBBox", &Amplifier::Builder::setRawPrescanBBox, "bbox"_a);
    cls.def("setRawSerialOverscanBBox", &Amplifier::Builder::setRawSerialOverscanBBox, "bbox"_a);
    cls.def("setRawParallelOverscanBBox", &Amplifier::Builder::setRawParallelOverscanBBox, "bbox"_a);
    cls.def("setRawSerialPrescanBBox", &Amplifier::Builder::setRawSerialPrescanBBox, "bbox"_a);
    cls.def("setRawHorizontalPrescanBBox", &Amplifier::Builder::setRawHorizontalPrescanBBox, "bbox"_a);
}

WRAP(Amplifier) {
    auto cls = declarePyAmplifier(mod);
    declarePyAmplifierBuilder(cls);
}

/**
@internal Declare methods common to CameraSysPrefix and CameraSys

@tparam CppClass  C++ class; one of CameraSysPrefix or CameraSys
@tparam PyClass  pybind11 class corresponding to `CppClass`
*/
template <typename CppClass, typename PyClass>
void declareCommonSysMethods(PyClass &cls) {
    /* Operators */
    cls.def("__eq__", [](CppClass const &self, CppClass const &other) { return self == other; },
            py::is_operator());
    cls.def("__ne__", [](CppClass const &self, CppClass const &other) { return self != other; },
            py::is_operator());
    utils::python::addOutputOp(cls, "__str__");
    utils::python::addOutputOp(cls, "__repr__");
    utils::python::addHash(cls);

    /* Methods */
    cls.def("getSysName", &CppClass::getSysName);
}

WRAP(CameraSys) {
    /* Module level */
    py::class_<CameraSysPrefix> clsCameraSysPrefix(mod, "CameraSysPrefix");
    py::class_<CameraSys> clsCameraSys(mod, "CameraSys");

    // The following must come after the associated pybind11 class is declared
    // (e.g. FOCAL_PLANE is a CameraSys, so clsCameraSys must have been declared
    mod.attr("FOCAL_PLANE") = py::cast(FOCAL_PLANE);
    mod.attr("FIELD_ANGLE") = py::cast(FIELD_ANGLE);
    mod.attr("PIXELS") = py::cast(PIXELS);
    mod.attr("TAN_PIXELS") = py::cast(TAN_PIXELS);
    mod.attr("ACTUAL_PIXELS") = py::cast(ACTUAL_PIXELS);

    /* Member types and enums */
    declareCommonSysMethods<CameraSysPrefix>(clsCameraSysPrefix);
    declareCommonSysMethods<CameraSys>(clsCameraSys);

    /* Constructors */
    clsCameraSysPrefix.def(py::init<std::string const &>(), "sysName"_a);
    clsCameraSys.def(py::init<std::string const &, std::string const &>(), "sysName"_a,
                     "detectorName"_a = "");
    clsCameraSys.def(py::init<CameraSysPrefix const &, std::string const &>(), "sysPrefix"_a,
                     "detectorName"_a = "");

    /* Operators */

    /* Members */
    clsCameraSys.def("getDetectorName", &CameraSys::getDetectorName);
    clsCameraSys.def("hasDetectorName", &CameraSys::hasDetectorName);
}

WRAP(Orientation) {
    //py::module::import("lsst.geom");

    /* Module level */
    py::class_<Orientation> cls(mod, "Orientation");

    /* Member types and enums */

    /* Constructors */
    cls.def(py::init<lsst::geom::Point2D, lsst::geom::Point2D, lsst::geom::Angle, lsst::geom::Angle, lsst::geom::Angle>(),
            "fpPosition"_a = lsst::geom::Point2D(0, 0), "refPoint"_a = lsst::geom::Point2D(-0.5, -0.5),
            "yaw"_a = lsst::geom::Angle(0), "pitch"_a = lsst::geom::Angle(0), "roll"_a = lsst::geom::Angle(0));

    /* Operators */

    /* Members */
    cls.def("getFpPosition", &Orientation::getFpPosition);
    cls.def("getReferencePoint", &Orientation::getReferencePoint);
    cls.def("getYaw", &Orientation::getYaw);
    cls.def("getPitch", &Orientation::getPitch);
    cls.def("getRoll", &Orientation::getRoll);
    cls.def("getNQuarter", &Orientation::getNQuarter);
    cls.def("makePixelFpTransform", &Orientation::makePixelFpTransform, "pixelSizeMm"_a);
    cls.def("makeFpPixelTransform", &Orientation::makeFpPixelTransform, "pixelSizeMm"_a);
    cls.def("getFpPosition", &Orientation::getFpPosition);
    cls.def("getFpPosition", &Orientation::getFpPosition);
    cls.def("getFpPosition", &Orientation::getFpPosition);
    cls.def("getFpPosition", &Orientation::getFpPosition);
}

void declareTransformMap(py::module & mod) {
    PyTransformMap cls(mod, "TransformMap");

    PyTransformMapConnection clsConnection(cls, "Connection");
    clsConnection.def(py::init<std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const>,
                               CameraSys const &, CameraSys const &>(),
                      "transform"_a, "fromSys"_a, "toSys"_a);
    clsConnection.def_readwrite("transform", &TransformMap::Connection::transform);
    clsConnection.def_readwrite("fromSys", &TransformMap::Connection::fromSys);
    clsConnection.def_readwrite("toSys", &TransformMap::Connection::toSys);
    utils::python::addOutputOp(clsConnection, "__repr__");

    cls.def(
        py::init([](
            CameraSys const &reference,
            TransformMap::Transforms const & transforms
        ) {
            // An apparent pybind11 bug: it's usually happy to cast away constness, but won't do it here.
            return std::const_pointer_cast<TransformMap>(TransformMap::make(reference, transforms));
        }),
        "reference"_a, "transforms"_a
    );
    cls.def(
        py::init([](
            CameraSys const &reference,
            std::vector<TransformMap::Connection> const & connections
        ) {
            // An apparent pybind11 bug: it's usually happy to cast away constness, but won't do it here.
            return std::const_pointer_cast<TransformMap>(TransformMap::make(reference, connections));
        }),
        "reference"_a, "connections"_a
    );
    cls.def("__len__", &TransformMap::size);
    cls.def("__contains__", &TransformMap::contains);
    cls.def("__iter__", [](TransformMap const &self) { return py::make_iterator(self.begin(), self.end()); },
            py::keep_alive<0, 1>()); /* Essential: keep object alive while iterator exists */

    cls.def(
        "transform",
        py::overload_cast<lsst::geom::Point2D const &, CameraSys const &, CameraSys const &>(
            &TransformMap::transform,
            py::const_
        ),
        "point"_a, "fromSys"_a, "toSys"_a
    );
    cls.def(
        "transform",
        py::overload_cast<std::vector<lsst::geom::Point2D> const &, CameraSys const &, CameraSys const &>(
            &TransformMap::transform,
            py::const_
        ),
        "pointList"_a, "fromSys"_a, "toSys"_a
    );
    cls.def("getTransform", &TransformMap::getTransform, "fromSys"_a, "toSys"_a);
    cls.def("getConnections", &TransformMap::getConnections);

    table::io::python::addPersistableMethods(cls);

}

WRAP(TransformMap) {
    declareTransformMap(mod);
}

// Declare Detector methods overloaded on one coordinate system class
template <typename SysT, typename PyClass>
void declare1SysMethods(PyClass &cls) {
    cls.def("getCorners",
            (std::vector<lsst::geom::Point2D>(Detector::*)(SysT const &) const) & Detector::getCorners,
            "cameraSys"_a);
    cls.def("getCenter", (lsst::geom::Point2D(Detector::*)(SysT const &) const) & Detector::getCenter,
            "cameraSys"_a);
    cls.def("hasTransform", (bool (Detector::*)(SysT const &) const) & Detector::hasTransform, "cameraSys"_a);

}

// Declare Detector methods templated on two coordinate system classes
template <typename FromSysT, typename ToSysT, typename PyClass>
void declare2SysMethods(PyClass &cls) {
    cls.def("getTransform",
            (std::shared_ptr<geom::TransformPoint2ToPoint2>(Detector::*)(FromSysT const &, ToSysT const &)
                     const) &
                    Detector::getTransform,
            "fromSys"_a, "toSys"_a);
    cls.def("transform",
            (lsst::geom::Point2D(Detector::*)(lsst::geom::Point2D const &, FromSysT const &, ToSysT const &) const) &
                    Detector::transform,
            "point"_a, "fromSys"_a, "toSys"_a);
    cls.def("transform",
            (std::vector<lsst::geom::Point2D>(Detector::*)(std::vector<lsst::geom::Point2D> const &, FromSysT const &,
                                                     ToSysT const &) const) &
                    Detector::transform,
            "points"_a, "fromSys"_a, "toSys"_a);
}

void declareDetectorBase(py::module & mod) {
    PyDetectorBase cls(mod, "DetectorBase");
    cls.def("getName", &DetectorBase::getName);
    cls.def("getId", &DetectorBase::getId);
    cls.def("getType", &DetectorBase::getType);
    cls.def("getPhysicalType", &DetectorBase::getPhysicalType);
    cls.def("getSerial", &DetectorBase::getSerial);
    cls.def("getBBox", &DetectorBase::getBBox);
    cls.def("getOrientation", &DetectorBase::getOrientation);
    cls.def("getPixelSize", &DetectorBase::getPixelSize);
    cls.def("hasCrosstalk", &DetectorBase::hasCrosstalk);
    cls.def("getCrosstalk", &DetectorBase::getCrosstalk);
    cls.def("getNativeCoordSys", &DetectorBase::getNativeCoordSys);
    cls.def("makeCameraSys",
            py::overload_cast<CameraSys const &>(&DetectorBase::makeCameraSys, py::const_),
            "cameraSys"_a);
    cls.def("makeCameraSys",
            py::overload_cast<CameraSysPrefix const &>(&DetectorBase::makeCameraSys, py::const_),
            "cameraSysPrefix"_a);
}

void declareDetectorBuilder(PyDetector & parent);
void declareDetectorPartialRebuilder(PyDetector & parent);
void declareDetectorInCameraBuilder(PyDetector & parent);

void declareDetector(py::module & mod) {
    PyDetector cls(mod, "Detector");
    declareDetectorBuilder(cls);
    declareDetectorPartialRebuilder(cls);
    declareDetectorInCameraBuilder(cls);
    cls.def("rebuild", &Detector::rebuild);
    declare1SysMethods<CameraSys>(cls);
    declare1SysMethods<CameraSysPrefix>(cls);
    declare2SysMethods<CameraSys, CameraSys>(cls);
    declare2SysMethods<CameraSys, CameraSysPrefix>(cls);
    declare2SysMethods<CameraSysPrefix, CameraSys>(cls);
    declare2SysMethods<CameraSysPrefix, CameraSysPrefix>(cls);
    cls.def("getTransformMap", &Detector::getTransformMap);
    cls.def("getAmplifiers", &Detector::getAmplifiers);
    // __iter__ defined in pure-Python extension
    cls.def("__getitem__",
            [](Detector const & self, std::ptrdiff_t i) {
                return self[utils::python::cppIndex(self.size(), i)];
            },
            "i"_a);
    cls.def("__getitem__",
            py::overload_cast<std::string const &>(&Detector::operator[], py::const_),
            "name"_a);
    cls.def("__len__", &Detector::size);
    table::io::python::addPersistableMethods(cls);
}

void declareDetectorBuilder(PyDetector & parent) {
    PyDetectorBuilder cls(parent, "Builder");
    cls.def("setBBox", &Detector::Builder::setBBox);
    cls.def("setType", &Detector::Builder::setType);
    cls.def("setSerial", &Detector::Builder::setSerial);
    cls.def("setPhysicalType", &Detector::Builder::setPhysicalType);
    cls.def("setCrosstalk", &Detector::Builder::setCrosstalk);
    cls.def("unsetCrosstalk", &Detector::Builder::unsetCrosstalk);
    cls.def("getAmplifiers", &Detector::Builder::getAmplifiers);
    // TODO: __iter__ defined in pure-Python extension
    cls.def("__getitem__",
            [](Detector::Builder const & self, std::ptrdiff_t i) {
                return self[utils::python::cppIndex(self.size(), i)];
            },
            "i"_a);
    cls.def("__getitem__",
            py::overload_cast<std::string const &>(&Detector::Builder::operator[], py::const_),
            "name"_a);
    cls.def("append", &Detector::Builder::append);
    cls.def("clear", &Detector::Builder::clear);
    cls.def("__len__", &Detector::Builder::size);
}

void declareDetectorPartialRebuilder(PyDetector & parent) {
    PyDetectorPartialRebuilder cls(parent, "PartialRebuilder");
    cls.def(py::init<Detector const &>(), "detector"_a);
    cls.def("finish", &Detector::PartialRebuilder::finish);
}

void declareDetectorInCameraBuilder(PyDetector & parent) {
    PyDetectorInCameraBuilder cls(parent, "InCameraBuilder");
    cls.def("setOrientation", &Detector::InCameraBuilder::setOrientation);
    cls.def("setPixelSize", &Detector::InCameraBuilder::setPixelSize);
    cls.def("setTransformFromPixelsTo",
            py::overload_cast<CameraSysPrefix const &,
                              std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const>>(
                &Detector::InCameraBuilder::setTransformFromPixelsTo
            ),
            "toSys"_a, "transform"_a);
    cls.def("setTransformFromPixelsTo",
            py::overload_cast<CameraSys const &,
                              std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const>>(
                &Detector::InCameraBuilder::setTransformFromPixelsTo
            ),
            "toSys"_a, "transform"_a);
    cls.def("discardTransformFromPixelsTo",
            py::overload_cast<CameraSysPrefix const &>(
                &Detector::InCameraBuilder::discardTransformFromPixelsTo
            ),
            "toSys"_a);
    cls.def("discardTransformFromPixelsTo",
            py::overload_cast<CameraSys const &>(
                &Detector::InCameraBuilder::discardTransformFromPixelsTo
            ),
            "toSys"_a);
}

WRAP(Detector) {
    py::enum_<DetectorType>(mod, "DetectorType")
            .value("SCIENCE", DetectorType::SCIENCE)
            .value("FOCUS", DetectorType::FOCUS)
            .value("GUIDER", DetectorType::GUIDER)
            .value("WAVEFRONT", DetectorType::WAVEFRONT)
            ;
    declareDetectorBase(mod);
    declareDetector(mod);
}

template <typename T>
using PyDetectorCollectionBase = py::class_<DetectorCollectionBase<T>,
                                            std::shared_ptr<DetectorCollectionBase<T>>>;

using PyDetectorCollection = py::class_<DetectorCollection, DetectorCollectionBase<Detector const>,
                                        std::shared_ptr<DetectorCollection>>;

template <typename T>
void declareDetectorCollectionBase(PyDetectorCollectionBase<T> & cls) {
    cls.def("getNameMap", &DetectorCollectionBase<T>::getNameMap);
    cls.def("getIdMap", &DetectorCollectionBase<T>::getIdMap);
    cls.def("__len__", &DetectorCollectionBase<T>::size);
    cls.def(
        "get",
        py::overload_cast<std::string const &, std::shared_ptr<T>>(
            &DetectorCollectionBase<T>::get, py::const_
        ),
        "name"_a, "default"_a=nullptr
    );
    cls.def(
        "get",
        py::overload_cast<int, std::shared_ptr<T>>(
            &DetectorCollectionBase<T>::get, py::const_
        ),
        "id"_a, "default"_a=nullptr
    );
    cls.def(
        "__contains__",
        [](DetectorCollectionBase<T> const & self, std::string const &name) {
            return self.get(name) != nullptr;
        }
    );
    cls.def(
        "__contains__",
        [](DetectorCollectionBase<T> const & self, int id) {
            return self.get(id) != nullptr;
        }
    );
}

void declareDetectorCollection(py::module & mod) {
    PyDetectorCollectionBase<Detector const> base(mod, "DetectorCollectionDetectorBase");
    declareDetectorCollectionBase(base);
    PyDetectorCollection cls(mod, "DetectorCollection");
    cls.def(py::init<DetectorCollection::List>());
    cls.def("getFpBBox", &DetectorCollection::getFpBBox);
    table::io::python::addPersistableMethods(cls);
}


WRAP(DetectorCollection) {
    declareDetectorCollection(mod);

    PyDetectorCollectionBase<Detector::InCameraBuilder> cameraBuilderBase(
        mod,
        "DetectorCollectionBuilderBase"
    );
    declareDetectorCollectionBase(cameraBuilderBase);
}
} // namespace

WRAP(CameraGeom) {
    auto cameramod = mod.def_submodule("cameraGeom");
    wrapTransformMap(cameramod);
    wrapDetector(cameramod);
    wrapDetectorCollection(cameramod);
    wrapCamera(cameramod);
    wrapAmplifier(cameramod);
    wrapCameraSys(cameramod);
    wrapOrientation(cameramod);
}
}
}
}
