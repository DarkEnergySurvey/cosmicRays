#include "pybind/meas_bind.h"
#include "pybind11/eigen.h"

#include <memory>
#include "lsst/meas/base/Transform.h"
#include "lsst/meas/base/SincCoeffs.h"
#include "ndarray/pybind11.h"

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/afw/table/FunctorKey.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/SdssShape.h"
#include "lsst/meas/base/SdssCentroid.h"
#include "lsst/meas/base/ScaledApertureFlux.h"
#include "lsst/meas/base/PsfFlux.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/PixelFlags.h"
#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"
#include "lsst/meas/base/PeakLikelihoodFlux.h"
#include "lsst/meas/base/LocalBackground.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/GaussianFlux.h"
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/utils/python.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/exceptions/python/Exception.h"
#include "lsst/meas/base/CircularApertureFlux.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/Blendedness.h"
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/NaiveCentroid.h"

namespace lsst {
namespace meas {
namespace base {
namespace sinc {
namespace {
template <typename T>
void declareSincCoeffs(py::module& mod, std::string const& suffix) {
    py::class_<SincCoeffs<T> >(mod, ("SincCoeffs" + suffix).c_str())
            .def_static("cache", &SincCoeffs<T>::cache, "rInner"_a, "rOuter"_a)
            .def_static("get", &SincCoeffs<T>::get, "outerEllipse"_a, "innerRadiusFactor"_a);
}
WRAP(SincCoeffs) {
    declareSincCoeffs<float>(mod, "F");
    declareSincCoeffs<double>(mod, "D");
}

}
}

namespace shape {
namespace {
void declareShapeResult(py::module &mod) {
    py::class_<ShapeResult, std::shared_ptr<ShapeResult> >(mod, "ShapeResult")
            .def(py::init<>())
            .def(py::init<ShapeElement, ShapeElement, ShapeElement, ShapeCov const &>(), "xx"_a, "yy"_a, "xy"_a,
                 "matrix"_a)
            .def(py::init<ShapeElement, ShapeElement, ShapeElement, ErrElement, ErrElement, ErrElement>(), "xx"_a,
                 "yy"_a, "xy"_a, "xxErr"_a, "yyErr"_a, "xyErr"_a)

            .def("getShape", &ShapeResult::getShape)
            .def("getQuadrupole", &ShapeResult::getQuadrupole)
            .def("setShape", &ShapeResult::setShape, "shape"_a)
            .def("getShapeErr", &ShapeResult::getShapeErr)
            .def("setShapeErr", (void (ShapeResult::*)(ShapeCov const &)) & ShapeResult::setShapeErr, "matrix"_a)
            .def("setShapeErr",
                 (void (ShapeResult::*)(ErrElement, ErrElement, ErrElement)) & ShapeResult::setShapeErr,
                 "xxErr"_a, "yyErr"_a, "xyErr"_a)

            .def_readwrite("xx", &ShapeResult::xx)
            .def_readwrite("yy", &ShapeResult::yy)
            .def_readwrite("xy", &ShapeResult::xy)
            .def_readwrite("xxErr", &ShapeResult::xxErr)
            .def_readwrite("yyErr", &ShapeResult::yyErr)
            .def_readwrite("xyErr", &ShapeResult::xyErr)
            .def_readwrite("xx_yy_Cov", &ShapeResult::xx_yy_Cov)
            .def_readwrite("xx_xy_Cov", &ShapeResult::xx_xy_Cov)
            .def_readwrite("yy_xy_Cov", &ShapeResult::yy_xy_Cov);
}

void declareShapeResultKey(py::module &mod) {
    py::class_<ShapeResultKey, std::shared_ptr<ShapeResultKey> >(mod, "ShapeResultKey")
            .def_static("addFields", &ShapeResultKey::addFields, "schema"_a, "name"_a, "doc"_a, "uncertainty"_a,
                   "coordType"_a = afw::table::CoordinateType::PIXEL)

            .def(py::init<>())
            .def(py::init<afw::table::QuadrupoleKey const &,
                 afw::table::CovarianceMatrixKey<ErrElement, 3> const &>(),
                 "shape"_a, "shapeErr"_a)
            .def(py::init<afw::table::SubSchema const &>(), "subSchema"_a)

            .def("__eq__", &ShapeResultKey::operator==, py::is_operator())
            .def("__ne__", &ShapeResultKey::operator!=, py::is_operator())

            .def("get", &ShapeResultKey::get, "record"_a)
            .def("set", &ShapeResultKey::set, "record"_a, "value"_a)
            .def("isValid", &ShapeResultKey::isValid)
            .def("getShape", &ShapeResultKey::getShape)
            .def("getShapeErr", &ShapeResultKey::getShapeErr)
            .def("getIxx", &ShapeResultKey::getIxx)
            .def("getIyy", &ShapeResultKey::getIyy)
            .def("getIxy", &ShapeResultKey::getIxy);
}

WRAP(ShapeUtilities) {
    declareShapeResult(mod);
    declareShapeResultKey(mod);

    mod.def("makeShapeTransformMatrix", &makeShapeTransformMatrix, "xform"_a);
}

}
}

namespace sdssShape {
namespace {
using PyShapeControl = py::class_<SdssShapeControl>;

using PyShapeAlgorithm = py::class_<SdssShapeAlgorithm, std::shared_ptr<SdssShapeAlgorithm>, SimpleAlgorithm>;
using PyShapeTransform = py::class_<SdssShapeTransform, std::shared_ptr<SdssShapeTransform>, BaseTransform>;

PyShapeControl declareShapeControl(py::module &mod) {
    PyShapeControl cls(mod, "SdssShapeControl");

    LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, background);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, maxIter);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, maxShift);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, tol1);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, tol2);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssShapeControl, doMeasurePsf);

    cls.def(py::init<>());

    return cls;
}

void declareShapeResultKey(py::module &mod) {
    py::class_<SdssShapeResultKey, std::shared_ptr<SdssShapeResultKey> >(mod, "SdssShapeResultKey")

            // TODO decide whether to wrap default constructor and do it or document why not
            .def(py::init<afw::table::SubSchema const &>(), "subSchema"_a)

            .def_static("addFields", &FluxResultKey::addFields, "schema"_a, "name"_a, "doMeasurePsf"_a)

            .def("__eq__", &SdssShapeResultKey::operator==, py::is_operator())
            .def("__ne__", &SdssShapeResultKey::operator!=, py::is_operator())

            .def("get", &SdssShapeResultKey::get, "record"_a)
            .def("set", &SdssShapeResultKey::set, "record"_a, "value"_a)
            .def("getPsfShape", &SdssShapeResultKey::getPsfShape, "record"_a)
            .def("setPsfShape", &SdssShapeResultKey::setPsfShape, "record"_a, "value"_a)
            .def("isValid", &SdssShapeResultKey::isValid)
            .def("getFlagHandler", &SdssShapeResultKey::getFlagHandler);
}

template <typename ImageT>
static void declareComputeMethods(PyShapeAlgorithm &cls) {
    cls.def_static(
            "computeAdaptiveMoments",
            (SdssShapeResult(*)(ImageT const &, geom::Point2D const &, bool, SdssShapeControl const &)) &
                    SdssShapeAlgorithm::computeAdaptiveMoments,
            "image"_a, "position"_a, "negative"_a = false, "ctrl"_a = SdssShapeControl());
    cls.def_static(
            "computeFixedMomentsFlux",
            (FluxResult(*)(ImageT const &, afw::geom::ellipses::Quadrupole const &, geom::Point2D const &)) &
                    SdssShapeAlgorithm::computeFixedMomentsFlux,
            "image"_a, "shape"_a, "position"_a);
}

PyShapeAlgorithm declareShapeAlgorithm(py::module &mod) {
    PyShapeAlgorithm cls(mod, "SdssShapeAlgorithm");

    cls.attr("FAILURE") = py::cast(SdssShapeAlgorithm::FAILURE);
    cls.attr("UNWEIGHTED_BAD") = py::cast(SdssShapeAlgorithm::UNWEIGHTED_BAD);
    cls.attr("UNWEIGHTED") = py::cast(SdssShapeAlgorithm::UNWEIGHTED);
    cls.attr("SHIFT") = py::cast(SdssShapeAlgorithm::SHIFT);
    cls.attr("MAXITER") = py::cast(SdssShapeAlgorithm::MAXITER);
    cls.attr("PSF_SHAPE_BAD") = py::cast(SdssShapeAlgorithm::PSF_SHAPE_BAD);

    cls.def(py::init<SdssShapeAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    declareComputeMethods<afw::image::Image<int>>(cls);
    declareComputeMethods<afw::image::Image<float>>(cls);
    declareComputeMethods<afw::image::Image<double>>(cls);
    declareComputeMethods<afw::image::MaskedImage<int>>(cls);
    declareComputeMethods<afw::image::MaskedImage<float>>(cls);
    declareComputeMethods<afw::image::MaskedImage<double>>(cls);

    cls.def("measure", &SdssShapeAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &SdssShapeAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

void declareShapeResult(py::module &mod) {
    py::class_<SdssShapeResult, std::shared_ptr<SdssShapeResult>, ShapeResult,
                                     CentroidResult, FluxResult>(mod, "SdssShapeResult")

            .def(py::init<>())

            .def_readwrite("instFlux_xx_Cov", &SdssShapeResult::instFlux_xx_Cov)
            .def_readwrite("instFlux_yy_Cov", &SdssShapeResult::instFlux_yy_Cov)
            .def_readwrite("instFlux_xy_Cov", &SdssShapeResult::instFlux_xy_Cov)
            .def_readwrite("flags", &SdssShapeResult::flags)

            // TODO this method says it's a workaround for Swig which doesn't understand std::bitset
            .def("getFlag", (bool (SdssShapeResult::*)(unsigned int) const) & SdssShapeResult::getFlag, "index"_a)
            .def("getFlag", (bool (SdssShapeResult::*)(std::string const &name) const) & SdssShapeResult::getFlag,
                 "name"_a);
}

PyShapeTransform declareShapeTransform(py::module &mod) {
    PyShapeTransform cls(mod, "SdssShapeTransform");

    cls.def(py::init<SdssShapeTransform::Control const &, std::string const &, afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    cls.def("__call__", &SdssShapeTransform::operator(), "inputCatalog"_a, "outputCatalog"_a, "wcs"_a,
            "photoCalib"_a);

    return cls;
}

WRAP(SdssShape) {
    auto clsShapeControl = declareShapeControl(mod);
    declareShapeResultKey(mod);
    auto clsShapeAlgorithm = declareShapeAlgorithm(mod);
    declareShapeResult(mod);
    auto clsShapeTransform = declareShapeTransform(mod);

    clsShapeAlgorithm.attr("Control") = clsShapeControl;
    clsShapeTransform.attr("Control") = clsShapeControl;

    python::declareAlgorithm<SdssShapeAlgorithm, SdssShapeControl, SdssShapeTransform>(
            clsShapeAlgorithm, clsShapeControl, clsShapeTransform);
}

}
}

namespace sdssCentroid {
namespace {
using PyCentroidAlgorithm =
        py::class_<SdssCentroidAlgorithm, std::shared_ptr<SdssCentroidAlgorithm>, SimpleAlgorithm>;
using PyCentroidControl = py::class_<SdssCentroidControl>;
using PyCentroidTransform =
        py::class_<SdssCentroidTransform, std::shared_ptr<SdssCentroidTransform>, BaseTransform>;

PyCentroidControl declareCentroidControl(py::module &mod) {
    PyCentroidControl cls(mod, "SdssCentroidControl");

    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, binmax);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, peakMin);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, wfac);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, doFootprintCheck);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, maxDistToPeak);

    cls.def(py::init<>());

    return cls;
}

PyCentroidAlgorithm declareCentroidAlgorithm(py::module &mod) {
    PyCentroidAlgorithm cls(mod, "SdssCentroidAlgorithm");

    cls.attr("FAILURE") = py::cast(SdssCentroidAlgorithm::FAILURE);
    cls.attr("EDGE") = py::cast(SdssCentroidAlgorithm::EDGE);
    cls.attr("NO_SECOND_DERIVATIVE") = py::cast(SdssCentroidAlgorithm::NO_SECOND_DERIVATIVE);
    cls.attr("ALMOST_NO_SECOND_DERIVATIVE") = py::cast(SdssCentroidAlgorithm::ALMOST_NO_SECOND_DERIVATIVE);
    cls.attr("NOT_AT_MAXIMUM") = py::cast(SdssCentroidAlgorithm::NOT_AT_MAXIMUM);

    cls.def(py::init<SdssCentroidAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("measure", &SdssCentroidAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &SdssCentroidAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyCentroidTransform declareCentroidTransform(py::module &mod) {
    PyCentroidTransform cls(mod, "SdssCentroidTransform");

    cls.def(py::init<SdssCentroidTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

WRAP(SdssCentroid) {
    auto clsCentroidControl = declareCentroidControl(mod);
    auto clsCentroidAlgorithm = declareCentroidAlgorithm(mod);
    auto clsCentroidTransform = declareCentroidTransform(mod);

    clsCentroidAlgorithm.attr("Control") = clsCentroidControl;
    clsCentroidTransform.attr("Control") = clsCentroidControl;

    python::declareAlgorithm<SdssCentroidAlgorithm, SdssCentroidControl, SdssCentroidTransform>(
            clsCentroidAlgorithm, clsCentroidControl, clsCentroidTransform);
}

}
}

namespace aflux {
namespace {
using PyFluxAlgorithm = py::class_<ScaledApertureFluxAlgorithm, std::shared_ptr<ScaledApertureFluxAlgorithm>,
                                   SimpleAlgorithm>;
using PyFluxControl = py::class_<ScaledApertureFluxControl>;
using PyFluxTransform =
        py::class_<ScaledApertureFluxTransform, std::shared_ptr<ScaledApertureFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "ScaledApertureFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, ScaledApertureFluxControl, scale);
    LSST_DECLARE_CONTROL_FIELD(cls, ScaledApertureFluxControl, shiftKernel);

    cls.def(py::init<>());

    return cls;
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "ScaledApertureFluxAlgorithm");

    cls.def(py::init<ScaledApertureFluxAlgorithm::Control const &, std::string const &,
                     afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("measure", &ScaledApertureFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &ScaledApertureFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "ScaledApertureFluxTransform");

    cls.def(py::init<ScaledApertureFluxTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

WRAP(ScaledApertureFlux) {
    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<ScaledApertureFluxAlgorithm, ScaledApertureFluxControl,
                             ScaledApertureFluxTransform>(clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}

}
}

namespace pflux {
namespace {

using PyFluxAlgorithm = py::class_<PsfFluxAlgorithm, std::shared_ptr<PsfFluxAlgorithm>, SimpleAlgorithm>;
using PyFluxControl = py::class_<PsfFluxControl>;
using PyFluxTransform = py::class_<PsfFluxTransform, std::shared_ptr<PsfFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "PsfFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, PsfFluxControl, badMaskPlanes);

    cls.def(py::init<>());

    return cls;
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "PsfFluxAlgorithm");

    cls.attr("FAILURE") = py::cast(PsfFluxAlgorithm::FAILURE);
    cls.attr("NO_GOOD_PIXELS") = py::cast(PsfFluxAlgorithm::NO_GOOD_PIXELS);
    cls.attr("EDGE") = py::cast(PsfFluxAlgorithm::EDGE);

    cls.def(py::init<PsfFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def(py::init<PsfFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &,
                     std::string const &>(),
            "ctrl"_a, "name"_a, "schema"_a, "logName"_a);
    return cls;
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "PsfFluxTransform");

    cls.def(py::init<PsfFluxTransform::Control const &, std::string const &, afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

WRAP(PsfFlux) {
    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<PsfFluxAlgorithm, PsfFluxControl, PsfFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}


}
}

namespace pkflux {
namespace {

using PyFluxAlgorithm = py::class_<PeakLikelihoodFluxAlgorithm, std::shared_ptr<PeakLikelihoodFluxAlgorithm>,
                                   SimpleAlgorithm>;
using PyFluxControl = py::class_<PeakLikelihoodFluxControl>;
using PyFluxTransform =
        py::class_<PeakLikelihoodFluxTransform, std::shared_ptr<PeakLikelihoodFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "PeakLikelihoodFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, PeakLikelihoodFluxControl, warpingKernelName);

    return cls;
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "PeakLikelihoodFluxAlgorithm");

    cls.attr("FAILURE") = py::cast(PeakLikelihoodFluxAlgorithm::FAILURE);

    cls.def(py::init<PeakLikelihoodFluxAlgorithm::Control const &, std::string const &,
                     afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("measure", &PeakLikelihoodFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &PeakLikelihoodFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "PeakLikelihoodFluxTransform");

    cls.def(py::init<PeakLikelihoodFluxTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

WRAP(PeakLikelihoodFlux) {
    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<PeakLikelihoodFluxAlgorithm, PeakLikelihoodFluxControl,
                             PeakLikelihoodFluxTransform>(clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}

}
}

namespace background {
namespace {

using PyAlgorithm =
        py::class_<LocalBackgroundAlgorithm, std::shared_ptr<LocalBackgroundAlgorithm>, SimpleAlgorithm>;
using PyControl = py::class_<LocalBackgroundControl>;
using PyTransform =
        py::class_<LocalBackgroundTransform, std::shared_ptr<LocalBackgroundTransform>, BaseTransform>;

PyControl declareControl(py::module &mod) {
    PyControl cls(mod, "LocalBackgroundControl");

    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, badMaskPlanes);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, annulusInner);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, annulusOuter);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, bgRej);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, bgIter);

    cls.def(py::init<>());

    return cls;
}

PyAlgorithm declareAlgorithm(py::module &mod) {
    PyAlgorithm cls(mod, "LocalBackgroundAlgorithm");

    cls.attr("FAILURE") = py::cast(LocalBackgroundAlgorithm::FAILURE);
    cls.attr("NO_GOOD_PIXELS") = py::cast(LocalBackgroundAlgorithm::NO_GOOD_PIXELS);

    cls.def(py::init<LocalBackgroundAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def(py::init<LocalBackgroundAlgorithm::Control const &, std::string const &, afw::table::Schema &,
                     std::string const &>(),
            "ctrl"_a, "name"_a, "schema"_a, "logName"_a);
    return cls;
}

PyTransform declareTransform(py::module &mod) {
    PyTransform cls(mod, "LocalBackgroundTransform");

    cls.def(py::init<LocalBackgroundTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

WRAP(LocalBackground) {
    auto clsControl = declareControl(mod);
    auto clsAlgorithm = declareAlgorithm(mod);
    auto clsTransform = declareTransform(mod);

    clsAlgorithm.attr("Control") = clsControl;
    clsTransform.attr("Control") = clsControl;

    python::declareAlgorithm<LocalBackgroundAlgorithm, LocalBackgroundControl, LocalBackgroundTransform>(
            clsAlgorithm, clsControl, clsTransform);
}

}
}

namespace gflux {
namespace {

using PyFluxAlgorithm =
        py::class_<GaussianFluxAlgorithm, std::shared_ptr<GaussianFluxAlgorithm>, SimpleAlgorithm>;
using PyFluxControl = py::class_<GaussianFluxControl>;
using PyFluxTransform =
        py::class_<GaussianFluxTransform, std::shared_ptr<GaussianFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "GaussianFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, GaussianFluxControl, background);

    return cls;
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "GaussianFluxAlgorithm");

    cls.def(py::init<GaussianFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.attr("FAILURE") = py::cast(GaussianFluxAlgorithm::FAILURE);

    cls.def("measure", &GaussianFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &GaussianFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "GaussianFluxTransform");

    cls.def(py::init<GaussianFluxTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

WRAP(GaussianFlux) {
    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<GaussianFluxAlgorithm, GaussianFluxControl, GaussianFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}

}
}

namespace fluxUtils {
namespace {

using PyFluxResult = py::class_<FluxResult, std::shared_ptr<FluxResult>>;
using PyFluxResultKey = py::class_<FluxResultKey, std::shared_ptr<FluxResultKey>>;
using PyMagResult = py::class_<MagResult, std::shared_ptr<MagResult>>;
using PyMagResultKey = py::class_<MagResultKey, std::shared_ptr<MagResultKey>>;

void declareFluxResult(py::module &mod) {
    PyFluxResult cls(mod, "FluxResult");

    cls.def_readwrite("instFlux", &FluxResult::instFlux);
    cls.def_readwrite("instFluxErr", &FluxResult::instFluxErr);
}

void declareFluxResultKey(py::module &mod) {
    PyFluxResultKey cls(mod, "FluxResultKey");

    cls.def(py::init<>());
    cls.def(py::init<afw::table::Key<meas::base::Flux> const &,
                     afw::table::Key<meas::base::FluxErrElement> const &>(),
            "instFlux"_a, "instFluxErr"_a);
    cls.def(py::init<afw::table::SubSchema const &>());

    cls.def("__eq__", &FluxResultKey::operator==, py::is_operator());
    cls.def("__ne__", &FluxResultKey::operator!=, py::is_operator());

    cls.def("get", &FluxResultKey::get);
    cls.def("set", &FluxResultKey::set);
    cls.def_static("addFields", &FluxResultKey::addFields, "schema"_a, "name"_a, "doc"_a);
    cls.def("isValid", &FluxResultKey::isValid);
    cls.def("getInstFlux", &FluxResultKey::getInstFlux);
    cls.def("getInstFluxErr", &FluxResultKey::getInstFluxErr);
}

void declareMagResult(py::module &mod) {
    PyMagResult cls(mod, "MagResult");

    cls.def_readwrite("mag", &MagResult::mag);
    cls.def_readwrite("magErr", &MagResult::magErr);
}

void declareMagResultKey(py::module &mod) {
    PyMagResultKey cls(mod, "MagResultKey");

    cls.def(py::init<>());
    cls.def(py::init<afw::table::SubSchema const &>());

    cls.def("get", &MagResultKey::get);
    cls.def("set",
            (void (MagResultKey::*)(afw::table::BaseRecord &, MagResult const &) const) & MagResultKey::set);
    cls.def("set", (void (MagResultKey::*)(afw::table::BaseRecord &, afw::image::Measurement const &) const) &
                           MagResultKey::set);
    cls.def_static("addFields", &MagResultKey::addFields, "schema"_a, "name"_a);
}

WRAP(FluxUtilities) {
    declareFluxResult(mod);
    declareFluxResultKey(mod);
    declareMagResult(mod);
    declareMagResultKey(mod);
}

}
}

namespace flags {
namespace {

void declareFlagDefinition(py::module &mod) {
    py::class_<FlagDefinition, std::shared_ptr<FlagDefinition>> cls(mod, "FlagDefinition");

    cls.def(py::init<>());
    cls.def(py::init<std::string, std::string, std::size_t>(), "name"_a, "doc"_a,
            "number"_a = FlagDefinition::number_undefined);

    cls.def_readwrite("name", &FlagDefinition::name);
    cls.def_readwrite("doc", &FlagDefinition::doc);
    cls.def_readwrite("number", &FlagDefinition::number);

    cls.def("__eq__", &FlagDefinition::operator==, py::is_operator());
    cls.def("__ne__", &FlagDefinition::operator!=, py::is_operator());
}

void declareFlagDefinitionList(py::module &mod) {
    py::class_<FlagDefinitionList, std::shared_ptr<FlagDefinitionList>> cls(mod, "FlagDefinitionList");

    cls.def(py::init<>());
    cls.def(py::init<std::initializer_list<FlagDefinition> const &>());

    cls.def("__getitem__", [](FlagDefinitionList const &self, int i) {
        try {
            auto cind = utils::python::cppIndex(self.size(), i);
            return self[cind];
        } catch (pex::exceptions::OutOfRangeError &err) {
            // Python needs exception to be IndexError to generate __iter__; see DM-9715
            PyErr_SetString(PyExc_IndexError, err.what());
            throw py::error_already_set();
        }
    });
    cls.def("__len__", &FlagDefinitionList::size);

    cls.def("getEmptyList", &FlagDefinitionList::getEmptyList);
    cls.def("getDefinition",
            (FlagDefinition(FlagDefinitionList::*)(std::size_t) const) & FlagDefinitionList::getDefinition,
            "index"_a);
    cls.def("getDefinition",
            (FlagDefinition(FlagDefinitionList::*)(std::string const &) const) &
                    FlagDefinitionList::getDefinition,
            "name"_a);
    cls.def("hasDefinition", &FlagDefinitionList::hasDefinition, "name"_a);
    cls.def("addFailureFlag", &FlagDefinitionList::addFailureFlag, "doc"_a = "General Failure Flag");
    cls.def("add", &FlagDefinitionList::add, "name"_a, "doc"_a);
}

void declareFlagHandler(py::module &mod) {
    py::class_<FlagHandler, std::shared_ptr<FlagHandler>> cls(mod, "FlagHandler");

    cls.def(py::init<>());
    cls.def(py::init<afw::table::SubSchema const &, FlagDefinitionList const &, FlagDefinitionList const &>(),
            "s"_a, "flagDefs"_a, "exclDefs"_a = FlagDefinitionList::getEmptyList());

    cls.def_static("getFailureFlagName", &FlagHandler::getFailureFlagName);
    cls.def_static("addFields", &FlagHandler::addFields, "schema"_a, "prefix"_a, "flagDefs"_a,
                   "exclDefs"_a = FlagDefinitionList::getEmptyList());

    cls.def("getFlagNumber", &FlagHandler::getFlagNumber, "flagName"_a);
    cls.def("getFlagName", &FlagHandler::getFlagName, "i"_a);
    cls.def("getValue",
            (bool (FlagHandler::*)(afw::table::BaseRecord const &, std::size_t) const) &
                    FlagHandler::getValue,
            "record"_a, "i"_a);
    cls.def("getValue",
            (bool (FlagHandler::*)(afw::table::BaseRecord const &, std::string const &) const) &
                    FlagHandler::getValue,
            "record"_a, "flagName"_a);
    cls.def("setValue",
            (void (FlagHandler::*)(afw::table::BaseRecord &, std::size_t, bool) const) &
                    FlagHandler::setValue,
            "record"_a, "i"_a, "value"_a);
    cls.def("setValue",
            (void (FlagHandler::*)(afw::table::BaseRecord &, std::string const &, bool) const) &
                    FlagHandler::setValue,
            "record"_a, "flagName"_a, "value"_a);
    cls.def("getFailureFlagNumber", &FlagHandler::getFailureFlagNumber);
    cls.def("handleFailure", &FlagHandler::handleFailure, "record"_a, "error"_a = nullptr);
}

WRAP(FlagHandler) {
    declareFlagDefinition(mod);
    declareFlagDefinitionList(mod);
    declareFlagHandler(mod);
}

}
}

namespace centroidUtils {
namespace {

void declareCentroidResult(py::module &mod) {
    py::class_<CentroidResult, std::shared_ptr<CentroidResult>> cls(mod, "CentroidResult");

    cls.def_readwrite("x", &CentroidResult::x);
    cls.def_readwrite("y", &CentroidResult::y);
    cls.def_readwrite("xErr", &CentroidResult::xErr);
    cls.def_readwrite("yErr", &CentroidResult::yErr);
    cls.def_readwrite("x_y_Cov", &CentroidResult::x_y_Cov);

    cls.def(py::init<>());
    cls.def(py::init<CentroidElement, CentroidElement, CentroidCov const &>(), "x"_a, "y"_a, "matrix"_a);
    cls.def(py::init<CentroidElement, CentroidElement, ErrElement, ErrElement>(), "x"_a, "y"_a, "xErr"_a,
            "yErr"_a);

    cls.def("getCentroid", &CentroidResult::getCentroid);
    cls.def("setCentroid", &CentroidResult::setCentroid, "centroid"_a);
    cls.def("getPoint", &CentroidResult::getPoint);
    cls.def("getCentroidErr", &CentroidResult::getCentroidErr);
    cls.def("setCentroidErr",
            (void (CentroidResult::*)(CentroidCov const &)) & CentroidResult::setCentroidErr, "matrix"_a);
    cls.def("setCentroidErr",
            (void (CentroidResult::*)(ErrElement, ErrElement)) & CentroidResult::setCentroidErr, "xErr"_a,
            "yErr"_a);
}

void declareCentroidResultKey(py::module &mod) {
    py::class_<CentroidResultKey> cls(mod, "CentroidResultKey");

    cls.def(py::init<>());
    cls.def(py::init<afw::table::PointKey<CentroidElement> const &,
                     afw::table::CovarianceMatrixKey<ErrElement, 2> const &>(),
            "centroid"_a, "uncertainty"_a);
    cls.def(py::init<afw::table::SubSchema const &>(), "subSchema"_a);

    cls.def_static("addFields", &CentroidResultKey::addFields, "schema"_a, "name"_a, "doc"_a, "uncertainty"_a);

    cls.def("__eq__", &CentroidResultKey::operator==, py::is_operator());
    cls.def("__nq__", &CentroidResultKey::operator!=, py::is_operator());

    cls.def("get", &CentroidResultKey::get, "record"_a);
    cls.def("set", &CentroidResultKey::set, "record"_a, "value"_a);
    cls.def("isValid", &CentroidResultKey::isValid);
    cls.def("getCentroid", &CentroidResultKey::getCentroid);
    cls.def("getCentroidErr", &CentroidResultKey::getCentroidErr);
    cls.def("getX", &CentroidResultKey::getX);
    cls.def("getY", &CentroidResultKey::getY);
}

void declareUCentroidTransform(py::module &mod) {
    py::class_<CentroidTransform, std::shared_ptr<CentroidTransform>, BaseTransform> cls(mod, "CentroidTransform");

    cls.def(py::init<std::string const &, afw::table::SchemaMapper &>(), "name"_a, "mapper"_a);

    cls.def("__call__", &CentroidTransform::operator(), "inputCatalog"_a, "outputCatalog"_a, "wcs"_a,
            "photoCalib"_a);
}

void declareCentroidChecker(py::module &mod) {
    py::class_<CentroidChecker> cls(mod, "CentroidChecker");

    cls.def(py::init<afw::table::Schema &, std::string const &, bool, double>(), "schema"_a, "name"_a,
            "inside"_a = true, "maxDistFromPeak"_a = -1.0);

    cls.def("__call__", &CentroidChecker::operator(), "record"_a);
}

void declareUncertaintyEnum(py::module &mod) {
    py::enum_<UncertaintyEnum> enm(mod, "UncertaintyEnum");

    enm.value("NO_UNCERTAINTY", UncertaintyEnum::NO_UNCERTAINTY);
    enm.value("SIGMA_ONLY", UncertaintyEnum::SIGMA_ONLY);
    enm.value("FULL_COVARIANCE", UncertaintyEnum::FULL_COVARIANCE);
}

WRAP(CentroidUtilities) {
    declareCentroidResult(mod);
    declareCentroidResultKey(mod);
    declareUCentroidTransform(mod);
    declareCentroidChecker(mod);
    declareUncertaintyEnum(mod);
}

}
}

namespace blend {
namespace {

using PyBlendenessAlgorithm =
        py::class_<BlendednessAlgorithm, std::shared_ptr<BlendednessAlgorithm>, SimpleAlgorithm>;
using PyBlendenessControl = py::class_<BlendednessControl>;

PyBlendenessControl declareBlendednessControl(py::module &mod) {
    PyBlendenessControl cls(mod, "BlendednessControl");

    LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, doOld);
    LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, doFlux);
    LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, doShape);
    LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, nSigmaWeightMax);

    cls.def(py::init<>());

    return cls;
}

PyBlendenessAlgorithm declareBlendednessAlgorithm(py::module &mod) {
    PyBlendenessAlgorithm cls(mod, "BlendednessAlgorithm");

    cls.def(py::init<BlendednessAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.attr("FAILURE") = py::cast(BlendednessAlgorithm::FAILURE);
    cls.attr("NO_CENTROID") = py::cast(BlendednessAlgorithm::NO_CENTROID);
    cls.attr("NO_SHAPE") = py::cast(BlendednessAlgorithm::NO_SHAPE);

    cls.def_static("computeAbsExpectation", &BlendednessAlgorithm::computeAbsExpectation, "data"_a,
                   "variance"_a);
    cls.def_static("computeAbsBias", &BlendednessAlgorithm::computeAbsBias, "mu"_a, "variance"_a);
    cls.def("measureChildPixels", &BlendednessAlgorithm::measureChildPixels, "image"_a, "child"_a);
    cls.def("measureParentPixels", &BlendednessAlgorithm::measureParentPixels, "image"_a, "child"_a);
    cls.def("measure", &BlendednessAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &BlendednessAlgorithm::measure, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

WRAP(Blendedness) {

    auto clsBlendednessControl = declareBlendednessControl(mod);
    auto clsBlendednessAlgorithm = declareBlendednessAlgorithm(mod);

    clsBlendednessAlgorithm.attr("Control") = clsBlendednessControl;

    python::declareAlgorithm<BlendednessAlgorithm, BlendednessControl>(clsBlendednessAlgorithm,
                                                                       clsBlendednessControl);
}

}
}

namespace apflux {
namespace {

using PyFluxAlgorithm =
        py::class_<ApertureFluxAlgorithm, std::shared_ptr<ApertureFluxAlgorithm>, SimpleAlgorithm>;
using PyFluxControl = py::class_<ApertureFluxControl>;
using PyFluxResult = py::class_<ApertureFluxResult, std::shared_ptr<ApertureFluxResult>, FluxResult>;
using PyFluxTransform =
        py::class_<ApertureFluxTransform, std::shared_ptr<ApertureFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "ApertureFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, ApertureFluxControl, radii);
    LSST_DECLARE_CONTROL_FIELD(cls, ApertureFluxControl, maxSincRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, ApertureFluxControl, shiftKernel);

    cls.def(py::init<>());

    return cls;
}

template <typename Image, class PyClass>
void declareComputeFluxes(PyClass &cls) {
    using Control = ApertureFluxAlgorithm::Control;
    using Result = ApertureFluxAlgorithm::Result;

    cls.def_static("computeSincFlux",
                   (Result(*)(Image const &, afw::geom::ellipses::Ellipse const &, Control const &)) &
                           ApertureFluxAlgorithm::computeSincFlux,
                   "image"_a, "ellipse"_a, "ctrl"_a = Control());
    cls.def_static("computeNaiveFlux",
                   (Result(*)(Image const &, afw::geom::ellipses::Ellipse const &, Control const &)) &
                           ApertureFluxAlgorithm::computeNaiveFlux,
                   "image"_a, "ellipse"_a, "ctrl"_a = Control());
    cls.def_static("computeFlux",
                   (Result(*)(Image const &, afw::geom::ellipses::Ellipse const &, Control const &)) &
                           ApertureFluxAlgorithm::computeFlux,
                   "image"_a, "ellipse"_a, "ctrl"_a = Control());
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "ApertureFluxAlgorithm");

    cls.attr("FAILURE") = py::cast(ApertureFluxAlgorithm::FAILURE);
    cls.attr("APERTURE_TRUNCATED") = py::cast(ApertureFluxAlgorithm::APERTURE_TRUNCATED);
    cls.attr("SINC_COEFFS_TRUNCATED") = py::cast(ApertureFluxAlgorithm::SINC_COEFFS_TRUNCATED);

    // constructor not wrapped because class is abstract

    declareComputeFluxes<afw::image::Image<double>>(cls);
    declareComputeFluxes<afw::image::MaskedImage<double>>(cls);
    declareComputeFluxes<afw::image::Image<float>>(cls);
    declareComputeFluxes<afw::image::MaskedImage<float>>(cls);

    cls.def("measure", &ApertureFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &ApertureFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);
    cls.def_static("makeFieldPrefix", &ApertureFluxAlgorithm::makeFieldPrefix, "name"_a, "radius"_a);

    return cls;
}

void declareAFluxResult(py::module &mod) {
    PyFluxResult cls(mod, "ApertureFluxResult");

    cls.def("getFlag", (bool (ApertureFluxResult::*)(unsigned int) const) & ApertureFluxResult::getFlag,
            "bit"_a);
    cls.def("getFlag",
            (bool (ApertureFluxResult::*)(std::string const &name) const) & ApertureFluxResult::getFlag,
            "name"_a);
    cls.def("setFlag", &ApertureFluxResult::setFlag, "index"_a, "value"_a);
    cls.def("unsetFlag", &ApertureFluxResult::unsetFlag, "index"_a);
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "ApertureFluxTransform");

    cls.def(py::init<ApertureFluxTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    cls.def("__call__", &ApertureFluxTransform::operator(), "inputCatalog"_a, "outputCatalog"_a, "wcs"_a,
            "photoCalib"_a);

    return cls;
}

WRAP(ApertureFlux) {
    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    declareAFluxResult(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    // no need to make ApertureFluxControl::Result visible to Python
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<ApertureFluxAlgorithm, ApertureFluxControl, ApertureFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);
}

}
}

namespace alg {
namespace {

WRAP(Algorithm) {
    py::class_<BaseAlgorithm, std::shared_ptr<BaseAlgorithm>> clsBaseAlgorithm(mod, "BaseAlgorithm");
    py::class_<SingleFrameAlgorithm, std::shared_ptr<SingleFrameAlgorithm>, BaseAlgorithm>
            clsSingleFrameAlgorithm(mod, "SingleFrameAlgorithm");
    py::class_<SimpleAlgorithm, std::shared_ptr<SimpleAlgorithm>, SingleFrameAlgorithm> clsSimpleAlgorithm(
            mod, "SimpleAlgorithm", py::multiple_inheritance());

    clsBaseAlgorithm.def("fail", &BaseAlgorithm::fail, "measRecord"_a, "error"_a = NULL);
    //clsBaseAlgorithm.def("getLogName", &SimpleAlgorithm::getLogName);

    clsSingleFrameAlgorithm.def("measure", &SingleFrameAlgorithm::measure, "record"_a, "exposure"_a);

    clsSimpleAlgorithm.def("measureForced", &SimpleAlgorithm::measureForced, "measRecord"_a, "exposure"_a,
                           "refRecord"_a, "refWcs"_a);
}

}
}

namespace caflux {
namespace {

WRAP(CircularApertureFlux) {
    py::class_<CircularApertureFluxAlgorithm,
            std::shared_ptr<CircularApertureFluxAlgorithm>, ApertureFluxAlgorithm> cls(mod, "CircularApertureFluxAlgorithm");

    cls.def(py::init<CircularApertureFluxAlgorithm::Control const &, std::string const &,
                     afw::table::Schema &, daf::base::PropertySet &>(),
            "ctrl"_a, "name"_a, "schema"_a, "metadata"_a);

    cls.def("measure", &CircularApertureFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
}
}
}

namespace exceptions {
namespace {

WRAP(Exceptions) {
    using pex::exceptions::DomainError;
    using pex::exceptions::RuntimeError;

    auto clsFatalAlgorithmError =
            py::class_<FatalAlgorithmError, RuntimeError>(mod, "FatalAlgorithmError");
    py::register_exception<FatalAlgorithmError>(mod, "FatalAlgorithmException", PyExc_RuntimeError);
    auto clsMeasurementError =
            py::class_<MeasurementError, RuntimeError>(mod, "MeasurementError");
    py::register_exception<MeasurementError>(mod, "MeasurementException", PyExc_RuntimeError);

    auto clsPixelValueError =
            py::class_<PixelValueError, DomainError>(mod, "PixelValueError");
    py::register_exception<PixelValueError>(mod, "PixelValueException", PyExc_RuntimeError);

    clsMeasurementError.def(py::init<std::string const &, std::size_t>(), "message"_a, "flagBit"_a);
    clsFatalAlgorithmError.def(py::init<std::string const &>(), "message"_a);
    clsPixelValueError.def(py::init<std::string const &>(), "message"_a);

    clsMeasurementError.def("getFlagBit", &MeasurementError::getFlagBit);
}

}
}

namespace inputUtils {
namespace {

WRAP(InputUtilities) {
    py::class_<SafeCentroidExtractor> clsSafeCentroidExtractor(mod, "SafeCentroidExtractor");

    clsSafeCentroidExtractor.def(py::init<afw::table::Schema &, std::string const &, bool>(), "schema"_a,
                                 "name"_a, "isCentroider"_a = false);

    clsSafeCentroidExtractor.def("__call__",
                                 [](SafeCentroidExtractor const &self, afw::table::SourceRecord &record,
                                    FlagHandler const &flags) { return self(record, flags); },
                                 "record"_a, "flags"_a);
}

}
}

namespace pixel {
namespace {

WRAP(PixelFlags) {

    py::class_<PixelFlagsAlgorithm, std::shared_ptr<PixelFlagsAlgorithm>, SimpleAlgorithm>
            clsPixelFlagsAlgorithm(mod, "PixelFlagsAlgorithm");
    py::class_<PixelFlagsControl> clsPixelFlagsControl(mod, "PixelFlagsControl");

    clsPixelFlagsAlgorithm.def(
            py::init<PixelFlagsAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    clsPixelFlagsControl.def(py::init<>());

    clsPixelFlagsAlgorithm.def("measure", &PixelFlagsAlgorithm::measure, "measRecord"_a, "exposure"_a);
    clsPixelFlagsAlgorithm.def("fail", &PixelFlagsAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    LSST_DECLARE_CONTROL_FIELD(clsPixelFlagsControl, PixelFlagsControl, masksFpAnywhere);
    LSST_DECLARE_CONTROL_FIELD(clsPixelFlagsControl, PixelFlagsControl, masksFpCenter);
}

}
}

namespace transform {
namespace {

WRAP(Transform) {
    py::class_<BaseTransform, std::shared_ptr<BaseTransform>> cls(mod, "BaseTransform");

    cls.def("__call__", &BaseTransform::operator(), "inputCatalog"_a, "outputCatalog"_a, "wcs"_a,
            "photoCalib"_a);
}

}
}

namespace nativeCentroid {
namespace {
using PyCentroidAlgorithm =
        py::class_<NaiveCentroidAlgorithm, std::shared_ptr<NaiveCentroidAlgorithm>, SimpleAlgorithm>;
using PyCentroidControl = py::class_<NaiveCentroidControl>;
using PyCentroidTransform =
        py::class_<NaiveCentroidTransform, std::shared_ptr<NaiveCentroidTransform>, CentroidTransform>;

PyCentroidControl declareCentroidControl(py::module &mod) {
    PyCentroidControl cls(mod, "NaiveCentroidControl");

    cls.def(py::init<>());

    LSST_DECLARE_CONTROL_FIELD(cls, NaiveCentroidControl, background);
    LSST_DECLARE_CONTROL_FIELD(cls, NaiveCentroidControl, doFootprintCheck);
    LSST_DECLARE_CONTROL_FIELD(cls, NaiveCentroidControl, maxDistToPeak);

    return cls;
}

PyCentroidAlgorithm declareCentroidAlgorithm(py::module &mod) {
    PyCentroidAlgorithm cls(mod, "NaiveCentroidAlgorithm");

    cls.attr("FAILURE") = py::cast(NaiveCentroidAlgorithm::FAILURE);
    cls.attr("NO_COUNTS") = py::cast(NaiveCentroidAlgorithm::NO_COUNTS);
    cls.attr("EDGE") = py::cast(NaiveCentroidAlgorithm::EDGE);

    cls.def(py::init<NaiveCentroidAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("measure", &NaiveCentroidAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &NaiveCentroidAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyCentroidTransform declareCentroidTransform(py::module &mod) {
    PyCentroidTransform cls(mod, "NaiveCentroidTransform");

    cls.def(py::init<NaiveCentroidTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}
WRAP(NaiveCentroid) {

    auto clsCentroidControl = declareCentroidControl(mod);
    auto clsCentroidAlgorithm = declareCentroidAlgorithm(mod);
    auto clsCentroidTransform = declareCentroidTransform(mod);

    clsCentroidAlgorithm.attr("Control") = clsCentroidControl;
    clsCentroidTransform.attr("Control") = clsCentroidControl;

    python::declareAlgorithm<NaiveCentroidAlgorithm, NaiveCentroidControl, NaiveCentroidTransform>(
            clsCentroidAlgorithm, clsCentroidControl, clsCentroidTransform);
}

}
}
WRAP(Base) {
    transform::wrapTransform(mod);
    sinc::wrapSincCoeffs(mod);
    shape::wrapShapeUtilities(mod);
    inputUtils::wrapInputUtilities(mod);
    fluxUtils::wrapFluxUtilities(mod);
    flags::wrapFlagHandler(mod);
    exceptions::wrapExceptions(mod);
    alg::wrapAlgorithm(mod);
    pixel::wrapPixelFlags(mod);
    centroidUtils::wrapCentroidUtilities(mod);
    nativeCentroid::wrapNaiveCentroid(mod);
    apflux::wrapApertureFlux(mod);
    blend::wrapBlendedness(mod);
    caflux::wrapCircularApertureFlux(mod);
    gflux::wrapGaussianFlux(mod);
    background::wrapLocalBackground(mod);
    pkflux::wrapPeakLikelihoodFlux(mod);
    pflux::wrapPsfFlux(mod);
    aflux::wrapScaledApertureFlux(mod);
    sdssCentroid::wrapSdssCentroid(mod);
    sdssShape::wrapSdssShape(mod);
}
}
}
}
