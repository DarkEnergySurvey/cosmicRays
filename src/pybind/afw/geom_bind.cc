#include "pybind/afw_bind.h"

#include <cstdint>
#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>

#include "Eigen/Core"
#include "astshim.h"
#include "ndarray/eigen.h"
#include "ndarray/pybind11.h"
#include "pybind11/eigen.h"

#include "lsst/afw/fits.h"
#include "lsst/afw/geom/Endpoint.h"
#include "lsst/afw/geom/SipApproximation.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/afw/geom/Transform.h"
#include "lsst/afw/geom/detail/frameSetUtils.h"
#include "lsst/afw/geom/ellipses/BaseCore.h"
#include "lsst/afw/geom/ellipses/ConformalShear.h"
#include "lsst/afw/geom/ellipses/Convolution.h"
#include "lsst/afw/geom/ellipses/Distortion.h"
#include "lsst/afw/geom/ellipses/EllipticityBase.h"
#include "lsst/afw/geom/ellipses/PixelRegion.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/geom/ellipses/ReducedShear.h"
#include "lsst/afw/geom/ellipses/Separable.h"
#include "lsst/afw/geom/ellipses/Transformer.h"
#include "lsst/afw/geom/ellipses/radii.h"
#include "lsst/afw/geom/polygon/Polygon.h"
#include "lsst/afw/geom/transformFactory.h"
#include "lsst/afw/geom/wcsUtils.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/afw/typehandling/Storable.h"
#include "lsst/daf/base.h"
#include "lsst/geom.h"
#include "lsst/geom/AffineTransform.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/Box.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/SpherePoint.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/exceptions/python/Exception.h"
#include "lsst/utils/python.h"


namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace afw {
namespace geom {

namespace {
// Return a string consisting of "_pythonClassName_[_fromNAxes_->_toNAxes_]",
// for example "TransformGenericToPoint2[4->2]"
template <class Class>
std::string formatStr(Class const &self, std::string const &pyClassName) {
    std::ostringstream os;
    os << pyClassName;
    os << "[" << self.getFromEndpoint().getNAxes() << "->" << self.getToEndpoint().getNAxes() << "]";
    return os.str();
}

template <class FromEndpoint, class ToEndpoint, class NextToEndpoint, class PyClass>
void declareMethodTemplates(PyClass &cls) {
    using ThisTransform = Transform<FromEndpoint, ToEndpoint>;
    using NextTransform = Transform<ToEndpoint, NextToEndpoint>;
    using SeriesTransform = Transform<FromEndpoint, NextToEndpoint>;
    // Need Python-specific logic to give sensible errors for mismatched Transform types
    cls.def("_then",
            (std::shared_ptr<SeriesTransform>(ThisTransform::*)(NextTransform const &, bool) const) &
                    ThisTransform::template then<NextToEndpoint>,
            "next"_a, "simplify"_a = true);
}

// Declare Transform<FromEndpoint, ToEndpoint> using python class name Transform<X>To<Y>
// where <X> and <Y> are the prefix of the from endpoint and to endpoint class, respectively,
// for example TransformGenericToPoint2
template <class FromEndpoint, class ToEndpoint>
void declareTransform(py::module &mod) {
    using Class = Transform<FromEndpoint, ToEndpoint>;
    using ToPoint = typename ToEndpoint::Point;
    using ToArray = typename ToEndpoint::Array;
    using FromPoint = typename FromEndpoint::Point;
    using FromArray = typename FromEndpoint::Array;

    std::string const pyClassName = Class::getShortClassName();

    py::class_<Class, std::shared_ptr<Class>> cls(mod, pyClassName.c_str());

    cls.def(py::init<ast::FrameSet const &, bool>(), "frameSet"_a, "simplify"_a = true);
    cls.def(py::init<ast::Mapping const &, bool>(), "mapping"_a, "simplify"_a = true);

    cls.def_property_readonly("hasForward", &Class::hasForward);
    cls.def_property_readonly("hasInverse", &Class::hasInverse);
    cls.def_property_readonly("fromEndpoint", &Class::getFromEndpoint);
    cls.def_property_readonly("toEndpoint", &Class::getToEndpoint);

    // Return a copy of the contained Mapping in order to assure changing the returned Mapping
    // will not affect the contained Mapping (since Python ignores constness)
    cls.def("getMapping", [](Class const &self) { return self.getMapping()->copy(); });

    cls.def("applyForward", py::overload_cast<FromArray const &>(&Class::applyForward, py::const_),
            "array"_a);
    cls.def("applyForward", py::overload_cast<FromPoint const &>(&Class::applyForward, py::const_),
            "point"_a);
    cls.def("applyInverse", py::overload_cast<ToArray const &>(&Class::applyInverse, py::const_), "array"_a);
    cls.def("applyInverse", py::overload_cast<ToPoint const &>(&Class::applyInverse, py::const_), "point"_a);
    cls.def("inverted", &Class::inverted);
    /* Need some extra handling of ndarray return type in Python to prevent dimensions
     * of length 1 from being deleted */
    cls.def("_getJacobian", &Class::getJacobian);
    // Do not wrap getShortClassName because it returns the name of the class;
    // use `<class>.__name__` or `type(<instance>).__name__` instead.
    // Do not wrap readStream or writeStream because C++ streams are not easy to wrap.
    cls.def_static("readString", &Class::readString);
    cls.def("writeString", &Class::writeString);

    declareMethodTemplates<FromEndpoint, ToEndpoint, GenericEndpoint>(cls);
    declareMethodTemplates<FromEndpoint, ToEndpoint, Point2Endpoint>(cls);
    declareMethodTemplates<FromEndpoint, ToEndpoint, SpherePointEndpoint>(cls);

    // str(self) = "<Python class name>[<nIn>-><nOut>]"
    cls.def("__str__", [pyClassName](Class const &self) { return formatStr(self, pyClassName); });
    // repr(self) = "lsst.afw.geom.<Python class name>[<nIn>-><nOut>]"
    cls.def("__repr__",
            [pyClassName](Class const &self) { return "lsst.afw.geom." + formatStr(self, pyClassName); });

    table::io::python::addPersistableMethods<Class>(cls);
}

WRAP(Transform) {
    declareTransform<GenericEndpoint, GenericEndpoint>(mod);
    declareTransform<GenericEndpoint, Point2Endpoint>(mod);
    declareTransform<GenericEndpoint, SpherePointEndpoint>(mod);
    declareTransform<Point2Endpoint, GenericEndpoint>(mod);
    declareTransform<Point2Endpoint, Point2Endpoint>(mod);
    declareTransform<Point2Endpoint, SpherePointEndpoint>(mod);
    declareTransform<SpherePointEndpoint, GenericEndpoint>(mod);
    declareTransform<SpherePointEndpoint, Point2Endpoint>(mod);
    declareTransform<SpherePointEndpoint, SpherePointEndpoint>(mod);
}

WRAP(SkyWcs) {
    //py::module::import("lsst.geom");
    //py::module::import("lsst.afw.geom.transform");
    //py::module::import("lsst.afw.typehandling");

    mod.def("makeCdMatrix", makeCdMatrix, "scale"_a, "orientation"_a = 0 * lsst::geom::degrees,
            "flipX"_a = false);
    mod.def("makeFlippedWcs", makeFlippedWcs, "wcs"_a, "flipLR"_a, "flipTB"_a, "center"_a);
    mod.def("makeModifiedWcs", makeModifiedWcs, "pixelTransform"_a, "wcs"_a, "modifyActualPixels"_a);
    mod.def("makeSkyWcs",
            (std::shared_ptr<SkyWcs>(*)(lsst::geom::Point2D const &, lsst::geom::SpherePoint const &,
                                        Eigen::Matrix2d const &, std::string const &))makeSkyWcs,
            "crpix"_a, "crval"_a, "cdMatrix"_a, "projection"_a = "TAN");
    mod.def("makeSkyWcs", (std::shared_ptr<SkyWcs>(*)(daf::base::PropertySet &, bool))makeSkyWcs,
            "metadata"_a, "strip"_a = false);
    mod.def("makeSkyWcs",
            (std::shared_ptr<SkyWcs>(*)(TransformPoint2ToPoint2 const &, lsst::geom::Angle const &, bool,
                                        lsst::geom::SpherePoint const &, std::string const &))makeSkyWcs,
            "pixelsToFieldAngle"_a, "orientation"_a, "flipX"_a, "boresight"_a, "projection"_a = "TAN");
    mod.def("makeTanSipWcs",
            (std::shared_ptr<SkyWcs>(*)(lsst::geom::Point2D const &, lsst::geom::SpherePoint const &,
                                        Eigen::Matrix2d const &, Eigen::MatrixXd const &,
                                        Eigen::MatrixXd const &))makeTanSipWcs,
            "crpix"_a, "crval"_a, "cdMatrix"_a, "sipA"_a, "sipB"_a);
    mod.def("makeTanSipWcs",
            (std::shared_ptr<SkyWcs>(*)(lsst::geom::Point2D const &, lsst::geom::SpherePoint const &,
                                        Eigen::Matrix2d const &, Eigen::MatrixXd const &,
                                        Eigen::MatrixXd const &, Eigen::MatrixXd const &,
                                        Eigen::MatrixXd const &))makeTanSipWcs,
            "crpix"_a, "crval"_a, "cdMatrix"_a, "sipA"_a, "sipB"_a, "sipAp"_a, "sipBp"_a);
    mod.def("makeWcsPairTransform", makeWcsPairTransform, "src"_a, "dst"_a);
    mod.def("getIntermediateWorldCoordsToSky", getIntermediateWorldCoordsToSky, "wcs"_a, "simplify"_a = true);
    mod.def("getPixelToIntermediateWorldCoords", getPixelToIntermediateWorldCoords, "wcs"_a,
            "simplify"_a = true);

    py::class_<SkyWcs, std::shared_ptr<SkyWcs>, typehandling::Storable> cls(mod, "SkyWcs");

    cls.def(py::init<daf::base::PropertySet &, bool>(), "metadata"_a, "strip"_a = false);
    cls.def(py::init<ast::FrameDict const &>(), "frameDict"_a);

    cls.def("__eq__", &SkyWcs::operator==, py::is_operator());
    cls.def("__ne__", &SkyWcs::operator!=, py::is_operator());

    table::io::python::addPersistableMethods<SkyWcs>(cls);

    cls.def("copyAtShiftedPixelOrigin", &SkyWcs::copyAtShiftedPixelOrigin, "shift"_a);
    cls.def("getFitsMetadata", &SkyWcs::getFitsMetadata, "precise"_a = false);
    cls.def("getPixelScale",
            (lsst::geom::Angle(SkyWcs::*)(lsst::geom::Point2D const &) const) & SkyWcs::getPixelScale,
            "pixel"_a);
    cls.def("getPixelScale", (lsst::geom::Angle(SkyWcs::*)() const) & SkyWcs::getPixelScale);
    cls.def("getPixelOrigin", &SkyWcs::getPixelOrigin);
    cls.def("getSkyOrigin", &SkyWcs::getSkyOrigin);
    cls.def("getCdMatrix",
            (Eigen::Matrix2d(SkyWcs::*)(lsst::geom::Point2D const &) const) & SkyWcs::getCdMatrix, "pixel"_a);
    cls.def("getCdMatrix", (Eigen::Matrix2d(SkyWcs::*)() const) & SkyWcs::getCdMatrix);
    cls.def("getTanWcs", &SkyWcs::getTanWcs, "pixel"_a);
    cls.def("getFrameDict", [](SkyWcs const &self) { return self.getFrameDict()->copy(); });
    cls.def("getTransform", &SkyWcs::getTransform);

    cls.def_property_readonly("isFits", &SkyWcs::isFits);
    cls.def_property_readonly("isFlipped", &SkyWcs::isFlipped);
    cls.def("linearizePixelToSky",
            (lsst::geom::AffineTransform(SkyWcs::*)(lsst::geom::SpherePoint const &,
                                                    lsst::geom::AngleUnit const &) const) &
                    SkyWcs::linearizePixelToSky,
            "coord"_a, "skyUnit"_a);
    cls.def("linearizePixelToSky",
            (lsst::geom::AffineTransform(SkyWcs::*)(lsst::geom::Point2D const &,
                                                    lsst::geom::AngleUnit const &) const) &
                    SkyWcs::linearizePixelToSky,
            "coord"_a, "skyUnit"_a);
    cls.def("linearizeSkyToPixel",
            (lsst::geom::AffineTransform(SkyWcs::*)(lsst::geom::SpherePoint const &,
                                                    lsst::geom::AngleUnit const &) const) &
                    SkyWcs::linearizeSkyToPixel,
            "coord"_a, "skyUnit"_a);
    cls.def("linearizeSkyToPixel",
            (lsst::geom::AffineTransform(SkyWcs::*)(lsst::geom::Point2D const &,
                                                    lsst::geom::AngleUnit const &) const) &
                    SkyWcs::linearizeSkyToPixel,
            "coord"_a, "skyUnit"_a);
    cls.def("pixelToSky",
            (lsst::geom::SpherePoint(SkyWcs::*)(lsst::geom::Point2D const &) const) & SkyWcs::pixelToSky,
            "pixel"_a);
    cls.def("pixelToSky", (lsst::geom::SpherePoint(SkyWcs::*)(double, double) const) & SkyWcs::pixelToSky,
            "x"_a, "y"_a);
    cls.def("pixelToSky",
            (std::vector<lsst::geom::SpherePoint>(SkyWcs::*)(std::vector<lsst::geom::Point2D> const &)
                     const) &
                    SkyWcs::pixelToSky,
            "pixel"_a);
    cls.def("skyToPixel",
            (lsst::geom::Point2D(SkyWcs::*)(lsst::geom::SpherePoint const &) const) & SkyWcs::skyToPixel,
            "sky"_a);
    cls.def("skyToPixel",
            (std::vector<lsst::geom::Point2D>(SkyWcs::*)(std::vector<lsst::geom::SpherePoint> const &)
                     const) &
                    SkyWcs::skyToPixel,
            "sky"_a);
    // Do not wrap getShortClassName because it returns the name of the class;
    // use `<class>.__name__` or `type(<instance>).__name__` instead.
    // Do not wrap readStream or writeStream because C++ streams are not easy to wrap.
    cls.def_static("readString", &SkyWcs::readString);
    cls.def("writeString", &SkyWcs::writeString);

    utils::python::addOutputOp(cls, "__str__");
    // For repr, we could instead call writeString for the very long AST Frame/Mapping output.
    utils::python::addOutputOp(cls, "__repr__");
}

WRAP(WcsUtils) {
    mod.def("createTrivialWcsMetadata", createTrivialWcsMetadata, "wcsName"_a, "xy0"_a);
    mod.def("deleteBasicWcsMetadata", deleteBasicWcsMetadata, "metadata"_a, "wcsName"_a);
    mod.def("getCdMatrixFromMetadata", getCdMatrixFromMetadata, "metadata"_a);
    mod.def("getImageXY0FromMetadata", getImageXY0FromMetadata, "metadata"_a, "wcsName"_a, "strip"_a = false);
    // getSipMatrixFromMetadata requires a pure python wrapper to return a matrix when order=0
    mod.def("_getSipMatrixFromMetadata", getSipMatrixFromMetadata, "metadata"_a, "name"_a);
    mod.def("hasSipMatrix", hasSipMatrix, "metadata"_a, "name"_a);
    mod.def("makeSipMatrixMetadata", makeSipMatrixMetadata, "matrix"_a, "name"_a);
    mod.def("makeSimpleWcsMetadata", makeSimpleWcsMetadata, "crpix"_a, "crval"_a, "cdMatrix"_a,
            "projection"_a = "TAN");
    mod.def("makeTanSipMetadata",
            (std::shared_ptr<daf::base::PropertyList>(*)(
                 lsst::geom::Point2D const&, lsst::geom::SpherePoint const&, Eigen::Matrix2d const&,
                 Eigen::MatrixXd const&, Eigen::MatrixXd const&))makeTanSipMetadata,
            "crpix"_a, "crval"_a, "cdMatrix"_a, "sipA"_a, "sipB"_a);
    mod.def("makeTanSipMetadata",
            (std::shared_ptr<daf::base::PropertyList>(*)(
                 lsst::geom::Point2D const&, lsst::geom::SpherePoint const&, Eigen::Matrix2d const&,
                 Eigen::MatrixXd const&, Eigen::MatrixXd const&, Eigen::MatrixXd const&,
                 Eigen::MatrixXd const&))makeTanSipMetadata,
            "crpix"_a, "crval"_a, "cdMatrix"_a, "sipA"_a, "sipB"_a, "sipAp"_a, "sipBp"_a);
}


WRAP(TransformFactory) {
    //py::module::import("lsst.afw.geom.transform");

    mod.def("linearizeTransform",
            (lsst::geom::AffineTransform(*)(TransformPoint2ToPoint2 const &, lsst::geom::Point2D const &)) & linearizeTransform,
            "original"_a, "point"_a);
    mod.def("makeTransform",
            (std::shared_ptr<TransformPoint2ToPoint2>(*)(lsst::geom::AffineTransform const &)) & makeTransform,
            "affine"_a);
    mod.def("makeRadialTransform",
            (std::shared_ptr<TransformPoint2ToPoint2>(*)(std::vector<double> const &)) & makeRadialTransform,
            "coeffs"_a);
    mod.def("makeRadialTransform",
            (std::shared_ptr<TransformPoint2ToPoint2>(*)(std::vector<double> const &,
                                                         std::vector<double> const &)) &
                    makeRadialTransform,
            "forwardCoeffs"_a, "inverseCoeffs"_a);
    mod.def("makeIdentityTransform", &makeIdentityTransform);
}


template <typename Pixel, typename PyClass>
void declareFlattenMethod(PyClass &cls) {
    cls.def("flatten",
            (ndarray::Array<Pixel, 1, 1>(SpanSet::*)(ndarray::Array<Pixel, 2, 0> const &,
                                                     lsst::geom::Point2I const &) const) &
                    SpanSet::flatten<Pixel, 2, 0>,
            "input"_a, "xy0"_a = lsst::geom::Point2I());
    cls.def("flatten",
            (ndarray::Array<Pixel, 2, 2>(SpanSet::*)(ndarray::Array<Pixel, 3, 0> const &,
                                                     lsst::geom::Point2I const &) const) &
                    SpanSet::flatten<Pixel, 3, 0>,
            "input"_a, "xy0"_a = lsst::geom::Point2I());
    cls.def("flatten",
            (void (SpanSet::*)(ndarray::Array<Pixel, 1, 0> const &, ndarray::Array<Pixel, 2, 0> const &,
                               lsst::geom::Point2I const &) const) &
                    SpanSet::flatten<Pixel, Pixel, 2, 0, 0>,
            "output"_a, "input"_a, "xy0"_a = lsst::geom::Point2I());
    cls.def("flatten",
            (void (SpanSet::*)(ndarray::Array<Pixel, 2, 0> const &, ndarray::Array<Pixel, 3, 0> const &,
                               lsst::geom::Point2I const &) const) &
                    SpanSet::flatten<Pixel, Pixel, 3, 0, 0>,
            "output"_a, "input"_a, "xy0"_a = lsst::geom::Point2I());
}

template <typename Pixel, typename PyClass>
void declareUnflattenMethod(PyClass &cls) {
    cls.def("unflatten",
            (ndarray::Array<Pixel, 2, 2>(SpanSet::*)(ndarray::Array<Pixel, 1, 0> const &input) const) &
                    SpanSet::unflatten<Pixel, 1, 0>);
    cls.def("unflatten",
            (ndarray::Array<Pixel, 3, 3>(SpanSet::*)(ndarray::Array<Pixel, 2, 0> const &input) const) &
                    SpanSet::unflatten<Pixel, 2, 0>);
    cls.def("unflatten",
            (void (SpanSet::*)(ndarray::Array<Pixel, 2, 0> const &, ndarray::Array<Pixel, 1, 0> const &,
                               lsst::geom::Point2I const &) const) &
                    SpanSet::unflatten<Pixel, Pixel, 1, 0, 0>,
            "output"_a, "input"_a, "xy0"_a = lsst::geom::Point2I());
    cls.def("unflatten",
            (void (SpanSet::*)(ndarray::Array<Pixel, 3, 0> const &, ndarray::Array<Pixel, 2, 0> const &,
                               lsst::geom::Point2I const &) const) &
                    SpanSet::unflatten<Pixel, Pixel, 2, 0, 0>,
            "output"_a, "input"_a, "xy0"_a = lsst::geom::Point2I());
}

template <typename Pixel, typename PyClass>
void declareSetMaskMethod(PyClass &cls) {
    cls.def("setMask", (void (SpanSet::*)(image::Mask<Pixel> &, Pixel) const) & SpanSet::setMask);
}

template <typename Pixel, typename PyClass>
void declareClearMaskMethod(PyClass &cls) {
    cls.def("clearMask", (void (SpanSet::*)(image::Mask<Pixel> &, Pixel) const) & SpanSet::clearMask);
}

template <typename Pixel, typename PyClass>
void declareIntersectMethod(PyClass &cls) {
    cls.def("intersect",
            (std::shared_ptr<SpanSet>(SpanSet::*)(image::Mask<Pixel> const &, Pixel) const) &
                    SpanSet::intersect,
            "other"_a, "bitmask"_a);
    // Default to compare any bit set
    cls.def("intersect",
            [](SpanSet const &self, image::Mask<Pixel> const &mask) {
                auto tempSpanSet = SpanSet::fromMask(mask);
                return self.intersect(*tempSpanSet);
            },
            "other"_a);
}

template <typename Pixel, typename PyClass>
void declareIntersectNotMethod(PyClass &cls) {
    cls.def("intersectNot",
            (std::shared_ptr<SpanSet>(SpanSet::*)(image::Mask<Pixel> const &, Pixel) const) &
                    SpanSet::intersectNot,
            "other"_a, "bitmask"_a);
    // Default to compare any bit set
    cls.def("intersectNot",
            [](SpanSet const &self, image::Mask<Pixel> const &mask) {
                auto tempSpanSet = SpanSet::fromMask(mask);
                return self.intersectNot(*tempSpanSet);
            },
            "other"_a);
}

template <typename Pixel, typename PyClass>
void declareUnionMethod(PyClass &cls) {
    cls.def("union",
            (std::shared_ptr<SpanSet>(SpanSet::*)(image::Mask<Pixel> const &, Pixel) const) & SpanSet::union_,
            "other"_a, "bitmask"_a);
    // Default to compare any bit set
    cls.def("union",
            [](SpanSet const &self, image::Mask<Pixel> const &mask) {
                auto tempSpanSet = SpanSet::fromMask(mask);
                return self.union_(*tempSpanSet);
            },
            "other"_a);
}

template <typename ImageT, typename PyClass>
void declareCopyImage(PyClass &cls) {
    cls.def("copyImage", &SpanSet::copyImage<ImageT>);
}

template <typename ImageT, typename PyClass>
void declareCopyMaskedImage(PyClass &cls) {
    using MaskPixel = image::MaskPixel;
    using VariancePixel = image::VariancePixel;
    cls.def("copyMaskedImage", &SpanSet::copyMaskedImage<ImageT, MaskPixel, VariancePixel>);
}

template <typename ImageT, typename PyClass>
void declareSetImage(PyClass &cls) {
    cls.def("setImage",
            (void (SpanSet::*)(image::Image<ImageT> &, ImageT, lsst::geom::Box2I const &, bool) const) &
                    SpanSet::setImage,
            "image"_a, "val"_a, "region"_a = lsst::geom::Box2I(), "doClip"_a = false);
}

template <typename MaskPixel, typename PyClass>
void declarefromMask(PyClass &cls) {
    cls.def_static("fromMask", [](image::Mask<MaskPixel> mask) { return SpanSet::fromMask(mask); });
    cls.def_static("fromMask", [](image::Mask<MaskPixel> mask, MaskPixel const &bitmask) {
        return SpanSet::fromMask(mask, bitmask);
    });
}

template <typename Pixel, typename PyClass>
void declareMaskMethods(PyClass &cls) {
    declareSetMaskMethod<Pixel>(cls);
    declareClearMaskMethod<Pixel>(cls);
    declareIntersectMethod<Pixel>(cls);
    declareIntersectNotMethod<Pixel>(cls);
    declareUnionMethod<Pixel>(cls);
}

template <typename Pixel, typename PyClass>
void declareImageTypes(PyClass &cls) {
    declareFlattenMethod<Pixel>(cls);
    declareUnflattenMethod<Pixel>(cls);
    declareCopyImage<Pixel>(cls);
    declareCopyMaskedImage<Pixel>(cls);
    declareSetImage<Pixel>(cls);
}

WRAP(SpanSet) {
    using MaskPixel = image::MaskPixel;

    //py::module::import("lsst.geom");
    //py::module::import("lsst.afw.geom.span");

    py::enum_<Stencil>(mod, "Stencil")
            .value("CIRCLE", Stencil::CIRCLE)
            .value("BOX", Stencil::BOX)
            .value("MANHATTAN", Stencil::MANHATTAN);

    py::class_<SpanSet, std::shared_ptr<SpanSet> > cls(mod, "SpanSet");

    /* SpanSet Constructors */
    cls.def(py::init<>());
    cls.def(py::init<lsst::geom::Box2I>(), "box"_a);
    cls.def(py::init<std::vector<Span>, bool>(), "spans"_a, "normalize"_a = true);

    table::io::python::addPersistableMethods<SpanSet>(cls);

    /* SpanSet Methods */
    cls.def("getArea", &SpanSet::getArea);
    cls.def("getBBox", &SpanSet::getBBox);
    cls.def("isContiguous", &SpanSet::isContiguous);
    cls.def("shiftedBy", (std::shared_ptr<SpanSet>(SpanSet::*)(int, int) const) & SpanSet::shiftedBy);
    cls.def("shiftedBy",
            (std::shared_ptr<SpanSet>(SpanSet::*)(lsst::geom::Extent2I const &) const) & SpanSet::shiftedBy);
    cls.def("clippedTo", &SpanSet::clippedTo);
    cls.def("transformedBy",
            (std::shared_ptr<SpanSet>(SpanSet::*)(lsst::geom::LinearTransform const &) const) &
                    SpanSet::transformedBy);
    cls.def("transformedBy",
            (std::shared_ptr<SpanSet>(SpanSet::*)(lsst::geom::AffineTransform const &) const) &
                    SpanSet::transformedBy);
    cls.def("transformedBy", (std::shared_ptr<SpanSet>(SpanSet::*)(TransformPoint2ToPoint2 const &) const) &
                                     SpanSet::transformedBy);
    cls.def("overlaps", &SpanSet::overlaps);
    cls.def("contains", (bool (SpanSet::*)(SpanSet const &) const) & SpanSet::contains);
    cls.def("contains", (bool (SpanSet::*)(lsst::geom::Point2I const &) const) & SpanSet::contains);
    cls.def("computeCentroid", &SpanSet::computeCentroid);
    cls.def("computeShape", &SpanSet::computeShape);
    cls.def("dilated", (std::shared_ptr<SpanSet>(SpanSet::*)(int, Stencil) const) & SpanSet::dilated,
            "radius"_a, "stencil"_a = Stencil::CIRCLE);
    cls.def("dilated", (std::shared_ptr<SpanSet>(SpanSet::*)(SpanSet const &) const) & SpanSet::dilated);
    cls.def("eroded", (std::shared_ptr<SpanSet>(SpanSet::*)(int, Stencil) const) & SpanSet::eroded,
            "radius"_a, "stencil"_a = Stencil::CIRCLE);
    cls.def("eroded", (std::shared_ptr<SpanSet>(SpanSet::*)(SpanSet const &) const) & SpanSet::eroded);
    cls.def("intersect", (std::shared_ptr<SpanSet>(SpanSet::*)(SpanSet const &) const) & SpanSet::intersect);
    cls.def("intersectNot",
            (std::shared_ptr<SpanSet>(SpanSet::*)(SpanSet const &) const) & SpanSet::intersectNot);
    cls.def("union", (std::shared_ptr<SpanSet>(SpanSet::*)(SpanSet const &) const) & SpanSet::union_);
    cls.def_static("fromShape",
                   (std::shared_ptr<SpanSet>(*)(int, Stencil, lsst::geom::Point2I)) & SpanSet::fromShape,
                   "radius"_a, "stencil"_a = Stencil::CIRCLE, "offset"_a = lsst::geom::Point2I());
    cls.def_static("fromShape",
                   [](int r, Stencil s, std::pair<int, int> point) {
                       return SpanSet::fromShape(r, s, lsst::geom::Point2I(point.first, point.second));
                   },
                   "radius"_a, "stencil"_a = Stencil::CIRCLE, "offset"_a = std::pair<int, int>(0, 0));
    cls.def_static("fromShape",
                   (std::shared_ptr<SpanSet>(*)(geom::ellipses::Ellipse const &)) & SpanSet::fromShape);
    cls.def("split", &SpanSet::split);
    cls.def("findEdgePixels", &SpanSet::findEdgePixels);
    cls.def("indices", [](SpanSet const &self) -> std::pair<std::vector<int>, std::vector<int>> {
        std::vector<int> yind;
        std::vector<int> xind;
        yind.reserve(self.getArea());
        xind.reserve(self.getArea());
        for (auto const &span : self) {
            auto y = span.getY();
            for (int x = span.getX0(); x <= span.getX1(); ++x) {
                yind.push_back(y);
                xind.push_back(x);
            }
        }
        return std::make_pair(yind, xind);
    });

    /* SpanSet Operators */
    cls.def("__eq__", [](SpanSet const &self, SpanSet const &other) -> bool { return self == other; },
            py::is_operator());
    cls.def("__ne__", [](SpanSet const &self, SpanSet const &other) -> bool { return self != other; },
            py::is_operator());
    cls.def("__iter__", [](SpanSet &self) { return py::make_iterator(self.begin(), self.end()); },
            py::keep_alive<0, 1>());
    cls.def("__len__", [](SpanSet const &self) -> decltype(self.size()) { return self.size(); });
    cls.def("__contains__", [](SpanSet &self, SpanSet const &other) -> bool { return self.contains(other); });
    cls.def("__contains__",
            [](SpanSet &self, lsst::geom::Point2I &other) -> bool { return self.contains(other); });
    cls.def("__repr__", [](SpanSet const &self) -> std::string {
        std::ostringstream os;
        image::Mask<MaskPixel> tempMask(self.getBBox());
        self.setMask(tempMask, static_cast<MaskPixel>(1));
        auto array = tempMask.getArray();
        auto dims = array.getShape();
        for (std::size_t i = 0; i < dims[0]; ++i) {
            os << "[";
            for (std::size_t j = 0; j < dims[1]; ++j) {
                os << array[i][j];
                if (j != dims[1] - 1) {
                    os << ", ";
                }
            }
            os << "]" << std::endl;
        }
        return os.str();
    });
    cls.def("__str__", [](SpanSet const &self) -> std::string {
        std::ostringstream os;
        for (auto const &span : self) {
            os << span.getY() << ": " << span.getMinX() << ".." << span.getMaxX() << std::endl;
        }
        return os.str();
    });
    // Instantiate all the templates

    declareMaskMethods<MaskPixel>(cls);

    declareImageTypes<std::uint16_t>(cls);
    declareImageTypes<std::uint64_t>(cls);
    declareImageTypes<int>(cls);
    declareImageTypes<float>(cls);
    declareImageTypes<double>(cls);

    // Extra instantiation for flatten unflatten methods
    declareFlattenMethod<long>(cls);
    declareUnflattenMethod<long>(cls);

    declarefromMask<MaskPixel>(cls);
}



// A thin wrapper around SpanPixelIterator.
// Unfortunately we cannot use py::make_iterator here, as we normally
// should, because for some reason the return values then all refer
// to the last element.
class SpanIterator {
public:
    SpanIterator(const Span &s) : _it{s.begin()}, _end{s.end()} {}
    lsst::geom::Point2I next() {
        if (_it == _end) {
            throw py::stop_iteration();
        }
        return *_it++;
    }

private:
    Span::Iterator _it;
    Span::Iterator _end;
};

static void declareSpanIterator(py::module &mod) {
    py::class_<SpanIterator> cls(mod, "SpanIterator");
    cls.def("__iter__", [](SpanIterator &it) -> SpanIterator & { return it; });
    cls.def("__next__", &SpanIterator::next);
}

WRAP(Span) {
    //py::module::import("lsst.geom");

    declareSpanIterator(mod);

    py::class_<Span, std::shared_ptr<Span>> cls(mod, "Span");
    cls.def(py::init<int, int, int>());
    cls.def(py::init<>());
    cls.def("__eq__", &Span::operator==, py::is_operator());
    cls.def("__ne__", &Span::operator!=, py::is_operator());
    cls.def("__lt__", &Span::operator<, py::is_operator());
    cls.def("__len__", &Span::getWidth);
    cls.def("__str__", &Span::toString);
    // unlike most iterators, SpanPixelIterator doesn't actually refer
    // back to its container (the Span), and there's no need to keep the
    // Span alive for the lifetime of the iterator.
    cls.def("__iter__", [](const Span &s) { return SpanIterator(s); });
    cls.def("getX0", (int (Span::*)() const) & Span::getX0);
    cls.def("getX1", (int (Span::*)() const) & Span::getX1);
    cls.def("getY", (int (Span::*)() const) & Span::getY);
    cls.def("getWidth", &Span::getWidth);
    cls.def("getMinX", &Span::getMinX);
    cls.def("getMaxX", &Span::getMaxX);
    cls.def("getBeginX", &Span::getBeginX);
    cls.def("getEndX", &Span::getEndX);
    cls.def("getMin", &Span::getMin);
    cls.def("getMax", &Span::getMax);
    cls.def("contains", (bool (Span::*)(int) const) & Span::contains);
    cls.def("contains", (bool (Span::*)(int, int) const) & Span::contains);
    cls.def("contains", (bool (Span::*)(lsst::geom::Point2I const &) const) & Span::contains);
    cls.def("isEmpty", &Span::isEmpty);
    cls.def("toString", &Span::toString);
    cls.def("shift", &Span::shift);
}


WRAP(SipApproximation) {
    //py::module::import("lsst.geom");
    //py::module::import("lsst.afw.geom.transform");

    py::class_<SipApproximation, std::shared_ptr<SipApproximation>> cls(mod, "SipApproximation");

    cls.def(
        py::init<
            std::shared_ptr<TransformPoint2ToPoint2>,
            lsst::geom::Point2D const &,
            Eigen::MatrixXd const &,
            lsst::geom::Box2D const &,
            lsst::geom::Extent2I const &,
            int,
            bool,
            double
        >(),
        "pixelToIwc"_a, "crpix"_a, "cd"_a, "bbox"_a, "gridShape"_a, "order"_a,
        "useInverse"_a=true, "svdThreshold"_a=-1
    );

    cls.def(
        py::init<
            std::shared_ptr<TransformPoint2ToPoint2>,
            lsst::geom::Point2D const &,
            Eigen::MatrixXd const &,
            lsst::geom::Box2D const &,
            lsst::geom::Extent2I const &,
            ndarray::Array<double const, 2> const &,
            ndarray::Array<double const, 2> const &,
            ndarray::Array<double const, 2> const &,
            ndarray::Array<double const, 2> const &,
            bool
        >(),
        "pixelToIwc"_a, "crpix"_a, "cd"_a, "bbox"_a, "gridShape"_a,
        "a"_a, "b"_a, "ap"_a, "bp"_a, "useInverse"_a=true
    );

    using ScalarTransform = lsst::geom::Point2D (SipApproximation::*)(lsst::geom::Point2D const &) const;
    using VectorTransform = std::vector<lsst::geom::Point2D> (SipApproximation::*)(
        std::vector<lsst::geom::Point2D> const &) const;

    cls.def("getOrder", &SipApproximation::getOrder);
    cls.def("getA", py::overload_cast<int, int>(&SipApproximation::getA, py::const_), "p"_a, "q"_a);
    cls.def("getB", py::overload_cast<int, int>(&SipApproximation::getB, py::const_), "p"_a, "q"_a);
    cls.def("getAP", py::overload_cast<int, int>(&SipApproximation::getAP, py::const_), "p"_a, "q"_a);
    cls.def("getBP", py::overload_cast<int, int>(&SipApproximation::getBP, py::const_), "p"_a, "q"_a);
    cls.def("getA", py::overload_cast<>(&SipApproximation::getA, py::const_));
    cls.def("getB", py::overload_cast<>(&SipApproximation::getB, py::const_));
    cls.def("getAP", py::overload_cast<>(&SipApproximation::getAP, py::const_));
    cls.def("getBP", py::overload_cast<>(&SipApproximation::getBP, py::const_));
    cls.def("applyForward", (ScalarTransform)&SipApproximation::applyForward);
    cls.def("applyForward", (VectorTransform)&SipApproximation::applyForward);
    cls.def("applyInverse", (ScalarTransform)&SipApproximation::applyInverse);
    cls.def("applyInverse", (VectorTransform)&SipApproximation::applyInverse);
    cls.def("getGridStep", &SipApproximation::getGridStep);
    cls.def("getGridShape", &SipApproximation::getGridShape);
    cls.def("getBBox", &SipApproximation::getBBox);
    cls.def("getPixelOrigin", &SipApproximation::getPixelOrigin);
    cls.def("getCdMatrix", &SipApproximation::getCdMatrix);
    cls.def("updateGrid", &SipApproximation::updateGrid, "shape"_a);
    cls.def("refineGrid", &SipApproximation::refineGrid, "factor"_a=2);
    cls.def("fit", &SipApproximation::fit, "order"_a, "svdThreshold"_a=-1);
    cls.def("computeMaxDeviation", &SipApproximation::computeMaxDeviation);
}

/*
Add `__str__`, `__repr__` and `getClassPrefix` methods to an Endpoint pybind11 wrapper

str(self) = "GenericEndpoint(_nAxes_)" for GenericEndpoint, e.g. "GenericEndpoint(4)";
            "_typeName_()" for all other Endpoint classes, e.g. "SpherePointEndpoint()",
repr(self) = "lsst.afw.geom." + str(self), e.g. "lsst.afw.geom.GenericEndpoint(4)"
*/
template <typename PyClass>
void addStrAndRepr(PyClass& cls) {
    using Class = typename PyClass::type;  // C++ class associated with pybind11 wrapper class
    utils::python::addOutputOp(cls, "__str__");
    cls.def("__repr__", [](Class const& self) {
        std::ostringstream os;
        os << "lsst.afw.geom." << self;
        return os.str();
    });
    cls.def_static("getClassPrefix", &Class::getClassPrefix);
}

/*
Add getNPoints, dataFromPoint, dataFromArray, pointFromData and arrayFromData
*/
template <typename PyClass>
void addDataConverters(PyClass& cls) {
    using Class = typename PyClass::type;  // C++ class associated with pybind11 wrapper class
    cls.def("getNPoints", &Class::getNPoints);
    cls.def("dataFromPoint", &Class::dataFromPoint);
    cls.def("dataFromArray", &Class::dataFromArray);
    cls.def("arrayFromData", &Class::arrayFromData);
    cls.def("pointFromData", &Class::pointFromData);
}

/*
Add makeFrame method
*/
template <typename PyClass>
void addMakeFrame(PyClass& cls) {
    using Class = typename PyClass::type;  // C++ class associated with pybind11 wrapper class
    // return a deep copy so Python cannot modify the internal state
    cls.def("makeFrame", [](Class const& self) {
        auto frame = self.makeFrame();
        return frame->copy();
    });
}

// Allow Python classes to be compared across different BaseEndpoints
template <typename SelfClass, typename OtherClass, typename PyClass>
std::enable_if_t<std::is_base_of<SelfClass, OtherClass>::value> addEquals(PyClass& cls) {
    cls.def("__eq__", &SelfClass::operator==);
    cls.def("__ne__", &SelfClass::operator!=);
}

template <typename SelfClass, typename OtherClass, typename PyClass>
std::enable_if_t<!std::is_base_of<SelfClass, OtherClass>::value> addEquals(PyClass& cls) {
    cls.def("__eq__", [](SelfClass const& self, OtherClass const& other) { return false; });
    cls.def("__ne__", [](SelfClass const& self, OtherClass const& other) { return true; });
}

template <typename SelfClass, typename PyClass>
void addAllEquals(PyClass& cls) {
    addEquals<SelfClass, GenericEndpoint>(cls);
    addEquals<SelfClass, Point2Endpoint>(cls);
    addEquals<SelfClass, SpherePointEndpoint>(cls);
}

/*
 * Declare BaseVectorEndpoint<Point, Array>;
 * this is meant to be called by other `declare...` functions;
 */
template <typename Point, typename Array>
void declareBaseEndpoint(py::module& mod, std::string const& suffix) {
    using Class = BaseEndpoint<Point, Array>;
    std::string const pyClassName = "_BaseEndpoint" + suffix;

    py::class_<Class, std::shared_ptr<Class>> cls(mod, pyClassName.c_str());

    cls.def_property_readonly("nAxes", &Class::getNAxes);
    addDataConverters(cls);
    addMakeFrame(cls);
    cls.def("normalizeFrame", &Class::normalizeFrame);
    addAllEquals<Class>(cls);
}

// Declare BaseVectorEndpoint and all subclasses (the corresponding BaseEndpoint)
// This is meant to be called by other `declare...` functions;
template <typename Point>
void declareBaseVectorEndpoint(py::module& mod, std::string const& suffix) {
    using Class = BaseVectorEndpoint<Point>;
    using Array = typename Class::Array;
    std::string const pyClassName = "_BaseVectorEndpoint" + suffix;

    declareBaseEndpoint<Point, Array>(mod, suffix);

    py::class_<Class, std::shared_ptr<Class>, BaseEndpoint<Point, Array>> cls(mod, pyClassName.c_str());

    addDataConverters(cls);
}

// Declare GenericEndpoint and all subclasses
void declareGenericEndpoint(py::module& mod) {
    using Class = GenericEndpoint;
    using Point = typename Class::Point;
    using Array = typename Class::Array;

    declareBaseEndpoint<Point, Array>(mod, "Generic");

    py::class_<Class, std::shared_ptr<Class>, BaseEndpoint<Point, Array>> cls(mod, "GenericEndpoint");

    cls.def(py::init<int>(), "nAxes"_a);

    addStrAndRepr(cls);
}

/// @internal declare PointNEndpoint (for N = 2 or 3) and all subclasses
void declarePoint2Endpoint(py::module& mod) {
    using Class = Point2Endpoint;
    using Point = typename Class::Point;
    std::string const pointNumStr = "Point2";
    std::string const pyClassName = pointNumStr + "Endpoint";

    declareBaseVectorEndpoint<Point>(mod, pointNumStr);

    py::class_<Class, std::shared_ptr<Class>, BaseVectorEndpoint<Point>> cls(mod, pyClassName.c_str());

    cls.def(py::init<>());
    // do not wrap the constructor that takes nAxes; it is an implementation detail

    cls.def("normalizeFrame", &Class::normalizeFrame);
    addStrAndRepr(cls);
}

/// @internal declare SpherePointEndpoint and all subclasses
void declareSpherePointEndpoint(py::module& mod) {
    using Class = SpherePointEndpoint;
    using Point = typename Class::Point;

    declareBaseVectorEndpoint<Point>(mod, "SpherePoint");

    py::class_<Class, std::shared_ptr<Class>, BaseVectorEndpoint<Point>> cls(mod, "SpherePointEndpoint");

    cls.def(py::init<>());
    // do not wrap the constructor that takes nAxes; it is an implementation detail

    addMakeFrame(cls);
    cls.def("normalizeFrame", &Class::normalizeFrame);
    addStrAndRepr(cls);
}

WRAP(Endpoint) {
    //py::module::import("lsst.geom");

    declareGenericEndpoint(mod);
    declarePoint2Endpoint(mod);
    declareSpherePointEndpoint(mod);
}

}  // namespace


namespace polygon {
namespace {
WRAP(Polygon) {
    //py::module::import("lsst.pex.exceptions");
    //py::module::import("lsst.afw.typehandling");

    // TODO: Commented-out code is waiting until needed and is untested.
    // Add tests for it and enable it or remove it before the final pybind11 merge.

    /* Module level */
    py::class_<Polygon, std::shared_ptr<Polygon>, typehandling::Storable> clsPolygon(mod, "Polygon");

    py::register_exception<SinglePolygonException>(mod, "SinglePolygonException", PyExc_RuntimeError);

    /* Member types and enums */

    /* Constructors */
    clsPolygon.def(py::init<Polygon::Box const &>());
    clsPolygon.def(py::init<Polygon::Box const &, TransformPoint2ToPoint2 const &>());
    clsPolygon.def(py::init<Polygon::Box const &, lsst::geom::AffineTransform const &>());
    clsPolygon.def(py::init<std::vector<Polygon::Point> const &>());

    table::io::python::addPersistableMethods<Polygon>(clsPolygon);

    /* Operators */
    clsPolygon.def("__eq__", [](Polygon const &self, Polygon const &other) { return self == other; },
                   py::is_operator());
    clsPolygon.def("__ne__", [](Polygon const &self, Polygon const &other) { return self != other; },
                   py::is_operator());

    /* Members */
    clsPolygon.def("getNumEdges", &Polygon::getNumEdges);
    clsPolygon.def("getBBox", &Polygon::getBBox);
    clsPolygon.def("calculateCenter", &Polygon::calculateCenter);
    clsPolygon.def("calculateArea", &Polygon::calculateArea);
    clsPolygon.def("calculatePerimeter", &Polygon::calculatePerimeter);
    clsPolygon.def("getVertices", &Polygon::getVertices);
    clsPolygon.def("getEdges", &Polygon::getEdges);
    clsPolygon.def("contains", &Polygon::contains);
    clsPolygon.def("overlaps", (bool (Polygon::*)(Polygon const &) const) & Polygon::overlaps);
    clsPolygon.def("overlaps", (bool (Polygon::*)(Polygon::Box const &) const) & Polygon::overlaps);
    clsPolygon.def("intersectionSingle", (std::shared_ptr<Polygon>(Polygon::*)(Polygon const &) const) &
                                                 Polygon::intersectionSingle);
    clsPolygon.def("intersectionSingle", (std::shared_ptr<Polygon>(Polygon::*)(Polygon::Box const &) const) &
                                                 Polygon::intersectionSingle);
    clsPolygon.def("intersection",
                   (std::vector<std::shared_ptr<Polygon>>(Polygon::*)(Polygon const &) const) &
                           Polygon::intersection);
    clsPolygon.def("intersection",
                   (std::vector<std::shared_ptr<Polygon>>(Polygon::*)(Polygon::Box const &) const) &
                           Polygon::intersection);
    clsPolygon.def("unionSingle",
                   (std::shared_ptr<Polygon>(Polygon::*)(Polygon const &) const) & Polygon::unionSingle);
    clsPolygon.def("unionSingle",
                   (std::shared_ptr<Polygon>(Polygon::*)(Polygon::Box const &) const) & Polygon::unionSingle);

    // Wrap Polygon::union_ (C++) as Polygon.union (Python)
    clsPolygon.def("union", (std::vector<std::shared_ptr<Polygon>>(Polygon::*)(Polygon const &) const) &
                                    Polygon::union_);
    clsPolygon.def("union", (std::vector<std::shared_ptr<Polygon>>(Polygon::*)(Polygon::Box const &) const) &
                                    Polygon::union_);
    clsPolygon.def("symDifference",
                   (std::vector<std::shared_ptr<Polygon>>(Polygon::*)(Polygon const &) const) &
                           Polygon::symDifference);
    clsPolygon.def("symDifference",
                   (std::vector<std::shared_ptr<Polygon>>(Polygon::*)(Polygon::Box const &) const) &
                           Polygon::symDifference);
    // clsPolygon.def("simplify", &Polygon::simplify);
    clsPolygon.def("convexHull", &Polygon::convexHull);
    clsPolygon.def("transform",
                   (std::shared_ptr<Polygon>(Polygon::*)(TransformPoint2ToPoint2 const &) const) &
                           Polygon::transform);
    clsPolygon.def("transform",
                   (std::shared_ptr<Polygon>(Polygon::*)(lsst::geom::AffineTransform const &) const) &
                           Polygon::transform);
    clsPolygon.def("subSample", (std::shared_ptr<Polygon>(Polygon::*)(size_t) const) & Polygon::subSample);
    clsPolygon.def("subSample", (std::shared_ptr<Polygon>(Polygon::*)(double) const) & Polygon::subSample);
    clsPolygon.def("createImage",
                   (std::shared_ptr<afw::image::Image<float>>(Polygon::*)(lsst::geom::Box2I const &) const) &
                           Polygon::createImage);
    clsPolygon.def(
            "createImage",
            (std::shared_ptr<afw::image::Image<float>>(Polygon::*)(lsst::geom::Extent2I const &) const) &
                    Polygon::createImage);
    // clsPolygon.def("isPersistable", &Polygon::isPersistable);
}
} //namespace
} // polygon

namespace ellipses {
namespace {
template <typename Ellipticity_, typename Radius_>
void declareSeparable(py::module &mod, const std::string &suffix) {
    using Class = Separable<Ellipticity_, Radius_>;

    py::class_<Class, std::shared_ptr<Class>, BaseCore> cls(mod, ("Separable" + suffix).c_str());

    //    py::enum_<typename Class::ParameterEnum>(mod, "ParameterEnum")
    //        .value("E0", Class::ParameterEnum::E0)
    //        .value("E1", Class::ParameterEnum::E1)
    //        .value("RADIUS", Class::ParameterEnum::RADIUS)
    //        .export_values();

    cls.def(py::init<double, double, double, bool>(), "e1"_a = 0.0, "e2"_a = 0.0, "radius"_a = Radius_(),
            "normalize"_a = true);
    //    cls.def(py::init<std::complex<double> const &, double, bool>(),
    //            "complex"_a, "radius"_a=Radius_(), "normalize"_a=true);
    //    cls.def(py::init<Ellipticity_ const &, double, bool>(),
    //            "ellipticity"_a, "radius"_a=Radius_(), "normalize"_a=true);
    //    cls.def(py::init<BaseCore::ParameterVector const &, bool>(),
    //            "vector"_a, "normalize"_a=true);
    cls.def(py::init<Class const &>());
    cls.def(py::init<BaseCore const &>());

    cls.def("getE1", &Class::getE1);
    cls.def("setE1", &Class::setE1);
    cls.def("getE2", &Class::getE2);
    cls.def("setE2", &Class::setE2);
    cls.def("getRadius", (Radius_ const &(Class::*)() const) & Class::getRadius);
    cls.def("setRadius", (void (Class::*)(double)) & Class::setRadius);
    cls.def("setRadius", (void (Class::*)(Radius_ const &)) & Class::setRadius);
    cls.def("getEllipticity", (Ellipticity_ const &(Class::*)() const) & Class::getEllipticity);
    cls.def("clone", &Class::clone);
    cls.def("getName", &Class::getName);
    cls.def("normalize", &Class::normalize);
    cls.def("assign", [](Class &self, Class &other) { self = other; });
    cls.def("assign", [](Class &self, BaseCore &other) { self = other; });
    cls.def("transform", [](Class &self, lsst::geom::LinearTransform const &t) {
        return std::static_pointer_cast<Class>(self.transform(t).copy());
    });
    cls.def("transformInPlace",
            [](Class &self, lsst::geom::LinearTransform const &t) { self.transform(t).inPlace(); });
    cls.def("__str__",
            [](Class &self) { return py::str("({}, {})").format(self.getEllipticity(), self.getRadius()); });
    cls.def("__repr__", [](Class &self) {
        return py::str("Separable({}, {})").format(self.getEllipticity(), self.getRadius());
    });
}

WRAP(Separable) {
    declareSeparable<Distortion, DeterminantRadius>(mod, "DistortionDeterminantRadius");
    declareSeparable<Distortion, TraceRadius>(mod, "DistortionTraceRadius");
    declareSeparable<Distortion, LogDeterminantRadius>(mod, "DistortionLogDeterminantRadius");
    declareSeparable<Distortion, LogTraceRadius>(mod, "DistortionLogTraceRadius");

    declareSeparable<ConformalShear, DeterminantRadius>(mod, "ConformalShearDeterminantRadius");
    declareSeparable<ConformalShear, TraceRadius>(mod, "ConformalShearTraceRadius");
    declareSeparable<ConformalShear, LogDeterminantRadius>(mod, "ConformalShearLogDeterminantRadius");
    declareSeparable<ConformalShear, LogTraceRadius>(mod, "ConformalShearLogTraceRadius");

    declareSeparable<ReducedShear, DeterminantRadius>(mod, "ReducedShearDeterminantRadius");
    declareSeparable<ReducedShear, TraceRadius>(mod, "ReducedShearTraceRadius");
    declareSeparable<ReducedShear, LogDeterminantRadius>(mod, "ReducedShearLogDeterminantRadius");
    declareSeparable<ReducedShear, LogTraceRadius>(mod, "ReducedShearLogTraceRadius");
}

WRAP(ReducedShear) {
    py::class_<ReducedShear, detail::EllipticityBase> cls(mod, "ReducedShear");

    /* Constructors */
    cls.def(py::init<std::complex<double> const&>());
    cls.def(py::init<double, double>(), "e1"_a = 0.0, "e2"_a = 0.0);

    /* Members */
    //    cls.def("dAssign", (Jacobian (ReducedShear::*)(Distortion const &)) &ReducedShear::dAssign);
    //    cls.def("dAssign", (Jacobian (ReducedShear::*)(ReducedShear const &)) &ReducedShear::dAssign);
    cls.def("getAxisRatio", &ReducedShear::getAxisRatio);
    cls.def("normalize", &ReducedShear::normalize);
    cls.def("getName", &ReducedShear::getName);
    cls.def("__repr__", [](ReducedShear const& self) {
        return py::str("{}({}, {})").format(self.getName(), self.getE1(), self.getE2());
    });
}

WRAP(Radii) {
    py::class_<DeterminantRadius> clsDeterminantRadius(mod, "DeterminantRadius");

    clsDeterminantRadius.def(py::init<double>(), "value"_a = 1.0);
    clsDeterminantRadius.def("normalize", &DeterminantRadius::normalize);
    clsDeterminantRadius.def_static("getName", DeterminantRadius::getName);
    clsDeterminantRadius.def("__str__", [](DeterminantRadius const& self) { return std::to_string(self); });
    clsDeterminantRadius.def("__repr__", [](DeterminantRadius const& self) {
        return self.getName() + "(" + std::to_string(self) + ")";
    });

    py::class_<TraceRadius> clsTraceRadius(mod, "TraceRadius");

    clsTraceRadius.def(py::init<double>(), "value"_a = 1.0);
    clsTraceRadius.def("normalize", &TraceRadius::normalize);
    clsTraceRadius.def_static("getName", TraceRadius::getName);
    clsTraceRadius.def("__str__", [](TraceRadius const& self) { return std::to_string(self); });
    clsTraceRadius.def("__repr__", [](TraceRadius const& self) {
        return self.getName() + "(" + std::to_string(self) + ")";
    });

    py::class_<LogDeterminantRadius> clsLogDeterminantRadius(mod, "LogDeterminantRadius");

    clsLogDeterminantRadius.def(py::init<double>(), "value"_a = 0.0);
    clsLogDeterminantRadius.def("normalize", &LogDeterminantRadius::normalize);
    clsLogDeterminantRadius.def_static("getName", LogDeterminantRadius::getName);
    clsLogDeterminantRadius.def("__str__",
                                [](LogDeterminantRadius const& self) { return std::to_string(self); });
    clsLogDeterminantRadius.def("__repr__", [](LogDeterminantRadius const& self) {
        return self.getName() + "(" + std::to_string(self) + ")";
    });

    py::class_<LogTraceRadius> clsLogTraceRadius(mod, "LogTraceRadius");

    clsLogTraceRadius.def(py::init<double>(), "value"_a = 0.0);
    clsLogTraceRadius.def("normalize", &LogTraceRadius::normalize);
    clsLogTraceRadius.def_static("getName", LogTraceRadius::getName);
    clsLogTraceRadius.def("__str__", [](LogTraceRadius const& self) { return std::to_string(self); });
    clsLogTraceRadius.def("__repr__", [](LogTraceRadius const& self) {
        return self.getName() + "(" + std::to_string(self) + ")";
    });
}

WRAP(PixelRegion) {
    py::class_<PixelRegion> clsPixelRegion(mod, "PixelRegion");

    /* Constructors */
    clsPixelRegion.def(py::init<Ellipse const &>());

    /* Members */
    clsPixelRegion.def("getBBox", &PixelRegion::getBBox, py::return_value_policy::copy);
    clsPixelRegion.def("getSpanAt", &PixelRegion::getSpanAt);
    clsPixelRegion.def("__iter__",
                       [](const PixelRegion &self) { return py::make_iterator(self.begin(), self.end()); },
                       py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);
}

WRAP(EllipticityBase) {
    py::class_<detail::EllipticityBase> cls(mod, "EllipticityBase");

    /* Member types and enums */
    py::enum_<detail::EllipticityBase::ParameterEnum>(cls, "ParameterEnum")
            .value("E1", detail::EllipticityBase::ParameterEnum::E1)
            .value("E2", detail::EllipticityBase::ParameterEnum::E2)
            .export_values();

    /* Members */
    cls.def("getComplex",
            (std::complex<double> & (detail::EllipticityBase::*)()) & detail::EllipticityBase::getComplex);
    cls.def("setComplex", &detail::EllipticityBase::setComplex);
    cls.def("getE1", &detail::EllipticityBase::getE1);
    cls.def("setE1", &detail::EllipticityBase::setE1);
    cls.def("getE2", &detail::EllipticityBase::getE2);
    cls.def("setE2", &detail::EllipticityBase::setE2);
    cls.def("getTheta", &detail::EllipticityBase::getTheta);
    cls.def("__str__", [](detail::EllipticityBase const& self) {
        return py::str("({}, {})").format(self.getE1(), self.getE2());
    });
}

WRAP(Distortion) {
    py::class_<Distortion, detail::EllipticityBase> cls(mod, "Distortion");

    /* Constructors */
    cls.def(py::init<std::complex<double> const&>());
    cls.def(py::init<double, double>(), "e1"_a = 0.0, "e2"_a = 0.0);

    /* Members */
    //    cls.def("dAssign", (Jacobian (Distortion::*)(Distortion const &)) &Distortion::dAssign);
    //    cls.def("dAssign", (Jacobian (Distortion::*)(ReducedShear const &)) &Distortion::dAssign);
    cls.def("getAxisRatio", &Distortion::getAxisRatio);
    cls.def("normalize", &Distortion::normalize);
    cls.def("getName", &Distortion::getName);
    cls.def("__repr__", [](Distortion const& self) {
        return py::str("{}({}, {})").format(self.getName(), self.getE1(), self.getE2());
    });
}

WRAP(ConformalShear) {
    py::class_<ConformalShear, detail::EllipticityBase> cls(mod, "ConformalShear");

    /* Constructors */
    cls.def(py::init<std::complex<double> const&>());
    cls.def(py::init<double, double>(), "e1"_a = 0.0, "e2"_a = 0.0);

    /* Members */
    //    cls.def("dAssign", (Jacobian (ConformalShear::*)(Distortion const &)) &ConformalShear::dAssign);
    //    cls.def("dAssign", (Jacobian (ConformalShear::*)(ReducedShear const &)) &ConformalShear::dAssign);
    cls.def("getAxisRatio", &ConformalShear::getAxisRatio);
    cls.def("normalize", &ConformalShear::normalize);
    cls.def("getName", &ConformalShear::getName);
    cls.def("__repr__", [](ConformalShear const& self) {
        return py::str("{}({}, {})").format(self.getName(), self.getE1(), self.getE2());
    });
}

WRAP(BaseCore) {
    /* Module level */
    py::class_<BaseCore, std::shared_ptr<BaseCore>> clsBaseCore(mod, "BaseCore");

    /* Member types and enums */
    py::class_<BaseCore::Convolution> clsBaseCoreConvolution(clsBaseCore, "Convolution");
    py::class_<BaseCore::Transformer> clsBaseCoreTransformer(clsBaseCore, "Transformer");

    //    clsBaseCoreTransformer.def(py::init<BaseCore &, lsst::geom::LinearTransform const &>());
    //
    //    clsBaseCoreTransformer.def("inPlace", &BaseCore::Transformer::inPlace);
    //    clsBaseCoreTransformer.def("apply", &BaseCore::Transformer::apply);
    //    clsBaseCoreTransformer.def("d", &BaseCore::Transformer::d);
    //    clsBaseCoreTransformer.def("dTransform", &BaseCore::Transformer::dTransform);
    //
    //    clsBaseCoreTransformer.def_readwrite("input", &BaseCore::Transformer::input);
    //    clsBaseCoreTransformer.def_readonly("transform", &BaseCore::Transformer::transform);

    /* Constructors */

    /* Operators */
    clsBaseCore.def("__eq__", &BaseCore::operator==, py::is_operator());
    clsBaseCore.def("__nq__", &BaseCore::operator!=, py::is_operator());

    /* Members */
    clsBaseCore.def("getName", &BaseCore::getName);
    clsBaseCore.def("clone", &BaseCore::clone);
    clsBaseCore.def("normalize", &BaseCore::normalize);
    clsBaseCore.def("grow", &BaseCore::grow);
    clsBaseCore.def("scale", &BaseCore::scale);
    clsBaseCore.def("getArea", &BaseCore::getArea);
    clsBaseCore.def("getDeterminantRadius", &BaseCore::getDeterminantRadius);
    clsBaseCore.def("getTraceRadius", &BaseCore::getTraceRadius);
    //    clsBaseCore.def("transform", (typename BaseCore::Transformer const
    //    (BaseCore::*)(lsst::geom::LinearTransform const &) const) &BaseCore::transform);
    //    clsBaseCore.def("getGridTransform", &BaseCore::getGridTransform);
    clsBaseCore.def("convolve", (BaseCore::Convolution (BaseCore::*)(BaseCore const &)) & BaseCore::convolve);
    clsBaseCore.def("computeDimensions", &BaseCore::computeDimensions);
    clsBaseCore.def("getParameterVector", &BaseCore::getParameterVector);
    clsBaseCore.def("setParameterVector", &BaseCore::setParameterVector);
    //    clsBaseCore.def("transformInPlace", [](BaseCore & self, lsst::geom::LinearTransform const & t)
    //    {
    //       self.transform(t).inPlace();
    //    });
}

WRAP(Quadrupole) {
    /* Module level */
    py::class_<Quadrupole, std::shared_ptr<Quadrupole>, BaseCore> clsQuadrupole(mod, "Quadrupole");

    /* Member types and enums */
    typedef Eigen::Matrix<double, 2, 2, Eigen::DontAlign> Matrix;

    /* Constructors */
    clsQuadrupole.def(py::init<double, double, double, bool>(), "ixx"_a = 1.0, "iyy"_a = 1.0, "ixy"_a = 0.0,
                      "normalize"_a = false);
    clsQuadrupole.def(py::init<BaseCore::ParameterVector const &, bool>(), "vector"_a, "normalize"_a = false);
    clsQuadrupole.def(py::init<Matrix const &, bool>(), "matrix"_a, "normalize"_a = true);
    clsQuadrupole.def(py::init<Quadrupole const &>());
    clsQuadrupole.def(py::init<BaseCore const &>());
    clsQuadrupole.def(py::init<BaseCore::Convolution const &>());

    py::implicitly_convertible<BaseCore::Convolution, Quadrupole>();

    /* Operators */
    clsQuadrupole.def("__eq__", [](Quadrupole &self, Quadrupole &other) { return self == other; },
    py::is_operator());
    clsQuadrupole.def("__ne__", [](Quadrupole &self, Quadrupole &other) { return self != other; },
    py::is_operator());

    /* Members */
    clsQuadrupole.def("getIxx", &Quadrupole::getIxx);
    clsQuadrupole.def("getIyy", &Quadrupole::getIyy);
    clsQuadrupole.def("getIxy", &Quadrupole::getIxy);
    clsQuadrupole.def("setIxx", &Quadrupole::setIxx);
    clsQuadrupole.def("setIyy", &Quadrupole::setIyy);
    clsQuadrupole.def("setIxy", &Quadrupole::setIxy);
    clsQuadrupole.def("assign", [](Quadrupole &self, Quadrupole &other) { self = other; });
    clsQuadrupole.def("assign", [](Quadrupole &self, BaseCore &other) { self = other; });
    clsQuadrupole.def("transform", [](Quadrupole &self, lsst::geom::LinearTransform const &t) {
        return std::static_pointer_cast<Quadrupole>(self.transform(t).copy());
    });
    clsQuadrupole.def("transformInPlace", [](Quadrupole &self, lsst::geom::LinearTransform const &t) {
        self.transform(t).inPlace();
    });
    // TODO: pybind11 based on swig wrapper for now. Hopefully can be removed once pybind11 gets smarter
    // handling of implicit conversions
    clsQuadrupole.def("convolve", [](Quadrupole &self, BaseCore const &other) {
        return Quadrupole(self.convolve(other));
    });
}

} //namespace
} //ellipses

namespace detail {
namespace {

WRAP(FrameSetUtils) {
    //py::module::import("lsst.daf.base");
    //py::module::import("lsst.geom");

    mod.def("readFitsWcs", readFitsWcs, "metadata"_a, "strip"_a = true);
    mod.def("readLsstSkyWcs", readLsstSkyWcs, "metadata"_a, "strip"_a = true);
    mod.def("getPropertyListFromFitsChan", getPropertyListFromFitsChan, "fitsChan"_a);
}

}  // namespace
}  // namespace detail



void wrapGeom(py::module_ &mod) {
    wrapWcsUtils(mod);
    wrapEndpoint(mod);
    wrapTransform(mod);
    wrapSkyWcs(mod);
    polygon::wrapPolygon(mod);
    ellipses::wrapBaseCore(mod);
    ellipses::wrapRadii(mod);
    ellipses::wrapPixelRegion(mod);
    ellipses::wrapEllipticityBase(mod);
    ellipses::wrapDistortion(mod);
    ellipses::wrapReducedShear(mod);
    ellipses::wrapConformalShear(mod);
    ellipses::wrapSeparable(mod);
    ellipses::wrapQuadrupole(mod);
    detail::wrapFrameSetUtils(mod);
    wrapTransformFactory(mod);
    wrapSpan(mod);
    wrapSpanSet(mod);
    wrapSipApproximation(mod);
}

}
}
}
