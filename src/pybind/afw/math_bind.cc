#include "pybind/afw_bind.h"

#include <memory>
#include <vector>

#include "ndarray/pybind11.h"

#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/math/Approximate.h"
#include "lsst/afw/math/Background.h"
#include "lsst/afw/math/BoundedField.h"
#include "lsst/afw/math/ChebyshevBoundedField.h"
#include "lsst/afw/math/ConvolveImage.h"
#include "lsst/afw/math/Function.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/math/GaussianProcess.h"
#include "lsst/afw/math/Interpolate.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/afw/math/PixelAreaBoundedField.h"
#include "lsst/afw/math/ProductBoundedField.h"
#include "lsst/afw/math/Random.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/afw/math/Stack.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/math/TransformBoundedField.h"
#include "lsst/afw/math/detail/Convolve.h"
#include "lsst/afw/math/detail/Spline.h"
#include "lsst/afw/math/minimize.h"
#include "lsst/afw/math/offsetImage.h"
#include "lsst/afw/math/warpExposure.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/geom/Box.h"
#include "lsst/geom/Point.h"
#include "lsst/pex/config/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_DECLARE_HOLDER_TYPE(MyType, std::shared_ptr<MyType>);

namespace lsst {
namespace afw {
namespace math {
using ChebyClsField = py::class_<ChebyshevBoundedField, std::shared_ptr<ChebyshevBoundedField>,
                            BoundedField>;

namespace {

using CandidateList = std::vector<std::shared_ptr<SpatialCellCandidate>>;
using ClsField = py::class_<TransformBoundedField, std::shared_ptr<TransformBoundedField>, BoundedField>;

// Wrap SpatialCellCandidate (an abstract class so no constructor is wrapped)
WRAP(SpatialCellCandidate) {
    py::class_<SpatialCellCandidate, std::shared_ptr<SpatialCellCandidate>> cls(mod, "SpatialCellCandidate");

    py::enum_<SpatialCellCandidate::Status>(cls, "Status")
            .value("BAD", SpatialCellCandidate::Status::BAD)
            .value("GOOD", SpatialCellCandidate::Status::GOOD)
            .value("UNKNOWN", SpatialCellCandidate::Status::UNKNOWN)
            .export_values();

    cls.def("getXCenter", &SpatialCellCandidate::getXCenter);
    cls.def("getYCenter", &SpatialCellCandidate::getYCenter);
    cls.def("instantiate", &SpatialCellCandidate::instantiate);
    cls.def("getCandidateRating", &SpatialCellCandidate::getCandidateRating);
    cls.def("setCandidateRating", &SpatialCellCandidate::setCandidateRating);
    cls.def("getId", &SpatialCellCandidate::getId);
    cls.def("getStatus", &SpatialCellCandidate::getStatus);
    cls.def("setStatus", &SpatialCellCandidate::setStatus);
    cls.def("isBad", &SpatialCellCandidate::isBad);
}

// Wrap SpatialCellCandidateIterator
WRAP(SpatialCellCandidateIterator) {
    py::class_<SpatialCellCandidateIterator> cls(mod, "SpatialCellCandidateIterator");
    cls.def("__incr__", &SpatialCellCandidateIterator::operator++, py::is_operator());
    cls.def("__deref__",
            [](SpatialCellCandidateIterator &it) -> std::shared_ptr<SpatialCellCandidate> { return *it; },
            py::is_operator());
    cls.def("__eq__", &SpatialCellCandidateIterator::operator==, py::is_operator());
    cls.def("__ne__", &SpatialCellCandidateIterator::operator!=, py::is_operator());
    cls.def("__sub__", &SpatialCellCandidateIterator::operator-, py::is_operator());
}

// Wrap SpatialCell
WRAP(SpatialCell) {
    py::class_<SpatialCell, std::shared_ptr<SpatialCell>> cls(mod, "SpatialCell");

    cls.def(py::init<std::string const &, lsst::geom::Box2I const &, CandidateList const &>(), "label"_a,
            "bbox"_a = lsst::geom::Box2I(), "candidateList"_a = CandidateList());

    cls.def("empty", &SpatialCell::empty);
    cls.def("size", &SpatialCell::size);
    cls.def("__len__", &SpatialCell::size);
    cls.def("getLabel", &SpatialCell::getLabel);
    cls.def("begin", (SpatialCellCandidateIterator (SpatialCell::*)()) & SpatialCell::begin);
    cls.def("begin", (SpatialCellCandidateIterator (SpatialCell::*)(bool)) & SpatialCell::begin);
    cls.def("end", (SpatialCellCandidateIterator (SpatialCell::*)()) & SpatialCell::end);
    cls.def("end", (SpatialCellCandidateIterator (SpatialCell::*)(bool)) & SpatialCell::end);
    cls.def("insertCandidate", &SpatialCell::insertCandidate);
    cls.def("removeCandidate", &SpatialCell::removeCandidate);
    cls.def("setIgnoreBad", &SpatialCell::setIgnoreBad, "ignoreBad"_a);
    cls.def("getIgnoreBad", &SpatialCell::getIgnoreBad);
    cls.def("getCandidateById", &SpatialCell::getCandidateById, "id"_a, "noThrow"_a = false);
    cls.def("getLabel", &SpatialCell::getLabel);
    cls.def("getBBox", &SpatialCell::getBBox);
    cls.def("sortCandidates", &SpatialCell::sortCandidates);
    cls.def("visitCandidates",
            (void (SpatialCell::*)(CandidateVisitor *, int const, bool const, bool const)) &
                    SpatialCell::visitCandidates,
            "visitor"_a, "nMaxPerCell"_a = -1, "ignoreExceptions"_a = false, "reset"_a = true);
    cls.def("visitAllCandidates", (void (SpatialCell::*)(CandidateVisitor *, bool const, bool const)) &
                                          SpatialCell::visitAllCandidates,
            "visitor"_a, "ignoreExceptions"_a = false, "reset"_a = true);
}

// Wrap SpatialCellSet
WRAP(SpatialCellSet) {
    py::class_<SpatialCellSet, std::shared_ptr<SpatialCellSet>> cls(mod, "SpatialCellSet");

    cls.def(py::init<lsst::geom::Box2I const &, int, int>(), "region"_a, "xSize"_a, "ySize"_a = 0);

    cls.def("getCellList", &SpatialCellSet::getCellList);
    cls.def("getBBox", &SpatialCellSet::getBBox);
    cls.def("insertCandidate", &SpatialCellSet::insertCandidate);
    cls.def("sortCandidates", &SpatialCellSet::sortCandidates);
    cls.def("visitCandidates", (void (SpatialCellSet::*)(CandidateVisitor *, int const, bool const)) &
                                       SpatialCellSet::visitCandidates,
            "visitor"_a, "nMaxPerCell"_a = -1, "ignoreExceptions"_a = false);
    cls.def("visitAllCandidates",
            (void (SpatialCellSet::*)(CandidateVisitor *, bool const)) & SpatialCellSet::visitAllCandidates,
            "visitor"_a, "ignoreExceptions"_a = false);
    cls.def("getCandidateById", &SpatialCellSet::getCandidateById, "id"_a, "noThrow"_a = false);
    cls.def("setIgnoreBad", &SpatialCellSet::setIgnoreBad, "ignoreBad"_a);
}

// Wrap CandidateVisitor
WRAP(CandidateVisitor) {
    py::class_<CandidateVisitor, std::shared_ptr<CandidateVisitor>> cls(mod, "CandidateVisitor");

    cls.def(py::init<>());

    cls.def("reset", &CandidateVisitor::reset);
    cls.def("processCandidate", &CandidateVisitor::processCandidate);
}

// Wrap class SpatialCellImageCandidate (an abstract class, so no constructor is wrapped)
WRAP(SpatialCellImageCandidate) {
    py::class_<SpatialCellImageCandidate, std::shared_ptr<SpatialCellImageCandidate>, SpatialCellCandidate>
            cls(mod, "SpatialCellImageCandidate");

    cls.def_static("setWidth", &SpatialCellImageCandidate::setWidth, "width"_a);
    cls.def_static("getWidth", &SpatialCellImageCandidate::getWidth);
    cls.def_static("setHeight", &SpatialCellImageCandidate::setHeight, "height"_a);
    cls.def_static("getHeight", &SpatialCellImageCandidate::getHeight);
    cls.def("setChi2", &SpatialCellImageCandidate::setChi2, "chi2"_a);
    cls.def("getChi2", &SpatialCellImageCandidate::getChi2);
}

template <typename OutImageT, typename InImageT, typename KernelT>
void declareConvolve(py::module &mod) {
    mod.def("convolve", (void (*)(OutImageT &, InImageT const &, KernelT const &,
                                  ConvolutionControl const &))convolve<OutImageT, InImageT, KernelT>,
            "convolvedImage"_a, "inImage"_a, "kernel"_a, "convolutionControl"_a = ConvolutionControl());
    mod.def("convolve", (void (*)(OutImageT &, InImageT const &, KernelT const &, bool,
                                  bool))convolve<OutImageT, InImageT, KernelT>,
            "convolvedImage"_a, "inImage"_a, "kernel"_a, "doNormalize"_a, "doCopyEdge"_a = false);
}

template <typename ImageType1, typename ImageType2>
void declareScaledPlus(py::module &mod) {
    mod.def("scaledPlus",
            (void (*)(ImageType1 &, double, ImageType2 const &, double, ImageType2 const &))scaledPlus);
}

template <typename ImageType1, typename ImageType2>
void declareByType(py::module &mod) {
    declareConvolve<ImageType1, ImageType2, AnalyticKernel>(mod);
    declareConvolve<ImageType1, ImageType2, DeltaFunctionKernel>(mod);
    declareConvolve<ImageType1, ImageType2, FixedKernel>(mod);
    declareConvolve<ImageType1, ImageType2, LinearCombinationKernel>(mod);
    declareConvolve<ImageType1, ImageType2, SeparableKernel>(mod);
    declareConvolve<ImageType1, ImageType2, Kernel>(mod);
    declareScaledPlus<ImageType1, ImageType2>(mod);
}

template <typename PixelType1, typename PixelType2>
void declareAll(py::module &mod) {
    using M1 = image::MaskedImage<PixelType1, image::MaskPixel, image::VariancePixel>;
    using M2 = image::MaskedImage<PixelType2, image::MaskPixel, image::VariancePixel>;

    declareByType<image::Image<PixelType1>, image::Image<PixelType2>>(mod);
    declareByType<M1, M2>(mod);
}

template <typename ImageT>
void declareMakeBackground(py::module &mod) {
    mod.def("makeBackground", makeBackground<ImageT>, "img"_a, "bgCtrl"_a);
}

template <typename PixelT, typename PyClass>
void declareGetImage(PyClass &cls, std::string const &suffix) {
    cls.def(("getImage" + suffix).c_str(),
            (std::shared_ptr<lsst::afw::image::Image<PixelT>> (Background::*)(
                    Interpolate::Style const, UndersampleStyle const) const) &
                    Background::getImage<PixelT>,
            "interpStyle"_a, "undersampleStyle"_a = THROW_EXCEPTION);
    cls.def(("getImage" + suffix).c_str(),
            (std::shared_ptr<lsst::afw::image::Image<PixelT>> (Background::*)(std::string const &,
                                                                                std::string const &) const) &
                    Background::getImage<PixelT>,
            "interpStyle"_a, "undersampleStyle"_a = "THROW_EXCEPTION");
    cls.def(("getImage" + suffix).c_str(),
            (std::shared_ptr<lsst::afw::image::Image<PixelT>> (Background::*)(
                    lsst::geom::Box2I const &, Interpolate::Style const, UndersampleStyle const) const) &
                    Background::getImage<PixelT>,
            "bbox"_a, "interpStyle"_a, "undersampleStyle"_a = THROW_EXCEPTION);
    cls.def(("getImage" + suffix).c_str(),
            (std::shared_ptr<lsst::afw::image::Image<PixelT>> (Background::*)(
                    lsst::geom::Box2I const &, std::string const &, std::string const &) const) &
                    Background::getImage<PixelT>,
            "bbox"_a, "interpStyle"_a, "undersampleStyle"_a = "THROW_EXCEPTION");
    cls.def(("getImage" + suffix).c_str(),
            (std::shared_ptr<lsst::afw::image::Image<PixelT>> (Background::*)() const) &
                    Background::getImage<PixelT>);
}

/**
@internal Declare a warping kernel class with no constructor

@tparam KernelT  class of warping kernel, e.g. LanczosWarpingKernel
@param[in] mod  pybind11 module to which to add the kernel
@param[in] name  Python name for class, e.g. "LanczosWarpingKernel"
@param[in] addConstructor  If true then add a default constructor.
*/
template <typename KernelT>
py::class_<KernelT, std::shared_ptr<KernelT>, SeparableKernel> declareWarpingKernel(py::module &mod,
                                                                                    std::string const &name) {
    py::class_<KernelT, std::shared_ptr<KernelT>, SeparableKernel> cls(mod, name.c_str());

    cls.def("clone", &KernelT::clone);
    return cls;
}

/**
@internal Declare a warping kernel class with a defaut constructor

@tparam KernelT  class of warping kernel, e.g. LanczosWarpingKernel
@param[in] mod  pybind11 module to which to add the kernel
@param[in] name  Python name for class, e.g. "LanczosWarpingKernel"
@param[in] addConstructor  If true then add a default constructor.
*/
template <typename KernelT>
py::class_<KernelT, std::shared_ptr<KernelT>, SeparableKernel> declareSimpleWarpingKernel(
        py::module &mod, std::string const &name, bool addConstructor = true) {
    auto cls = declareWarpingKernel<KernelT>(mod, name);
    cls.def(py::init<>());
    return cls;
}

/**
@internal Declare wrappers for warpImage and warpCenteredImage
for a particular pair of image or masked image types

@tparam DestImageT  Desination image type, e.g. Image<int> or MaskedImage<float, MaskType, VarianceType>
@tparam SrcImageT  Source image type, e.g. Image<int> or MaskedImage<float, MaskType, VarianceType>
@param[in,out] mod  pybind11 module for which to declare the function wrappers
*/
template <typename DestImageT, typename SrcImageT>
void declareImageWarpingFunctions(py::module &mod) {
    auto const EdgePixel =
            edgePixel<DestImageT>(typename image::detail::image_traits<DestImageT>::image_category());
    mod.def("warpImage", (int (*)(DestImageT &, geom::SkyWcs const &, SrcImageT const &, geom::SkyWcs const &,
                                  WarpingControl const &, typename DestImageT::SinglePixel)) &
                                 warpImage<DestImageT, SrcImageT>,
            "destImage"_a, "destWcs"_a, "srcImage"_a, "srcWcs"_a, "control"_a, "padValue"_a = EdgePixel);

    mod.def("warpImage",
            (int (*)(DestImageT &, SrcImageT const &, geom::TransformPoint2ToPoint2 const &,
                     WarpingControl const &, typename DestImageT::SinglePixel)) &
                    warpImage<DestImageT, SrcImageT>,
            "destImage"_a, "srcImage"_a, "srcToDest"_a, "control"_a, "padValue"_a = EdgePixel);

    mod.def("warpCenteredImage", &warpCenteredImage<DestImageT, SrcImageT>, "destImage"_a, "srcImage"_a,
            "linearTransform"_a, "centerPoint"_a, "control"_a, "padValue"_a = EdgePixel);
}

/**
@internal Declare wrappers for warpExposure, warpImage and warpCenteredImage
for a particular pair of source and destination pixel types.

Declares both image and masked image variants of warpImage and warpCenteredImage.

@tparam DestPixelT  Desination pixel type, e.g. `int` or `float`
@tparam SrcPixelT  Source pixel type, e.g. `int` or `float`
@param[in,out] mod  pybind11 module for which to declare the function wrappers
*/
template <typename DestPixelT, typename SrcPixelT>
void declareWarpingFunctions(py::module &mod) {
    using DestExposureT = image::Exposure<DestPixelT, image::MaskPixel, image::VariancePixel>;
    using SrcExposureT = image::Exposure<SrcPixelT, image::MaskPixel, image::VariancePixel>;
    using DestImageT = image::Image<DestPixelT>;
    using SrcImageT = image::Image<SrcPixelT>;
    using DestMaskedImageT = image::MaskedImage<DestPixelT, image::MaskPixel, image::VariancePixel>;
    using SrcMaskedImageT = image::MaskedImage<SrcPixelT, image::MaskPixel, image::VariancePixel>;

    //mod.def("warpExposure", &warpExposure<DestExposureT, SrcExposureT>, "destExposure"_a, "srcExposure"_a,
    //        "control"_a, "padValue"_a = edgePixel<DestMaskedImageT>(
    //                             typename image::detail::image_traits<DestMaskedImageT>::image_category()));

    declareImageWarpingFunctions<DestImageT, SrcImageT>(mod);
    //declareImageWarpingFunctions<DestMaskedImageT, SrcMaskedImageT>(mod);
}

WRAP(TransformBoundedField) {
    ClsField cls(mod, "TransformBoundedField");

    cls.def(py::init<lsst::geom::Box2I const &, TransformBoundedField::Transform const &>(), "bbox"_a,
            "transform"_a);

    table::io::python::addPersistableMethods<TransformBoundedField>(cls);

    cls.def("__mul__", &TransformBoundedField::operator*, py::is_operator());
    cls.def("__eq__", &TransformBoundedField::operator==, py::is_operator());

    cls.def("getTransform", &TransformBoundedField::getTransform);
    cls.def("evaluate", (double (BoundedField::*)(double, double) const) & BoundedField::evaluate);
    cls.def("evaluate",
            (ndarray::Array<double, 1, 1>(TransformBoundedField::*)(
                    ndarray::Array<double const, 1> const &, ndarray::Array<double const, 1> const &) const) &
                    TransformBoundedField::evaluate);
    cls.def("evaluate", (double (TransformBoundedField::*)(lsst::geom::Point2D const &) const) &
                                TransformBoundedField::evaluate);

}
template <typename PixelT>
void declareStatisticsStack(py::module &mod) {
    mod.def("statisticsStack", (std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>(*)(
                                       lsst::afw::image::Image<PixelT> const &, Property, char,
                                       StatisticsControl const &))statisticsStack<PixelT>,
            "image"_a, "flags"_a, "dimensions"_a, "sctrl"_a = StatisticsControl());
    mod.def("statisticsStack", (std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>(*)(
                                       lsst::afw::image::MaskedImage<PixelT> const &, Property, char,
                                       StatisticsControl const &))statisticsStack<PixelT>,
            "image"_a, "flags"_a, "dimensions"_a, "sctrl"_a = StatisticsControl());
    mod.def("statisticsStack",
            (void (*)(lsst::afw::image::Image<PixelT> &,
                      std::vector<std::shared_ptr<lsst::afw::image::Image<PixelT>>> &, Property,
                      StatisticsControl const &,
                      std::vector<lsst::afw::image::VariancePixel> const &))statisticsStack<PixelT>,
            "out"_a, "images"_a, "flags"_a, "sctrl"_a = StatisticsControl(),
            "wvector"_a = std::vector<lsst::afw::image::VariancePixel>(0));
    mod.def("statisticsStack",
            (void (*)(lsst::afw::image::MaskedImage<PixelT> &,
                      std::vector<std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>> &, Property,
                      StatisticsControl const &,
                      std::vector<lsst::afw::image::VariancePixel> const &,
                      lsst::afw::image::MaskPixel, lsst::afw::image::MaskPixel))statisticsStack<PixelT>,
            "out"_a, "images"_a, "flags"_a, "sctrl"_a = StatisticsControl(),
            "wvector"_a = std::vector<lsst::afw::image::VariancePixel>(0),
            "clipped"_a=0, "excuse"_a=0);
    mod.def("statisticsStack",
            (void (*)(lsst::afw::image::MaskedImage<PixelT> &,
                      std::vector<std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>> &, Property,
                      StatisticsControl const &,
                      std::vector<lsst::afw::image::VariancePixel> const &,
                      lsst::afw::image::MaskPixel,
                      std::vector<std::pair<lsst::afw::image::MaskPixel, lsst::afw::image::MaskPixel>> const &
                ))statisticsStack<PixelT>,
            "out"_a, "images"_a, "flags"_a, "sctrl"_a, "wvector"_a, "clipped"_a, "maskMap"_a);
    mod.def("statisticsStack",
            (std::shared_ptr<lsst::afw::image::Image<PixelT>>(*)(
                    std::vector<std::shared_ptr<lsst::afw::image::Image<PixelT>>> &, Property,
                    StatisticsControl const &,
                    std::vector<lsst::afw::image::VariancePixel> const &))statisticsStack<PixelT>,
            "images"_a, "flags"_a, "sctrl"_a = StatisticsControl(),
            "wvector"_a = std::vector<lsst::afw::image::VariancePixel>(0));
    mod.def("statisticsStack",
            (std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>(*)(
                    std::vector<std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>> &, Property,
                    StatisticsControl const &,
                    std::vector<lsst::afw::image::VariancePixel> const &,
                    lsst::afw::image::MaskPixel, lsst::afw::image::MaskPixel))statisticsStack<PixelT>,
            "images"_a, "flags"_a, "sctrl"_a = StatisticsControl(),
            "wvector"_a = std::vector<lsst::afw::image::VariancePixel>(0),
            "clipped"_a=0, "excuse"_a=0);
    mod.def("statisticsStack",
            (std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>(*)(
                    std::vector<std::shared_ptr<lsst::afw::image::MaskedImage<PixelT>>> &, Property,
                    StatisticsControl const &,
                    std::vector<lsst::afw::image::VariancePixel> const &,
                    lsst::afw::image::MaskPixel,
                    std::vector<std::pair<lsst::afw::image::MaskPixel, lsst::afw::image::MaskPixel>> const &
                ))statisticsStack<PixelT>,
            "images"_a, "flags"_a, "sctrl"_a, "wvector"_a, "clipped"_a, "maskMap"_a);
    mod.def("statisticsStack",
            (std::vector<PixelT>(*)(
                    std::vector<std::vector<PixelT>> &, Property, StatisticsControl const &,
                    std::vector<lsst::afw::image::VariancePixel> const &))statisticsStack<PixelT>,
            "vectors"_a, "flags"_a, "sctrl"_a = StatisticsControl(),
            "wvector"_a = std::vector<lsst::afw::image::VariancePixel>(0));
}


template <typename ImageT>
void declareRandomImage(py::module &mod) {
    mod.def("randomUniformImage", (void (*)(ImageT *, Random &))randomUniformImage<ImageT>);
    mod.def("randomUniformPosImage", (void (*)(ImageT *, Random &))randomUniformPosImage<ImageT>);
    mod.def("randomUniformIntImage",
            (void (*)(ImageT *, Random &, unsigned long))randomUniformIntImage<ImageT>);
    mod.def("randomFlatImage",
            (void (*)(ImageT *, Random &, double const, double const))randomFlatImage<ImageT>);
    mod.def("randomGaussianImage", (void (*)(ImageT *, Random &))randomGaussianImage<ImageT>);
    mod.def("randomChisqImage", (void (*)(ImageT *, Random &, double const))randomChisqImage<ImageT>);
    mod.def("randomPoissonImage", (void (*)(ImageT *, Random &, double const))randomPoissonImage<ImageT>);
}

using PyClass = py::class_<PixelAreaBoundedField, std::shared_ptr<PixelAreaBoundedField>, BoundedField>;

WRAP(PixelAreaBound) {
    PyClass cls(mod, "PixelAreaBoundedField");
    cls.def(
        py::init<lsst::geom::Box2I const &, std::shared_ptr<afw::geom::SkyWcs const>,
                 lsst::geom::AngleUnit const &, double>(),
        "bbox"_a, "skyWcs"_a, "unit"_a=lsst::geom::radians, "scaling"_a=1.0
    );
    // All other operations are wrapped by the BoundedField base class.
}

template <typename ImageT>
static void declareOffsetImage(py::module& mod) {
    mod.def("offsetImage", offsetImage<ImageT>, "image"_a, "dx"_a, "dy"_a, "algorithmName"_a = "lanczos5",
            "buffer"_a = 0);
}

template <typename ImageT>
static void declareRotateImageBy90(py::module& mod) {
    mod.def("rotateImageBy90", rotateImageBy90<ImageT>, "image"_a, "nQuarter"_a);
}

template <typename ImageT>
static void declareFlipImage(py::module& mod) {
    mod.def("flipImage", flipImage<ImageT>, "inImage"_a, "flipLR"_a, "flipTB"_a);
}

template <typename ImageT>
static void declareBinImage(py::module& mod) {
    mod.def("binImage", (std::shared_ptr<ImageT>(*)(ImageT const&, int const, int const,
                                                    lsst::afw::math::Property const))binImage<ImageT>,
            "inImage"_a, "binX"_a, "binY"_a, "flags"_a = lsst::afw::math::MEAN);
    mod.def("binImage", (std::shared_ptr<ImageT>(*)(ImageT const&, int const,
                                                    lsst::afw::math::Property const))binImage<ImageT>,
            "inImage"_a, "binsize"_a, "flags"_a = lsst::afw::math::MEAN);
}


template <typename ReturnT>
void declareFunction(py::module &mod, std::string const &suffix) {
    auto const name = "Function" + suffix;

    py::class_<Function<ReturnT>, std::shared_ptr<Function<ReturnT>>> cls(mod, name.c_str());

    cls.def(py::init<unsigned int>(), "nParams"_a);
    cls.def(py::init<std::vector<double> const &>(), "params"_a);

    table::io::python::addPersistableMethods<Function<ReturnT>>(cls);

    cls.def("getNParameters", &Function<ReturnT>::getNParameters);
    cls.def("getParameters", &Function<ReturnT>::getParameters, py::return_value_policy::copy);
    cls.def("getParameter", &Function<ReturnT>::getParameter, "index"_a);
    cls.def("isLinearCombination", &Function<ReturnT>::isLinearCombination);
    cls.def("setParameter", &Function<ReturnT>::setParameter, "index"_a, "value"_a);
    cls.def("setParameters", &Function<ReturnT>::setParameters);
    cls.def("toString", &Function<ReturnT>::toString, "prefix"_a = "");
}

template <typename ReturnT>
void declareFunction1(py::module &mod, const std::string &suffix) {
    auto const name = "Function1" + suffix;

    py::class_<Function1<ReturnT>, std::shared_ptr<Function1<ReturnT>>, Function<ReturnT>> cls(mod,
                                                                                               name.c_str());

    table::io::python::addPersistableMethods<Function1<ReturnT>>(cls);

    cls.def("clone", &Function1<ReturnT>::clone);
    cls.def("__call__", &Function1<ReturnT>::operator(), "x"_a);
    cls.def("toString", &Function1<ReturnT>::toString, "prefix"_a = "");
    cls.def("computeCache", &Function1<ReturnT>::computeCache, "n"_a);
}

template <typename ReturnT>
void declareFunction2(py::module &mod, const std::string &suffix) {
    auto const name = "Function2" + suffix;

    py::class_<Function2<ReturnT>, std::shared_ptr<Function2<ReturnT>>, Function<ReturnT>> cls(mod,
                                                                                               name.c_str());

    table::io::python::addPersistableMethods<Function2<ReturnT>>(cls);

    cls.def("clone", &Function2<ReturnT>::clone);
    cls.def("__call__", &Function2<ReturnT>::operator(), "x"_a, "y"_a);
    cls.def("toString", &Function2<ReturnT>::toString, "prefix"_a = "");
    cls.def("getDFuncDParameters", &Function2<ReturnT>::getDFuncDParameters, "x"_a, "y"_a);
}

template <typename ReturnT>
void declareBasePolynomialFunction2(py::module &mod, const std::string &suffix) {
    auto const name = "BasePolynomialFunction2" + suffix;

    py::class_<BasePolynomialFunction2<ReturnT>, std::shared_ptr<BasePolynomialFunction2<ReturnT>>,
               Function2<ReturnT>>
            cls(mod, name.c_str());

    cls.def("getOrder", &BasePolynomialFunction2<ReturnT>::getOrder);
    cls.def("isLinearCombination", &BasePolynomialFunction2<ReturnT>::isLinearCombination);
    cls.def_static("nParametersFromOrder", &BasePolynomialFunction2<ReturnT>::nParametersFromOrder,
                   "order"_a);
    cls.def_static("orderFromNParameters", &BasePolynomialFunction2<ReturnT>::orderFromNParameters,
                   "nParameters"_a);
    cls.def("getDFuncDParameters", &BasePolynomialFunction2<ReturnT>::getDFuncDParameters, "x"_a, "y"_a);
}

template <typename ReturnT>
void declareNullFunction1(py::module &mod, const std::string &suffix) {
    auto const name = "NullFunction1" + suffix;

    py::class_<NullFunction1<ReturnT>, std::shared_ptr<NullFunction1<ReturnT>>, Function1<ReturnT>> cls(
            mod, name.c_str());

    cls.def(py::init<>());

    cls.def("clone", &NullFunction1<ReturnT>::clone);
}

template <typename ReturnT>
void declareNullFunction2(py::module &mod, const std::string &suffix) {
    auto const name = "NullFunction2" + suffix;

    py::class_<NullFunction2<ReturnT>, std::shared_ptr<NullFunction2<ReturnT>>, Function2<ReturnT>> cls(
            mod, name.c_str());

    cls.def(py::init<>());

    cls.def("clone", &NullFunction2<ReturnT>::clone);
}

template <typename ReturnT>
void declareAllFunctions(py::module &mod, const std::string &suffix) {
    declareFunction<ReturnT>(mod, suffix);
    declareFunction1<ReturnT>(mod, suffix);
    declareFunction2<ReturnT>(mod, suffix);
    declareBasePolynomialFunction2<ReturnT>(mod, suffix);
    declareNullFunction1<ReturnT>(mod, suffix);
    declareNullFunction2<ReturnT>(mod, suffix);
};

using PyBoundClass = py::class_<BoundedField, std::shared_ptr<BoundedField>>;

template <typename PixelT>
void declareTemplates(PyBoundClass &cls) {
    cls.def("fillImage", &BoundedField::fillImage<PixelT>,
            "image"_a, "overlapOnly"_a = false, "xStep"_a = 1, "yStep"_a = 1);
    cls.def("addToImage", &BoundedField::addToImage<PixelT>, "image"_a, "scaleBy"_a = 1.0,
            "overlapOnly"_a = false, "xStep"_a = 1, "yStep"_a = 1);
    cls.def("multiplyImage", &BoundedField::multiplyImage<PixelT>,
            "image"_a, "overlapOnly"_a = false, "xStep"_a = 1, "yStep"_a = 1);
    cls.def("divideImage", &BoundedField::divideImage<PixelT>,
            "image"_a, "overlapOnly"_a = false, "xStep"_a = 1, "yStep"_a = 1);
}

WRAP(BoundedField) {
    PyBoundClass cls(mod, "BoundedField");

    table::io::python::addPersistableMethods<BoundedField>(cls);

    cls.def("__rmul__", [](BoundedField &bf, double const scale) { return bf * scale; }, py::is_operator());
    cls.def("__mul__", &BoundedField::operator*, py::is_operator());
    cls.def("__truediv__", &BoundedField::operator/, py::is_operator());
    cls.def("__eq__", &BoundedField::operator==, py::is_operator());
    cls.def("__ne__", &BoundedField::operator!=, py::is_operator());

    cls.def("evaluate", (double (BoundedField::*)(double, double) const) & BoundedField::evaluate);
    cls.def("evaluate",
            (ndarray::Array<double, 1, 1> (BoundedField::*)(ndarray::Array<double const, 1> const &,
                                                            ndarray::Array<double const, 1> const &) const) &
                    BoundedField::evaluate);
    cls.def("evaluate",
            (double (BoundedField::*)(lsst::geom::Point2D const &) const) & BoundedField::evaluate);
    cls.def("integrate", &BoundedField::integrate);
    cls.def("mean", &BoundedField::mean);
    cls.def("getBBox", &BoundedField::getBBox);

    // Pybind11 resolves overloads by picking the first one that might work
    declareTemplates<double>(cls);
    declareTemplates<float>(cls);

    utils::python::addOutputOp(cls, "__str__");
    cls.def("__repr__", [](BoundedField const &self) {
        std::ostringstream os;
        os << "BoundedField(" << self << ")";
        return os.str();
    });
}

template <typename PixelT>
void declareApproximate(py::module &mod, std::string const &suffix) {
    using Class = Approximate<PixelT>;

    py::class_<Class, std::shared_ptr<Class>> cls(mod, ("Approximate" + suffix).c_str());

    cls.def("getImage", &Class::getImage, "orderX"_a = -1, "orderY"_a = -1);
    cls.def("getMaskedImage", &Class::getMaskedImage, "orderX"_a = -1, "orderY"_a = -1);

    mod.def("makeApproximate",
            (std::shared_ptr<Approximate<PixelT>>(*)(std::vector<double> const &, std::vector<double> const &,
                                                     image::MaskedImage<PixelT> const &, lsst::geom::Box2I const &,
                                                     ApproximateControl const &))makeApproximate<PixelT>,
            "x"_a, "y"_a, "im"_a, "bbox"_a, "ctrl"_a);
}

// PYBIND11_DECLARE_HOLDER_TYPE(MyType, std::shared_ptr<MyType>);
WRAP(Approximate) {
    py::class_<ApproximateControl, std::shared_ptr<ApproximateControl>> clsApproximateControl(
                mod, "ApproximateControl");

    py::enum_<ApproximateControl::Style>(clsApproximateControl, "Style")
            .value("UNKNOWN", ApproximateControl::Style::UNKNOWN)
            .value("CHEBYSHEV", ApproximateControl::Style::CHEBYSHEV)
            .value("NUM_STYLES", ApproximateControl::Style::NUM_STYLES)
            .export_values();

    clsApproximateControl.def(py::init<ApproximateControl::Style, int, int, bool>(), "style"_a, "orderX"_a,
                              "orderY"_a = -1, "weighting"_a = true);

    clsApproximateControl.def("getStyle", &ApproximateControl::getStyle);
    clsApproximateControl.def("setStyle", &ApproximateControl::setStyle);
    clsApproximateControl.def("getOrderX", &ApproximateControl::getOrderX);
    clsApproximateControl.def("setOrderX", &ApproximateControl::setOrderX);
    clsApproximateControl.def("getOrderY", &ApproximateControl::getOrderY);
    clsApproximateControl.def("setOrderY", &ApproximateControl::setOrderY);
    clsApproximateControl.def("getWeighting", &ApproximateControl::getWeighting);
    clsApproximateControl.def("setWeighting", &ApproximateControl::setWeighting);

    // Yes, really only float
    declareApproximate<float>(mod, "F");
}

template <typename PixelT>
void declareTemplates(ChebyClsField & cls) {
    cls.def_static("fit", (std::shared_ptr<ChebyshevBoundedField>(*)(lsst::afw::image::Image<PixelT> const &,
                                                                     ChebyshevBoundedFieldControl const &)) &
                                  ChebyshevBoundedField::fit);
}

WRAP(ChebyshexBoundedField) {
    /* Module level */
    py::class_<ChebyshevBoundedFieldControl> clsChebyshevBoundedFieldControl(mod,
                                                                             "ChebyshevBoundedFieldControl");
    ChebyClsField clsChebyshevBoundedField(mod, "ChebyshevBoundedField");

    /* Member types and enums */
    using Control = ChebyshevBoundedFieldControl;

    /* Constructors */
    clsChebyshevBoundedFieldControl.def(py::init<>());
    clsChebyshevBoundedField.def(
            py::init<lsst::geom::Box2I const &, ndarray::Array<double const, 2, 2> const &>());

    /* Operators are defined only in the BoundedField base class;
       we let Python inheritance provide them here. */

    /* Members */
    LSST_DECLARE_CONTROL_FIELD(clsChebyshevBoundedFieldControl, ChebyshevBoundedFieldControl, orderX);
    LSST_DECLARE_CONTROL_FIELD(clsChebyshevBoundedFieldControl, ChebyshevBoundedFieldControl, orderY);
    LSST_DECLARE_CONTROL_FIELD(clsChebyshevBoundedFieldControl, ChebyshevBoundedFieldControl, triangular);

    clsChebyshevBoundedFieldControl.def("computeSize", &ChebyshevBoundedFieldControl::computeSize);

    clsChebyshevBoundedField.def("getCoefficients", &ChebyshevBoundedField::getCoefficients);
    clsChebyshevBoundedField.def_static(
            "fit", (std::shared_ptr<ChebyshevBoundedField>(*)(
                           lsst::geom::Box2I const &, ndarray::Array<double const, 1> const &,
                           ndarray::Array<double const, 1> const &, ndarray::Array<double const, 1> const &,
                           Control const &)) &
                           ChebyshevBoundedField::fit);
    clsChebyshevBoundedField.def_static(
            "fit", (std::shared_ptr<ChebyshevBoundedField>(*)(
                           lsst::geom::Box2I const &, ndarray::Array<double const, 1> const &,
                           ndarray::Array<double const, 1> const &, ndarray::Array<double const, 1> const &,
                           ndarray::Array<double const, 1> const &, Control const &)) &
                           ChebyshevBoundedField::fit);

    clsChebyshevBoundedField.def("truncate", &ChebyshevBoundedField::truncate);

    // Pybind11 resolves overloads by picking the first one that might work
    declareTemplates<double>(clsChebyshevBoundedField);
    declareTemplates<float>(clsChebyshevBoundedField);
}

template <typename ReturnT>
void declareMinimize(py::module &mod) {
    mod.def("minimize", (FitResults(*)(lsst::afw::math::Function1<ReturnT> const &,
                                       std::vector<double> const &, std::vector<double> const &,
                                       std::vector<double> const &, std::vector<double> const &,
                                       std::vector<double> const &, double))minimize<ReturnT>);
    mod.def("minimize",
            (FitResults(*)(lsst::afw::math::Function2<ReturnT> const &, std::vector<double> const &,
                           std::vector<double> const &, std::vector<double> const &,
                           std::vector<double> const &, std::vector<double> const &,
                           std::vector<double> const &, double))minimize<ReturnT>);
};

WRAP(Minimize) {
    py::class_<FitResults> clsFitResults(mod, "FitResults");
    clsFitResults.def_readwrite("isValid", &FitResults::isValid);
    clsFitResults.def_readwrite("chiSq", &FitResults::chiSq);
    clsFitResults.def_readwrite("parameterList", &FitResults::parameterList);
    clsFitResults.def_readwrite("parameterErrorList", &FitResults::parameterErrorList);

    declareMinimize<double>(mod);
    declareMinimize<float>(mod);
}
template <typename T1, typename T2, int C1, int C2>
void declareLeastSquares(py::module &mod) {
    py::class_<LeastSquares> cls(mod, "LeastSquares");
    py::enum_<LeastSquares::Factorization>(cls, "Factorization")
            .value("NORMAL_EIGENSYSTEM", LeastSquares::Factorization::NORMAL_EIGENSYSTEM)
            .value("NORMAL_CHOLESKY", LeastSquares::Factorization::NORMAL_CHOLESKY)
            .value("DIRECT_SVD", LeastSquares::Factorization::DIRECT_SVD)
            .export_values();
    cls.def_static("fromDesignMatrix",
                   (LeastSquares(*)(ndarray::Array<T1, 2, C1> const &, ndarray::Array<T2, 1, C2> const &,
                                    LeastSquares::Factorization)) &
                           LeastSquares::fromDesignMatrix<T1, T2, C1, C2>,
                   "design"_a, "data"_a, "factorization"_a = LeastSquares::NORMAL_EIGENSYSTEM);
    cls.def_static("fromNormalEquations",
                   (LeastSquares(*)(ndarray::Array<T1, 2, C1> const &, ndarray::Array<T2, 1, C2> const &,
                                    LeastSquares::Factorization)) &
                           LeastSquares::fromNormalEquations<T1, T2, C1, C2>,
                   "fisher"_a, "rhs"_a, "factorization"_a = LeastSquares::NORMAL_EIGENSYSTEM);
    cls.def("getRank", &LeastSquares::getRank);
    cls.def("setDesignMatrix",
            (void (LeastSquares::*)(ndarray::Array<T1, 2, C1> const &, ndarray::Array<T2, 1, C2> const &)) &
                    LeastSquares::setDesignMatrix<T1, T2, C1, C2>);
    cls.def("getDimension", &LeastSquares::getDimension);
    cls.def("setNormalEquations",
            (void (LeastSquares::*)(ndarray::Array<T1, 2, C1> const &, ndarray::Array<T2, 1, C2> const &)) &
                    LeastSquares::setNormalEquations<T1, T2, C1, C2>);
    cls.def("getSolution", &LeastSquares::getSolution);
    cls.def("getFisherMatrix", &LeastSquares::getFisherMatrix);
    cls.def("getCovariance", &LeastSquares::getCovariance);
    cls.def("getFactorization", &LeastSquares::getFactorization);
    cls.def("getDiagnostic", &LeastSquares::getDiagnostic);
    cls.def("getThreshold", &LeastSquares::getThreshold);
    cls.def("setThreshold", &LeastSquares::setThreshold);
};

template <typename T>
void declareKdTree(py::module &mod, const std::string &suffix) {
    py::class_<KdTree<T>> clsKdTree(mod, ("KdTree" + suffix).c_str());
    clsKdTree.def(py::init<>());
    clsKdTree.def("Initialize", &KdTree<T>::Initialize);
    clsKdTree.def("removePoint", &KdTree<T>::removePoint);
    clsKdTree.def("getData", (T (KdTree<T>::*)(int, int) const) & KdTree<T>::getData);
    clsKdTree.def("getData", (ndarray::Array<T, 1, 1> (KdTree<T>::*)(int) const) & KdTree<T>::getData);
    clsKdTree.def("addPoint", &KdTree<T>::addPoint);
    clsKdTree.def("getNPoints", &KdTree<T>::getNPoints);
    clsKdTree.def("getTreeNode", &KdTree<T>::getTreeNode);
    clsKdTree.def("findNeighbors", &KdTree<T>::findNeighbors);
};

template <typename T>
void declareCovariograms(py::module &mod, const std::string &suffix) {
    /* Covariogram */
    py::class_<Covariogram<T>, std::shared_ptr<Covariogram<T>>> clsCovariogram(
            mod, ("Covariogram" + suffix).c_str());
    clsCovariogram.def(py::init<>());
    clsCovariogram.def("__call__", &Covariogram<T>::operator());

    /* SquaredExpCovariogram */
    py::class_<SquaredExpCovariogram<T>, std::shared_ptr<SquaredExpCovariogram<T>>, Covariogram<T>>
            clsSquaredExpCovariogram(mod, ("SquaredExpCovariogram" + suffix).c_str());
    clsSquaredExpCovariogram.def(py::init<>());
    clsSquaredExpCovariogram.def("__call__", &SquaredExpCovariogram<T>::operator());
    clsSquaredExpCovariogram.def("setEllSquared", &SquaredExpCovariogram<T>::setEllSquared);

    /* SquaredExpCovariogram */
    py::class_<NeuralNetCovariogram<T>, std::shared_ptr<NeuralNetCovariogram<T>>, Covariogram<T>>
            clsNeuralNetCovariogram(mod, ("NeuralNetCovariogram" + suffix).c_str());
    clsNeuralNetCovariogram.def(py::init<>());
    clsNeuralNetCovariogram.def("setSigma0", &NeuralNetCovariogram<T>::setSigma0);
    clsNeuralNetCovariogram.def("setSigma1", &NeuralNetCovariogram<T>::setSigma1);
};

template <typename T>
void declareGaussianProcess(py::module &mod, const std::string &suffix) {
    py::class_<GaussianProcess<T>> clsGaussianProcess(mod, ("GaussianProcess" + suffix).c_str());
    /* Constructors */
    clsGaussianProcess.def(py::init<ndarray::Array<T, 2, 2> const &, ndarray::Array<T, 1, 1> const &,
                                    std::shared_ptr<Covariogram<T>> const &>());
    clsGaussianProcess.def(py::init<ndarray::Array<T, 2, 2> const &, ndarray::Array<T, 1, 1> const &,
                                    ndarray::Array<T, 1, 1> const &, ndarray::Array<T, 1, 1> const &,
                                    std::shared_ptr<Covariogram<T>> const &>());
    clsGaussianProcess.def(py::init<ndarray::Array<T, 2, 2> const &, ndarray::Array<T, 2, 2> const &,
                                    std::shared_ptr<Covariogram<T>> const &>());
    clsGaussianProcess.def(py::init<ndarray::Array<T, 2, 2> const &, ndarray::Array<T, 1, 1> const &,
                                    ndarray::Array<T, 1, 1> const &, ndarray::Array<T, 2, 2> const &,
                                    std::shared_ptr<Covariogram<T>> const &>());
    /* Members */
    clsGaussianProcess.def(
            "interpolate",
            (T (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1>, ndarray::Array<T, 1, 1> const &, int) const) &
                    GaussianProcess<T>::interpolate);
    clsGaussianProcess.def("interpolate",
                           (void (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1>, ndarray::Array<T, 1, 1>,
                                                         ndarray::Array<T, 1, 1> const &, int) const) &
                                   GaussianProcess<T>::interpolate);
    clsGaussianProcess.def("selfInterpolate",
                           (T (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1>, int, int) const) &
                                   GaussianProcess<T>::selfInterpolate);
    clsGaussianProcess.def(
            "selfInterpolate",
            (void (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1>, ndarray::Array<T, 1, 1>, int, int) const) &
                    GaussianProcess<T>::selfInterpolate);
    clsGaussianProcess.def("setLambda", &GaussianProcess<T>::setLambda);
    clsGaussianProcess.def("setCovariogram", &GaussianProcess<T>::setCovariogram);
    clsGaussianProcess.def("addPoint", (void (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1> const &, T)) &
                                               GaussianProcess<T>::addPoint);
    clsGaussianProcess.def("addPoint", (void (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1> const &,
                                                                     ndarray::Array<T, 1, 1> const &)) &
                                               GaussianProcess<T>::addPoint);
    clsGaussianProcess.def("batchInterpolate",
                           (void (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1>, ndarray::Array<T, 1, 1>,
                                                         ndarray::Array<T, 2, 2> const &) const) &
                                   GaussianProcess<T>::batchInterpolate);
    clsGaussianProcess.def(
            "batchInterpolate",
            (void (GaussianProcess<T>::*)(ndarray::Array<T, 1, 1>, ndarray::Array<T, 2, 2> const &) const) &
                    GaussianProcess<T>::batchInterpolate);
    clsGaussianProcess.def("batchInterpolate",
                           (void (GaussianProcess<T>::*)(ndarray::Array<T, 2, 2>, ndarray::Array<T, 2, 2>,
                                                         ndarray::Array<T, 2, 2> const &) const) &
                                   GaussianProcess<T>::batchInterpolate);
    clsGaussianProcess.def(
            "batchInterpolate",
            (void (GaussianProcess<T>::*)(ndarray::Array<T, 2, 2>, ndarray::Array<T, 2, 2> const &) const) &
                    GaussianProcess<T>::batchInterpolate);
    clsGaussianProcess.def("setKrigingParameter", &GaussianProcess<T>::setKrigingParameter);
    clsGaussianProcess.def("removePoint", &GaussianProcess<T>::removePoint);
    clsGaussianProcess.def("getNPoints", &GaussianProcess<T>::getNPoints);
    clsGaussianProcess.def("getData",
                           (void (GaussianProcess<T>::*)(ndarray::Array<T, 2, 2>, ndarray::Array<T, 1, 1>,
                                                         ndarray::Array<int, 1, 1>) const) &
                                   GaussianProcess<T>::getData);
    clsGaussianProcess.def("getData",
                           (void (GaussianProcess<T>::*)(ndarray::Array<T, 2, 2>, ndarray::Array<T, 2, 2>,
                                                         ndarray::Array<int, 1, 1>) const) &
                                   GaussianProcess<T>::getData);
};

WRAP(Gaussian) {
    declareCovariograms<double>(mod, "D");
    declareGaussianProcess<double>(mod, "D");
    declareKdTree<double>(mod, "D");
}

WRAP(Function) {
    declareAllFunctions<float>(mod, "F");
    declareAllFunctions<double>(mod, "D");
}

template <typename ReturnT>
void declarePolynomialFunctions(py::module &mod, const std::string &suffix) {
    /* PolynomialFunction1 */
    py::class_<PolynomialFunction1<ReturnT>, std::shared_ptr<PolynomialFunction1<ReturnT>>,
               Function1<ReturnT>>
            clsPolynomialFunction1(mod, ("PolynomialFunction1" + suffix).c_str());

    clsPolynomialFunction1.def(py::init<std::vector<double> const &>(), "params"_a);
    clsPolynomialFunction1.def(py::init<unsigned int>(), "order"_a);

    clsPolynomialFunction1.def("__call__", &PolynomialFunction1<ReturnT>::operator(), "x"_a);
    clsPolynomialFunction1.def("clone", &PolynomialFunction1<ReturnT>::clone);
    clsPolynomialFunction1.def("isLinearCombination", &PolynomialFunction1<ReturnT>::isLinearCombination);
    clsPolynomialFunction1.def("getOrder", &PolynomialFunction1<ReturnT>::getOrder);
    clsPolynomialFunction1.def("toString", &PolynomialFunction1<ReturnT>::toString, "prefix"_a = "");

    /* PolynomialFunction2 */
    py::class_<PolynomialFunction2<ReturnT>, std::shared_ptr<PolynomialFunction2<ReturnT>>,
               BasePolynomialFunction2<ReturnT>>
            clsPolynomialFunction2(mod, ("PolynomialFunction2" + suffix).c_str());

    clsPolynomialFunction2.def(py::init<std::vector<double> const &>(), "params"_a);
    clsPolynomialFunction2.def(py::init<unsigned int>(), "order"_a);

    clsPolynomialFunction2.def("__call__", &PolynomialFunction2<ReturnT>::operator(), "x"_a, "y"_a);
    clsPolynomialFunction2.def("clone", &PolynomialFunction2<ReturnT>::clone);
    clsPolynomialFunction2.def("getOrder", &PolynomialFunction2<ReturnT>::getOrder);
    clsPolynomialFunction2.def("getDFuncDParameters", &PolynomialFunction2<ReturnT>::getDFuncDParameters);
    clsPolynomialFunction2.def("toString", &PolynomialFunction2<ReturnT>::toString, "prefix"_a = "");
    clsPolynomialFunction2.def("isPersistable", &PolynomialFunction2<ReturnT>::isPersistable);
};

template <typename ReturnT>
void declareChebyshevFunctions(py::module &mod, const std::string &suffix) {
    /* Chebyshev1Function1 */
    py::class_<Chebyshev1Function1<ReturnT>, std::shared_ptr<Chebyshev1Function1<ReturnT>>,
               Function1<ReturnT>>
            clsChebyshev1Function1(mod, ("Chebyshev1Function1" + suffix).c_str());

    clsChebyshev1Function1.def(py::init<std::vector<double>, double, double>(), "params"_a, "minX"_a = -1,
                               "maxX"_a = 1);
    clsChebyshev1Function1.def(py::init<unsigned int, double, double>(), "order"_a, "minX"_a = -1,
                               "maxX"_a = 1);

    clsChebyshev1Function1.def("__call__", &Chebyshev1Function1<ReturnT>::operator(), "x"_a);
    clsChebyshev1Function1.def("clone", &Chebyshev1Function1<ReturnT>::clone);
    clsChebyshev1Function1.def("getMinX", &Chebyshev1Function1<ReturnT>::getMinX);
    clsChebyshev1Function1.def("getMaxX", &Chebyshev1Function1<ReturnT>::getMaxX);
    clsChebyshev1Function1.def("getOrder", &Chebyshev1Function1<ReturnT>::getOrder);
    clsChebyshev1Function1.def("isLinearCombination", &Chebyshev1Function1<ReturnT>::isLinearCombination);
    clsChebyshev1Function1.def("toString", &Chebyshev1Function1<ReturnT>::toString, "prefix"_a = "");

    /* Chebyshev1Function2 */
    py::class_<Chebyshev1Function2<ReturnT>, std::shared_ptr<Chebyshev1Function2<ReturnT>>,
               BasePolynomialFunction2<ReturnT>>
            clsChebyshev1Function2(mod, ("Chebyshev1Function2" + suffix).c_str());

    clsChebyshev1Function2.def(py::init<std::vector<double>, lsst::geom::Box2D const &>(), "params"_a,
                               "xyRange"_a = lsst::geom::Box2D(lsst::geom::Point2D(-1.0, -1.0), lsst::geom::Point2D(1.0, 1.0)));
    clsChebyshev1Function2.def(py::init<unsigned int, lsst::geom::Box2D const &>(), "order"_a,
                               "xyRange"_a = lsst::geom::Box2D(lsst::geom::Point2D(-1.0, -1.0), lsst::geom::Point2D(1.0, 1.0)));

    clsChebyshev1Function2.def("__call__", &Chebyshev1Function2<ReturnT>::operator(), "x"_a, "y"_a);
    clsChebyshev1Function2.def("clone", &Chebyshev1Function2<ReturnT>::clone);
    clsChebyshev1Function2.def("getXYRange", &Chebyshev1Function2<ReturnT>::getXYRange);
    clsChebyshev1Function2.def("truncate", &Chebyshev1Function2<ReturnT>::truncate, "order"_a);
    clsChebyshev1Function2.def("toString", &Chebyshev1Function2<ReturnT>::toString, "prefix"_a = "");
    clsChebyshev1Function2.def("isPersistable", &Chebyshev1Function2<ReturnT>::isPersistable);
};

template <typename ReturnT>
void declareGaussianFunctions(py::module &mod, const std::string &suffix) {
    /* GaussianFunction1 */
    py::class_<GaussianFunction1<ReturnT>, std::shared_ptr<GaussianFunction1<ReturnT>>, Function1<ReturnT>>
            clsGaussianFunction1(mod, ("GaussianFunction1" + suffix).c_str());

    clsGaussianFunction1.def(py::init<double>(), "sigma"_a);

    clsGaussianFunction1.def("__call__", &GaussianFunction1<ReturnT>::operator(), "x"_a);
    clsGaussianFunction1.def("clone", &GaussianFunction1<ReturnT>::clone);
    clsGaussianFunction1.def("toString", &GaussianFunction1<ReturnT>::toString, "prefix"_a = "");

    /* GaussianFunction2 */
    py::class_<GaussianFunction2<ReturnT>, std::shared_ptr<GaussianFunction2<ReturnT>>, Function2<ReturnT>>
            clsGaussianFunction2(mod, ("GaussianFunction2" + suffix).c_str());

    clsGaussianFunction2.def(py::init<double, double, double>(), "sigma1"_a, "sigma2"_a, "angle"_a = 0.0);

    clsGaussianFunction2.def("__call__", &GaussianFunction2<ReturnT>::operator(), "x"_a, "y"_a);
    clsGaussianFunction2.def("clone", &GaussianFunction2<ReturnT>::clone);
    clsGaussianFunction2.def("toString", &GaussianFunction2<ReturnT>::toString, "prefix"_a = "");
    clsGaussianFunction2.def("isPersistable", &GaussianFunction2<ReturnT>::isPersistable);

    /* DoubleGaussianFunction2 */
    py::class_<DoubleGaussianFunction2<ReturnT>, std::shared_ptr<DoubleGaussianFunction2<ReturnT>>,
               Function2<ReturnT>>
            clsDoubleGaussianFunction2(mod, ("DoubleGaussianFunction2" + suffix).c_str());

    clsDoubleGaussianFunction2.def(py::init<double, double, double>(), "sigma1"_a, "sigma2"_a = 0,
                                   "ampl"_a = 0);

    clsDoubleGaussianFunction2.def("__call__", &DoubleGaussianFunction2<ReturnT>::operator(), "x"_a, "y"_a);
    clsDoubleGaussianFunction2.def("clone", &DoubleGaussianFunction2<ReturnT>::clone);
    clsDoubleGaussianFunction2.def("toString", &DoubleGaussianFunction2<ReturnT>::toString, "prefix"_a = "");
    clsDoubleGaussianFunction2.def("isPersistable", &DoubleGaussianFunction2<ReturnT>::isPersistable);
};

template <typename ReturnT>
void declareIntegerDeltaFunctions(py::module &mod, const std::string &suffix) {
    /* IntegerDeltaFunction1 */
    py::class_<IntegerDeltaFunction1<ReturnT>, std::shared_ptr<IntegerDeltaFunction1<ReturnT>>,
               Function1<ReturnT>>
            clsIntegerDeltaFunction1(mod, ("IntegerDeltaFunction1" + suffix).c_str());

    clsIntegerDeltaFunction1.def(py::init<double>(), "xo"_a);

    clsIntegerDeltaFunction1.def("__call__", &IntegerDeltaFunction1<ReturnT>::operator(), "x"_a);
    clsIntegerDeltaFunction1.def("clone", &IntegerDeltaFunction1<ReturnT>::clone);
    clsIntegerDeltaFunction1.def("toString", &IntegerDeltaFunction1<ReturnT>::toString, "prefix"_a = "");

    /* IntegerDeltaFunction2 */
    py::class_<IntegerDeltaFunction2<ReturnT>, std::shared_ptr<IntegerDeltaFunction2<ReturnT>>,
               Function2<ReturnT>>
            clsIntegerDeltaFunction2(mod, ("IntegerDeltaFunction2" + suffix).c_str());

    clsIntegerDeltaFunction2.def(py::init<double, double>(), "xo"_a, "yo"_a);

    clsIntegerDeltaFunction2.def("__call__", &IntegerDeltaFunction2<ReturnT>::operator(), "x"_a, "y"_a);
    clsIntegerDeltaFunction2.def("clone", &IntegerDeltaFunction2<ReturnT>::clone);
    clsIntegerDeltaFunction2.def("toString", &IntegerDeltaFunction2<ReturnT>::toString, "prefix"_a = "");
};

template <typename ReturnT>
void declareLanczosFunctions(py::module &mod, const std::string &suffix) {
    /* LanczosFunction1 */
    py::class_<LanczosFunction1<ReturnT>, std::shared_ptr<LanczosFunction1<ReturnT>>, Function1<ReturnT>>
            clsLanczosFunction1(mod, ("LanczosFunction1" + suffix).c_str());

    clsLanczosFunction1.def(py::init<unsigned int, double>(), "n"_a, "xOffset"_a = 0.0);

    clsLanczosFunction1.def("__call__", &LanczosFunction1<ReturnT>::operator(), "x"_a);
    clsLanczosFunction1.def("clone", &LanczosFunction1<ReturnT>::clone);
    clsLanczosFunction1.def("getOrder", &LanczosFunction1<ReturnT>::getOrder);
    clsLanczosFunction1.def("toString", &LanczosFunction1<ReturnT>::toString, "prefix"_a = "");

    /* LanczosFunction2 */
    py::class_<LanczosFunction2<ReturnT>, std::shared_ptr<LanczosFunction2<ReturnT>>, Function2<ReturnT>>
            clsLanczosFunction2(mod, ("LanczosFunction2" + suffix).c_str());

    clsLanczosFunction2.def(py::init<unsigned int, double, double>(), "n"_a, "xOffset"_a = 0.0,
                            "yOffset"_a = 0.0);

    clsLanczosFunction2.def("__call__", &LanczosFunction2<ReturnT>::operator(), "x"_a, "y"_a);
    clsLanczosFunction2.def("clone", &LanczosFunction2<ReturnT>::clone);
    clsLanczosFunction2.def("getOrder", &LanczosFunction2<ReturnT>::getOrder);
    clsLanczosFunction2.def("toString", &LanczosFunction2<ReturnT>::toString, "prefix"_a = "");
};

WRAP(FunctionLibrary) {

    declarePolynomialFunctions<float>(mod, "F");
    declareChebyshevFunctions<float>(mod, "F");
    declareGaussianFunctions<float>(mod, "F");
    declareIntegerDeltaFunctions<float>(mod, "F");
    declareLanczosFunctions<float>(mod, "F");

    declarePolynomialFunctions<double>(mod, "D");
    declareChebyshevFunctions<double>(mod, "D");
    declareGaussianFunctions<double>(mod, "D");
    declareIntegerDeltaFunctions<double>(mod, "D");
    declareLanczosFunctions<double>(mod, "D");
}

WRAP(LeastSquares) {
    declareLeastSquares<double, double, 0, 0>(mod);
}

WRAP(Interpolate) {
    py::class_<Interpolate, std::shared_ptr<Interpolate>> clsInterpolate(mod, "Interpolate");
    py::enum_<Interpolate::Style>(clsInterpolate, "Style")
            .value("UNKNOWN", Interpolate::Style::UNKNOWN)
            .value("CONSTANT", Interpolate::Style::CONSTANT)
            .value("LINEAR", Interpolate::Style::LINEAR)
            .value("NATURAL_SPLINE", Interpolate::Style::NATURAL_SPLINE)
            .value("CUBIC_SPLINE", Interpolate::Style::CUBIC_SPLINE)
            .value("CUBIC_SPLINE_PERIODIC", Interpolate::Style::CUBIC_SPLINE_PERIODIC)
            .value("AKIMA_SPLINE", Interpolate::Style::AKIMA_SPLINE)
            .value("AKIMA_SPLINE_PERIODIC", Interpolate::Style::AKIMA_SPLINE_PERIODIC)
            .value("NUM_STYLES", Interpolate::Style::NUM_STYLES)
            .export_values();

    clsInterpolate.def("interpolate", [](Interpolate &t, double const x) {
        /*
        We use a lambda function here because interpolate (with a double) is a virtual function and therefor
        cannot be wrapped directly.
        */
        return t.interpolate(x);
    });

    clsInterpolate.def("interpolate",
                       (std::vector<double> (Interpolate::*)(std::vector<double> const &) const) &
                               Interpolate::interpolate);
    clsInterpolate.def(
            "interpolate",
            (ndarray::Array<double, 1> (Interpolate::*)(ndarray::Array<double const, 1> const &) const) &
                    Interpolate::interpolate);

    mod.def("makeInterpolate",
            (std::shared_ptr<Interpolate>(*)(std::vector<double> const &, std::vector<double> const &,
                                             Interpolate::Style const))makeInterpolate,
            "x"_a, "y"_a, "style"_a = Interpolate::AKIMA_SPLINE);
    mod.def("makeInterpolate", (std::shared_ptr<Interpolate>(*)(ndarray::Array<double const, 1> const &,
                                                                ndarray::Array<double const, 1> const &y,
                                                                Interpolate::Style const))makeInterpolate,
            "x"_a, "y"_a, "style"_a = Interpolate::AKIMA_SPLINE);

    mod.def("stringToInterpStyle", stringToInterpStyle, "style"_a);
    mod.def("lookupMaxInterpStyle", lookupMaxInterpStyle, "n"_a);
    mod.def("lookupMinInterpPoints", lookupMinInterpPoints, "style"_a);

}
WRAP(OffsetImage) {
    using MaskPixel = lsst::afw::image::MaskPixel;

    /* Module level */
    declareOffsetImage<lsst::afw::image::Image<int>>(mod);
    declareOffsetImage<lsst::afw::image::Image<float>>(mod);
    declareOffsetImage<lsst::afw::image::Image<double>>(mod);
    declareOffsetImage<lsst::afw::image::MaskedImage<int>>(mod);
    declareOffsetImage<lsst::afw::image::MaskedImage<float>>(mod);
    declareOffsetImage<lsst::afw::image::MaskedImage<double>>(mod);

    declareRotateImageBy90<lsst::afw::image::Image<std::uint16_t>>(mod);
    declareRotateImageBy90<lsst::afw::image::Image<int>>(mod);
    declareRotateImageBy90<lsst::afw::image::Image<float>>(mod);
    declareRotateImageBy90<lsst::afw::image::Image<double>>(mod);
    declareRotateImageBy90<lsst::afw::image::MaskedImage<std::uint16_t>>(mod);
    declareRotateImageBy90<lsst::afw::image::MaskedImage<int>>(mod);
    declareRotateImageBy90<lsst::afw::image::MaskedImage<float>>(mod);
    declareRotateImageBy90<lsst::afw::image::MaskedImage<double>>(mod);
    declareRotateImageBy90<lsst::afw::image::Mask<MaskPixel>>(mod);

    declareFlipImage<lsst::afw::image::Image<std::uint16_t>>(mod);
    declareFlipImage<lsst::afw::image::Image<int>>(mod);
    declareFlipImage<lsst::afw::image::Image<float>>(mod);
    declareFlipImage<lsst::afw::image::Image<double>>(mod);
    declareFlipImage<lsst::afw::image::MaskedImage<std::uint16_t>>(mod);
    declareFlipImage<lsst::afw::image::MaskedImage<int>>(mod);
    declareFlipImage<lsst::afw::image::MaskedImage<float>>(mod);
    declareFlipImage<lsst::afw::image::MaskedImage<double>>(mod);
    declareFlipImage<lsst::afw::image::Mask<MaskPixel>>(mod);

    declareBinImage<lsst::afw::image::Image<std::uint16_t>>(mod);
    declareBinImage<lsst::afw::image::Image<int>>(mod);
    declareBinImage<lsst::afw::image::Image<float>>(mod);
    declareBinImage<lsst::afw::image::Image<double>>(mod);
    declareBinImage<lsst::afw::image::MaskedImage<std::uint16_t>>(mod);
    declareBinImage<lsst::afw::image::MaskedImage<int>>(mod);
    declareBinImage<lsst::afw::image::MaskedImage<float>>(mod);
    declareBinImage<lsst::afw::image::MaskedImage<double>>(mod);
}
WRAP(TestClasses) {
    /*
     * Test class for SpatialCellCandidate
     */
    class TestCandidate : public SpatialCellCandidate {
    public:
        TestCandidate(float const xCenter,  ///< @internal The object's column-centre
                      float const yCenter,  ///< @internal The object's row-centre
                      float const flux      ///< @internal The object's flux
                      )
                : SpatialCellCandidate(xCenter, yCenter), _flux(flux) {}

        /// @internal Return candidates rating
        virtual double getCandidateRating() const { return _flux; }
        virtual void setCandidateRating(double flux) { _flux = flux; }

    private:
        double _flux;
    };

    /// @internal A class to pass around to all our TestCandidates
    class TestCandidateVisitor : public CandidateVisitor {
    public:
        TestCandidateVisitor() : CandidateVisitor(), _n(0) {}

        // Called by SpatialCellSet::visitCandidates before visiting any Candidates
        void reset() { _n = 0; }

        // Called by SpatialCellSet::visitCandidates for each Candidate
        void processCandidate(SpatialCellCandidate *candidate) { ++_n; }

        int getN() const { return _n; }

    private:
        int _n;  // number of TestCandidates
    };

    class TestImageCandidate : public SpatialCellImageCandidate {
    public:
        typedef image::MaskedImage<float> MaskedImageT;

        TestImageCandidate(float const xCenter,  ///< @internal The object's column-centre
                           float const yCenter,  ///< @internal The object's row-centre
                           float const flux      ///< @internal The object's flux
                           )
                : SpatialCellImageCandidate(xCenter, yCenter), _flux(flux) {}

        /// @internal Return candidates rating
        double getCandidateRating() const { return _flux; }

        /// @internal Return the %image
        std::shared_ptr<MaskedImageT const> getMaskedImage() const {
            if (!_image) {
                _image = std::make_shared<MaskedImageT>(lsst::geom::ExtentI(getWidth(), getHeight()));
                *_image->getImage() = _flux;
            }
            return _image;
        }

    private:
        mutable std::shared_ptr<MaskedImageT> _image;
        double _flux;
    };

    py::class_<TestCandidate, std::shared_ptr<TestCandidate>, SpatialCellCandidate> clsTestCandidate(
            mod, ("TestCandidate"));
    clsTestCandidate.def(py::init<float const, float const, float const>());
    clsTestCandidate.def("getCandidateRating", &TestCandidate::getCandidateRating);
    clsTestCandidate.def("setCandidateRating", &TestCandidate::setCandidateRating);

    py::class_<TestCandidateVisitor, std::shared_ptr<TestCandidateVisitor>, CandidateVisitor>
            clsTestCandidateVisitor(mod, ("TestCandidateVisitor"));
    clsTestCandidateVisitor.def(py::init<>());
    clsTestCandidateVisitor.def("getN", &TestCandidateVisitor::getN);

    py::class_<TestImageCandidate, std::shared_ptr<TestImageCandidate>, SpatialCellImageCandidate>
            clsTestImageCandidate(mod, "TestImageCandidate");
    clsTestImageCandidate.def(py::init<float const, float const, float const>(), "xCenter"_a, "yCenter"_a,
                              "flux"_a);
    clsTestImageCandidate.def("getCandidateRating", &TestImageCandidate::getCandidateRating);
    clsTestImageCandidate.def("getMaskedImage", &TestImageCandidate::getMaskedImage);
};

WRAP(Spatial) {
    wrapSpatialCellCandidate(mod);
    wrapSpatialCellCandidateIterator(mod);
    wrapSpatialCell(mod);
    wrapSpatialCellSet(mod);
    wrapCandidateVisitor(mod);
    wrapSpatialCellImageCandidate(mod);

    /* Test Members */
    wrapTestClasses(mod);
}

WRAP(Kernel) {
    py::class_<Kernel, std::shared_ptr<Kernel>> clsKernel(mod, "Kernel");

    lsst::afw::table::io::python::addPersistableMethods<Kernel>(clsKernel);

    clsKernel.def("clone", &Kernel::clone);
    clsKernel.def("resized", &Kernel::resized, "width"_a, "height"_a);
    clsKernel.def("computeImage", &Kernel::computeImage, "image"_a, "doNormalize"_a, "x"_a = 0.0,
                  "y"_a = 0.0);
    clsKernel.def("getDimensions", &Kernel::getDimensions);
    clsKernel.def("setDimensions", &Kernel::setDimensions);
    clsKernel.def("setWidth", &Kernel::setWidth);
    clsKernel.def("setHeight", &Kernel::setHeight);
    clsKernel.def("getWidth", &Kernel::getWidth);
    clsKernel.def("getHeight", &Kernel::getHeight);
    clsKernel.def("getCtr", &Kernel::getCtr);
    clsKernel.def("getBBox", &Kernel::getBBox);
    clsKernel.def("getNKernelParameters", &Kernel::getNKernelParameters);
    clsKernel.def("getNSpatialParameters", &Kernel::getNSpatialParameters);
    clsKernel.def("getSpatialFunction", &Kernel::getSpatialFunction);
    clsKernel.def("getSpatialFunctionList", &Kernel::getSpatialFunctionList);
    clsKernel.def("getKernelParameter", &Kernel::getKernelParameter);
    clsKernel.def("getKernelParameters", &Kernel::getKernelParameters);
    clsKernel.def("growBBox", &Kernel::growBBox);
    clsKernel.def("shrinkBBox", &Kernel::shrinkBBox);
    clsKernel.def("setCtr", &Kernel::setCtr);
    clsKernel.def("getSpatialParameters", &Kernel::getSpatialParameters);
    clsKernel.def("isSpatiallyVarying", &Kernel::isSpatiallyVarying);
    clsKernel.def("setKernelParameters",
                  (void (Kernel::*)(std::vector<double> const &)) & Kernel::setKernelParameters);
    clsKernel.def("setKernelParameters",
                  (void (Kernel::*)(std::pair<double, double> const &)) & Kernel::setKernelParameters);
    clsKernel.def("setSpatialParameters", &Kernel::setSpatialParameters);
    clsKernel.def("computeKernelParametersFromSpatialModel",
                  &Kernel::computeKernelParametersFromSpatialModel);
    clsKernel.def("toString", &Kernel::toString, "prefix"_a = "");
    clsKernel.def("computeCache", &Kernel::computeCache);
    clsKernel.def("getCacheSize", &Kernel::getCacheSize);

    py::class_<FixedKernel, std::shared_ptr<FixedKernel>, Kernel> clsFixedKernel(mod, "FixedKernel");

    clsFixedKernel.def(py::init<>());
    clsFixedKernel.def(py::init<lsst::afw::image::Image<Kernel::Pixel> const &>(), "image"_a);
    clsFixedKernel.def(py::init<lsst::afw::math::Kernel const &, lsst::geom::Point2D const &>(), "kernel"_a,
                       "pos"_a);
    clsFixedKernel.def("clone", &FixedKernel::clone);
    clsFixedKernel.def("resized", &FixedKernel::resized, "width"_a, "height"_a);
    clsFixedKernel.def("toString", &FixedKernel::toString, "prefix"_a = "");
    clsFixedKernel.def("getSum", &FixedKernel::getSum);
    clsFixedKernel.def("isPersistable", &FixedKernel::isPersistable);

    py::class_<AnalyticKernel, std::shared_ptr<AnalyticKernel>, Kernel> clsAnalyticKernel(mod,
                                                                                          "AnalyticKernel");
    clsAnalyticKernel.def(py::init<>());
    // Workaround for NullSpatialFunction and py::arg not playing well with Citizen (TODO: no longer needed?)
    clsAnalyticKernel.def(py::init<int, int, AnalyticKernel::KernelFunction const &>(), "width"_a, "height"_a,
                          "kernelFunction"_a);
    clsAnalyticKernel.def(
            py::init<int, int, AnalyticKernel::KernelFunction const &, Kernel::SpatialFunction const &>(),
            "width"_a, "height"_a, "kernelFunction"_a, "spatialFunction"_a);
    clsAnalyticKernel.def(py::init<int, int, AnalyticKernel::KernelFunction const &,
                                   std::vector<Kernel::SpatialFunctionPtr> const &>(),
                          "width"_a, "height"_a, "kernelFunction"_a, "spatialFunctionList"_a);
    clsAnalyticKernel.def("clone", &AnalyticKernel::clone);
    clsAnalyticKernel.def("resized", &AnalyticKernel::resized, "width"_a, "height"_a);
    clsAnalyticKernel.def("computeImage", &AnalyticKernel::computeImage, "image"_a, "doNormalize"_a,
                          "x"_a = 0.0, "y"_a = 0.0);
    clsAnalyticKernel.def("getKernelParameters", &AnalyticKernel::getKernelParameters);
    clsAnalyticKernel.def("getKernelFunction", &AnalyticKernel::getKernelFunction);
    clsAnalyticKernel.def("toString", &AnalyticKernel::toString, "prefix"_a = "");
    clsAnalyticKernel.def("isPersistable", &AnalyticKernel::isPersistable);

    py::class_<DeltaFunctionKernel, std::shared_ptr<DeltaFunctionKernel>, Kernel> clsDeltaFunctionKernel(
            mod, "DeltaFunctionKernel");

    clsDeltaFunctionKernel.def(py::init<int, int, lsst::geom::Point2I const &>(), "width"_a, "height"_a,
                               "point"_a);
    clsDeltaFunctionKernel.def("clone", &DeltaFunctionKernel::clone);
    clsDeltaFunctionKernel.def("resized", &DeltaFunctionKernel::resized, "width"_a, "height"_a);
    clsDeltaFunctionKernel.def("getPixel", &DeltaFunctionKernel::getPixel);
    clsDeltaFunctionKernel.def("toString", &DeltaFunctionKernel::toString, "prefix"_a = "");
    clsDeltaFunctionKernel.def("isPersistable", &DeltaFunctionKernel::isPersistable);

    py::class_<LinearCombinationKernel, std::shared_ptr<LinearCombinationKernel>, Kernel>
            clsLinearCombinationKernel(mod, "LinearCombinationKernel");

    clsLinearCombinationKernel.def(py::init<>());
    clsLinearCombinationKernel.def(py::init<KernelList const &, std::vector<double> const &>(),
                                   "kernelList"_a, "kernelParameters"_a);
    clsLinearCombinationKernel.def(py::init<KernelList const &, Kernel::SpatialFunction const &>(),
                                   "kernelList"_a, "spatialFunction"_a);
    clsLinearCombinationKernel.def(
            py::init<KernelList const &, std::vector<Kernel::SpatialFunctionPtr> const &>(), "kernelList"_a,
            "spatialFunctionList"_a);
    clsLinearCombinationKernel.def("clone", &LinearCombinationKernel::clone);
    clsLinearCombinationKernel.def("resized", &LinearCombinationKernel::resized, "width"_a, "height"_a);
    clsLinearCombinationKernel.def("getKernelParameters", &LinearCombinationKernel::getKernelParameters);
    clsLinearCombinationKernel.def("getKernelList", &LinearCombinationKernel::getKernelList);
    clsLinearCombinationKernel.def("getKernelSumList", &LinearCombinationKernel::getKernelSumList);
    clsLinearCombinationKernel.def("getNBasisKernels", &LinearCombinationKernel::getNBasisKernels);
    clsLinearCombinationKernel.def("checkKernelList", &LinearCombinationKernel::checkKernelList);
    clsLinearCombinationKernel.def("isDeltaFunctionBasis", &LinearCombinationKernel::isDeltaFunctionBasis);
    clsLinearCombinationKernel.def("refactor", &LinearCombinationKernel::refactor);
    clsLinearCombinationKernel.def("toString", &LinearCombinationKernel::toString, "prefix"_a = "");
    clsLinearCombinationKernel.def("isPersistable", &LinearCombinationKernel::isPersistable);

    py::class_<SeparableKernel, std::shared_ptr<SeparableKernel>, Kernel> clsSeparableKernel(
            mod, "SeparableKernel");

    clsSeparableKernel.def(py::init<>());
    // Workaround for NullSpatialFunction and py::arg not playing well with Citizen (TODO: no longer needed?)
    clsSeparableKernel.def(py::init<int, int, SeparableKernel::KernelFunction const &,
                                    SeparableKernel::KernelFunction const &>(),
                           "width"_a, "height"_a, "kernelColFunction"_a, "kernelRowFunction"_a);
    clsSeparableKernel.def(
            py::init<int, int, SeparableKernel::KernelFunction const &,
                     SeparableKernel::KernelFunction const &, Kernel::SpatialFunction const &>(),
            "width"_a, "height"_a, "kernelColFunction"_a, "kernelRowFunction"_a, "spatialFunction"_a);
    clsSeparableKernel.def(py::init<int, int, SeparableKernel::KernelFunction const &,
                                    SeparableKernel::KernelFunction const &,
                                    std::vector<Kernel::SpatialFunctionPtr> const &>(),
                           "width"_a, "height"_a, "kernelColFunction"_a, "kernelRowFunction"_a,
                           "spatialFunctionList"_a);
    clsSeparableKernel.def("clone", &SeparableKernel::clone);
    clsSeparableKernel.def("resized", &SeparableKernel::resized, "width"_a, "height"_a);
    clsSeparableKernel.def("computeVectors", &SeparableKernel::computeVectors);
    clsSeparableKernel.def("getKernelParameter", &SeparableKernel::getKernelParameter);
    clsSeparableKernel.def("getKernelParameters", &SeparableKernel::getKernelParameters);
    clsSeparableKernel.def("getKernelColFunction", &SeparableKernel::getKernelColFunction);
    clsSeparableKernel.def("getKernelRowFunction", &SeparableKernel::getKernelRowFunction);
    clsSeparableKernel.def("toString", &SeparableKernel::toString, "prefix"_a = "");
    clsSeparableKernel.def("computeCache", &SeparableKernel::computeCache);
    clsSeparableKernel.def("getCacheSize", &SeparableKernel::getCacheSize);

}

WRAP(ConvolveImage) {
    py::class_<ConvolutionControl, std::shared_ptr<ConvolutionControl>> clsConvolutionControl(
            mod, "ConvolutionControl");

    clsConvolutionControl.def(py::init<bool, bool, int>(), "doNormalize"_a = true, "doCopyEdge"_a = false,
                              "maxInterpolationDistance"_a = 10);

    clsConvolutionControl.def("getDoNormalize", &ConvolutionControl::getDoNormalize);
    clsConvolutionControl.def("getDoCopyEdge", &ConvolutionControl::getDoCopyEdge);
    clsConvolutionControl.def("getMaxInterpolationDistance",
                              &ConvolutionControl::getMaxInterpolationDistance);
    clsConvolutionControl.def("setDoNormalize", &ConvolutionControl::setDoNormalize);
    clsConvolutionControl.def("setDoCopyEdge", &ConvolutionControl::setDoCopyEdge);
    clsConvolutionControl.def("setMaxInterpolationDistance",
                              &ConvolutionControl::setMaxInterpolationDistance);

    declareAll<double, double>(mod);
    declareAll<double, float>(mod);
    declareAll<double, int>(mod);
    declareAll<double, std::uint16_t>(mod);
    declareAll<float, float>(mod);
    declareAll<float, int>(mod);
    declareAll<float, std::uint16_t>(mod);
    declareAll<int, int>(mod);
    declareAll<std::uint16_t, std::uint16_t>(mod);
}

WRAP(Background) {

    /* Member types and enums */
    py::enum_<UndersampleStyle>(mod, "UndersampleStyle")
            .value("THROW_EXCEPTION", UndersampleStyle::THROW_EXCEPTION)
            .value("REDUCE_INTERP_ORDER", UndersampleStyle::REDUCE_INTERP_ORDER)
            .value("INCREASE_NXNYSAMPLE", UndersampleStyle::INCREASE_NXNYSAMPLE)
            .export_values();

    py::class_<BackgroundControl, std::shared_ptr<BackgroundControl>> clsBackgroundControl(
            mod, "BackgroundControl");

    /* Constructors */
    clsBackgroundControl.def(py::init<int const, int const, StatisticsControl const, Property const,
                                      ApproximateControl const>(),
                             "nxSample"_a, "nySample"_a, "sctrl"_a = StatisticsControl(), "prop"_a = MEANCLIP,
                             "actrl"_a = ApproximateControl(ApproximateControl::UNKNOWN, 1));
    clsBackgroundControl.def(py::init<int const, int const, StatisticsControl const, std::string const &,
                                      ApproximateControl const>(),
                             "nxSample"_a, "nySample"_a, "sctrl"_a, "prop"_a,
                             "actrl"_a = ApproximateControl(ApproximateControl::UNKNOWN, 1));
    clsBackgroundControl.def(py::init<Interpolate::Style const, int const, int const, UndersampleStyle const,
                                      StatisticsControl const, Property const, ApproximateControl const>(),
                             "style"_a, "nxSample"_a = 10, "nySample"_a = 10,
                             "undersampleStyle"_a = THROW_EXCEPTION, "sctrl"_a = StatisticsControl(),
                             "prop"_a = MEANCLIP,
                             "actrl"_a = ApproximateControl(ApproximateControl::UNKNOWN, 1));
    clsBackgroundControl.def(
            py::init<std::string const &, int const, int const, std::string const &, StatisticsControl const,
                     std::string const &, ApproximateControl const>(),
            "style"_a, "nxSample"_a = 10, "nySample"_a = 10, "undersampleStyle"_a = "THROW_EXCEPTION",
            "sctrl"_a = StatisticsControl(), "prop"_a = "MEANCLIP",
            "actrl"_a = ApproximateControl(ApproximateControl::UNKNOWN, 1));

    /* Members */
    clsBackgroundControl.def("setNxSample", &BackgroundControl::setNxSample);
    clsBackgroundControl.def("setNySample", &BackgroundControl::setNySample);
    clsBackgroundControl.def("setInterpStyle", (void (BackgroundControl::*)(Interpolate::Style const)) &
                                                       BackgroundControl::setInterpStyle);
    clsBackgroundControl.def("setInterpStyle", (void (BackgroundControl::*)(std::string const &)) &
                                                       BackgroundControl::setInterpStyle);
    clsBackgroundControl.def("setUndersampleStyle", (void (BackgroundControl::*)(UndersampleStyle const)) &
                                                            BackgroundControl::setUndersampleStyle);
    clsBackgroundControl.def("setUndersampleStyle", (void (BackgroundControl::*)(std::string const &)) &
                                                            BackgroundControl::setUndersampleStyle);
    clsBackgroundControl.def("getNxSample", &BackgroundControl::getNxSample);
    clsBackgroundControl.def("getNySample", &BackgroundControl::getNySample);
    clsBackgroundControl.def("getInterpStyle", &BackgroundControl::getInterpStyle);
    clsBackgroundControl.def("getUndersampleStyle", &BackgroundControl::getUndersampleStyle);
    clsBackgroundControl.def("getStatisticsControl",
                             (std::shared_ptr<StatisticsControl> (BackgroundControl::*)()) &
                                     BackgroundControl::getStatisticsControl);
    clsBackgroundControl.def("getStatisticsProperty", &BackgroundControl::getStatisticsProperty);
    clsBackgroundControl.def("setStatisticsProperty", (void (BackgroundControl::*)(Property)) &
                                                              BackgroundControl::setStatisticsProperty);
    clsBackgroundControl.def("setStatisticsProperty", (void (BackgroundControl::*)(std::string)) &
                                                              BackgroundControl::setStatisticsProperty);
    clsBackgroundControl.def("setApproximateControl", &BackgroundControl::setApproximateControl);
    clsBackgroundControl.def("getApproximateControl",
                             (std::shared_ptr<ApproximateControl> (BackgroundControl::*)()) &
                                     BackgroundControl::getApproximateControl);

    /* Note that, in this case, the holder type must be unique_ptr to enable usage
     * of py::nodelete, which in turn is needed because Background has a protected
     * destructor. Adding py::nodelete prevents pybind11 from calling the destructor
     * when the pointer is destroyed. Thus care needs to be taken to prevent leaks.
     * Basically Background should only ever be used as a base class (without data
     * members). */
    py::class_<Background, std::shared_ptr<Background> > clsBackground(mod, "Background");

    /* Members */
    declareGetImage<float>(clsBackground, "F");

    clsBackground.def("getAsUsedInterpStyle", &Background::getAsUsedInterpStyle);
    clsBackground.def("getAsUsedUndersampleStyle", &Background::getAsUsedUndersampleStyle);
    clsBackground.def("getApproximate", &Background::getApproximate, "actrl"_a,
                      "undersampleStyle"_a = THROW_EXCEPTION);
    clsBackground.def("getBackgroundControl", (std::shared_ptr<BackgroundControl> (Background::*)()) &
                                                      Background::getBackgroundControl);

    py::class_<BackgroundMI, std::shared_ptr<BackgroundMI>, Background> clsBackgroundMI(mod, "BackgroundMI");

    /* Constructors */
    clsBackgroundMI.def(
            py::init<lsst::geom::Box2I const, image::MaskedImage<typename Background::InternalPixelT> const &>(),
            "imageDimensions"_a, "statsImage"_a);

    /* Operators */
    clsBackgroundMI.def("__iadd__", &BackgroundMI::operator+=);
    clsBackgroundMI.def("__isub__", &BackgroundMI::operator-=);

    /* Members */
    clsBackgroundMI.def("getStatsImage", &BackgroundMI::getStatsImage);
    clsBackgroundMI.def("getImageBBox", &BackgroundMI::getImageBBox);

    // Yes, really only float
    declareMakeBackground<image::Image<float>>(mod);
    declareMakeBackground<image::MaskedImage<float>>(mod);

    mod.def("stringToUndersampleStyle", stringToUndersampleStyle, "style"_a);
}

WRAP(WarpExposure) {
    /* Module level */
    auto clsLanczosWarpingKernel = declareWarpingKernel<LanczosWarpingKernel>(mod, "LanczosWarpingKernel");
    declareSimpleWarpingKernel<BilinearWarpingKernel>(mod, "BilinearWarpingKernel");
    declareSimpleWarpingKernel<NearestWarpingKernel>(mod, "NearestWarpingKernel");

    py::class_<WarpingControl, std::shared_ptr<WarpingControl>> clsWarpingControl(mod, "WarpingControl");

    declareWarpingFunctions<double, double>(mod);
    declareWarpingFunctions<double, float>(mod);
    declareWarpingFunctions<double, int>(mod);
    declareWarpingFunctions<double, std::uint16_t>(mod);
    declareWarpingFunctions<float, float>(mod);
    declareWarpingFunctions<float, int>(mod);
    declareWarpingFunctions<float, std::uint16_t>(mod);
    declareWarpingFunctions<int, int>(mod);
    declareWarpingFunctions<std::uint16_t, std::uint16_t>(mod);

    /* Member types and enums */

    /* Constructors */
    clsLanczosWarpingKernel.def(py::init<int>(), "order"_a);

    clsWarpingControl.def(py::init<std::string, std::string, int, int, image::MaskPixel>(),
                          "warpingKernelName"_a, "maskWarpingKernelName"_a = "", "cacheSize"_a = 0,
                          "interpLength"_a = 0, "growFullMask"_a = 0);

    /* Operators */
    clsLanczosWarpingKernel.def("getOrder", &LanczosWarpingKernel::getOrder);

    clsWarpingControl.def("getCacheSize", &WarpingControl::getCacheSize);
    clsWarpingControl.def("setCacheSize", &WarpingControl::setCacheSize, "cacheSize"_a);
    clsWarpingControl.def("getInterpLength", &WarpingControl::getInterpLength);
    clsWarpingControl.def("setInterpLength", &WarpingControl::setInterpLength, "interpLength"_a);
    clsWarpingControl.def("setWarpingKernelName", &WarpingControl::setWarpingKernelName,
                          "warpingKernelName"_a);
    clsWarpingControl.def("getWarpingKernel", &WarpingControl::getWarpingKernel);
    clsWarpingControl.def("setWarpingKernel", &WarpingControl::setWarpingKernel, "warpingKernel"_a);
    clsWarpingControl.def("setMaskWarpingKernelName", &WarpingControl::setMaskWarpingKernelName,
                          "maskWarpingKernelName"_a);
    clsWarpingControl.def("getMaskWarpingKernel", &WarpingControl::getMaskWarpingKernel);
    clsWarpingControl.def("hasMaskWarpingKernel", &WarpingControl::hasMaskWarpingKernel);
    clsWarpingControl.def("setMaskWarpingKernelName", &WarpingControl::setMaskWarpingKernelName,
                          "maskWarpingKernelName"_a);
    clsWarpingControl.def("setMaskWarpingKernel", &WarpingControl::setMaskWarpingKernel,
                          "maskWarpingKernel"_a);
    clsWarpingControl.def("getGrowFullMask", &WarpingControl::getGrowFullMask);
    clsWarpingControl.def("setGrowFullMask", &WarpingControl::setGrowFullMask, "growFullMask"_a);

}

WRAP(Stack) {
    declareStatisticsStack<float>(mod);
    declareStatisticsStack<double>(mod);
}

template <typename Pixel>
void declareStatistics(py::module &mod) {
    mod.def("makeStatistics",
            (Statistics(*)(image::Image<Pixel> const &, image::Mask<image::MaskPixel> const &, int const,
                           StatisticsControl const &))makeStatistics<Pixel>,
            "img"_a, "msk"_a, "flags"_a, "sctrl"_a = StatisticsControl());
    mod.def("makeStatistics", (Statistics(*)(image::MaskedImage<Pixel> const &, int const,
                                             StatisticsControl const &))makeStatistics<Pixel>,
            "mimg"_a, "flags"_a, "sctrl"_a = StatisticsControl());
    mod.def("makeStatistics",
            (Statistics(*)(image::MaskedImage<Pixel> const &, image::Image<WeightPixel> const &, int const,
                           StatisticsControl const &))makeStatistics<Pixel>,
            "mimg"_a, "weights"_a, "flags"_a, "sctrl"_a = StatisticsControl());
    mod.def("makeStatistics",
            (Statistics(*)(image::Mask<image::MaskPixel> const &, int const, StatisticsControl const &))
                    makeStatistics,  // this is not a template, just a regular overload
            "msk"_a,
            "flags"_a, "sctrl"_a = StatisticsControl());
    mod.def("makeStatistics", (Statistics(*)(image::Image<Pixel> const &, int const,
                                             StatisticsControl const &))makeStatistics<Pixel>,
            "img"_a, "flags"_a, "sctrl"_a = StatisticsControl());
}

template <typename Pixel>
void declareStatisticsVectorOverloads(py::module &mod) {
    mod.def("makeStatistics", (Statistics(*)(std::vector<Pixel> const &, int const,
                                             StatisticsControl const &))makeStatistics<Pixel>,
            "v"_a, "flags"_a, "sctrl"_a = StatisticsControl());
    mod.def("makeStatistics", (Statistics(*)(std::vector<Pixel> const &, std::vector<WeightPixel> const &,
                                             int const, StatisticsControl const &))makeStatistics<Pixel>,
            "v"_a, "vweights"_a, "flags"_a, "sctrl"_a = StatisticsControl());
}

WRAP(Statistics) {
    py::enum_<Property>(mod, "Property", py::arithmetic())
            .value("NOTHING", Property::NOTHING)
            .value("ERRORS", Property::ERRORS)
            .value("NPOINT", Property::NPOINT)
            .value("MEAN", Property::MEAN)
            .value("STDEV", Property::STDEV)
            .value("VARIANCE", Property::VARIANCE)
            .value("MEDIAN", Property::MEDIAN)
            .value("IQRANGE", Property::IQRANGE)
            .value("MEANCLIP", Property::MEANCLIP)
            .value("STDEVCLIP", Property::STDEVCLIP)
            .value("VARIANCECLIP", Property::VARIANCECLIP)
            .value("MIN", Property::MIN)
            .value("MAX", Property::MAX)
            .value("SUM", Property::SUM)
            .value("MEANSQUARE", Property::MEANSQUARE)
            .value("ORMASK", Property::ORMASK)
            .value("NCLIPPED", Property::NCLIPPED)
            .value("NMASKED", Property::NMASKED)
            .export_values();

    mod.def("stringToStatisticsProperty", stringToStatisticsProperty);

    py::class_<StatisticsControl, std::shared_ptr<StatisticsControl>> clsStatisticsControl(
            mod, "StatisticsControl");

    py::enum_<StatisticsControl::WeightsBoolean>(clsStatisticsControl, "WeightsBoolean")
            .value("WEIGHTS_FALSE", StatisticsControl::WeightsBoolean::WEIGHTS_FALSE)
            .value("WEIGHTS_TRUE", StatisticsControl::WeightsBoolean::WEIGHTS_TRUE)
            .value("WEIGHTS_NONE", StatisticsControl::WeightsBoolean::WEIGHTS_NONE)
            .export_values();

    clsStatisticsControl.def(py::init<double, int, lsst::afw::image::MaskPixel, bool,
                                      typename StatisticsControl::WeightsBoolean>(),
                             "numSigmaClip"_a = 3.0, "numIter"_a = 3, "andMask"_a = 0x0, "isNanSafe"_a = true,
                             "useWeights"_a = StatisticsControl::WEIGHTS_NONE);

    clsStatisticsControl.def("getMaskPropagationThreshold", &StatisticsControl::getMaskPropagationThreshold);
    clsStatisticsControl.def("setMaskPropagationThreshold", &StatisticsControl::setMaskPropagationThreshold);
    clsStatisticsControl.def("getNumSigmaClip", &StatisticsControl::getNumSigmaClip);
    clsStatisticsControl.def("getNumIter", &StatisticsControl::getNumIter);
    clsStatisticsControl.def("getAndMask", &StatisticsControl::getAndMask);
    clsStatisticsControl.def("getNoGoodPixelsMask", &StatisticsControl::getNoGoodPixelsMask);
    clsStatisticsControl.def("getNanSafe", &StatisticsControl::getNanSafe);
    clsStatisticsControl.def("getWeighted", &StatisticsControl::getWeighted);
    clsStatisticsControl.def("getWeightedIsSet", &StatisticsControl::getWeightedIsSet);
    clsStatisticsControl.def("getCalcErrorFromInputVariance",
                             &StatisticsControl::getCalcErrorFromInputVariance);
    clsStatisticsControl.def("setNumSigmaClip", &StatisticsControl::setNumSigmaClip);
    clsStatisticsControl.def("setNumIter", &StatisticsControl::setNumIter);
    clsStatisticsControl.def("setAndMask", &StatisticsControl::setAndMask);
    clsStatisticsControl.def("setNoGoodPixelsMask", &StatisticsControl::setNoGoodPixelsMask);
    clsStatisticsControl.def("setNanSafe", &StatisticsControl::setNanSafe);
    clsStatisticsControl.def("setWeighted", &StatisticsControl::setWeighted);
    clsStatisticsControl.def("setCalcErrorFromInputVariance",
                             &StatisticsControl::setCalcErrorFromInputVariance);

    py::class_<Statistics> clsStatistics(mod, "Statistics");

    clsStatistics.def("getResult", &Statistics::getResult, "prop"_a = Property::NOTHING);
    clsStatistics.def("getError", &Statistics::getError, "prop"_a = Property::NOTHING);
    clsStatistics.def("getValue", &Statistics::getValue, "prop"_a = Property::NOTHING);
    clsStatistics.def("getOrMask", &Statistics::getOrMask);

    declareStatistics<unsigned short>(mod);
    declareStatistics<double>(mod);
    declareStatistics<float>(mod);
    declareStatistics<int>(mod);

    // Declare vector overloads separately to prevent casting errors
    // that otherwise (mysteriously) occur when overloads are tried
    // in order.
    declareStatisticsVectorOverloads<unsigned short>(mod);
    declareStatisticsVectorOverloads<double>(mod);
    declareStatisticsVectorOverloads<float>(mod);
    declareStatisticsVectorOverloads<int>(mod);
}

WRAP(Random) {
    py::class_<Random> clsRandom(mod, "Random");

    /* Member types and enums */
    py::enum_<Random::Algorithm>(clsRandom, "Algorithm")
            .value("MT19937", Random::Algorithm::MT19937)
            .value("RANLXS0", Random::Algorithm::RANLXS0)
            .value("RANLXS1", Random::Algorithm::RANLXS1)
            .value("RANLXS2", Random::Algorithm::RANLXS2)
            .value("RANLXD1", Random::Algorithm::RANLXD1)
            .value("RANLXD2", Random::Algorithm::RANLXD2)
            .value("RANLUX", Random::Algorithm::RANLUX)
            .value("RANLUX389", Random::Algorithm::RANLUX389)
            .value("CMRG", Random::Algorithm::CMRG)
            .value("MRG", Random::Algorithm::MRG)
            .value("TAUS", Random::Algorithm::TAUS)
            .value("TAUS2", Random::Algorithm::TAUS2)
            .value("GFSR4", Random::Algorithm::GFSR4)
            .value("NUM_ALGORITHMS", Random::Algorithm::NUM_ALGORITHMS)
            .export_values();

    /* Constructors */
    clsRandom.def(py::init<Random::Algorithm, unsigned long>(), "algorithm"_a = Random::Algorithm::MT19937,
                  "seed"_a = 1);
    clsRandom.def(py::init<std::string const &, unsigned long>(), "algorithm"_a, "seed"_a = 1);

    /* Members */
    clsRandom.def("deepCopy", &Random::deepCopy);
    clsRandom.def("getAlgorithm", &Random::getAlgorithm);
    clsRandom.def("getAlgorithmName", &Random::getAlgorithmName);
    clsRandom.def_static("getAlgorithmNames", &Random::getAlgorithmNames);
    clsRandom.def("getSeed", &Random::getSeed);
    clsRandom.def("uniform", &Random::uniform);
    clsRandom.def("uniformPos", &Random::uniformPos);
    clsRandom.def("uniformInt", &Random::uniformInt);
    clsRandom.def("flat", &Random::flat);
    clsRandom.def("gaussian", &Random::gaussian);
    clsRandom.def("chisq", &Random::chisq);
    clsRandom.def("poisson", &Random::poisson);

    // getState and setState are special, their std::string cannot
    // be converted to a Python string (needs to go to bytes instead)
    // thus use the same solution as employed with Swig
    clsRandom.def("getState", [](Random &self) -> py::object {
        std::string state = self.getState();
        return py::reinterpret_steal<py::object>(PyBytes_FromStringAndSize(state.data(), state.size()));
    });
    clsRandom.def("setState", [](Random &self, py::bytes const &state) { self.setState(state); });

    /* Module level */
    declareRandomImage<lsst::afw::image::Image<double>>(mod);
    declareRandomImage<lsst::afw::image::Image<float>>(mod);
}

WRAP(ProductBoundField) {
    py::class_<ProductBoundedField, std::shared_ptr<ProductBoundedField>, BoundedField> cls(mod, "ProductBoundedField");
    cls.def(py::init<std::vector<std::shared_ptr<BoundedField const>>>());
}


}  // namespace lsst::afw::math::<anonymous>

namespace detail {

namespace {
template <typename OutImageT, typename InImageT>
void declareByType(py::module &mod) {
    mod.def("basicConvolve",
            (void (*)(OutImageT &, InImageT const &, lsst::afw::math::Kernel const &,
                      lsst::afw::math::ConvolutionControl const &))basicConvolve<OutImageT, InImageT>);
    mod.def("basicConvolve",
            (void (*)(OutImageT &, InImageT const &, lsst::afw::math::DeltaFunctionKernel const &,
                      lsst::afw::math::ConvolutionControl const &))basicConvolve<OutImageT, InImageT>);
    mod.def("basicConvolve",
            (void (*)(OutImageT &, InImageT const &, lsst::afw::math::LinearCombinationKernel const &,
                      lsst::afw::math::ConvolutionControl const &))basicConvolve<OutImageT, InImageT>);
    mod.def("basicConvolve",
            (void (*)(OutImageT &, InImageT const &, lsst::afw::math::SeparableKernel const &,
                      lsst::afw::math::ConvolutionControl const &))basicConvolve<OutImageT, InImageT>);
    mod.def("convolveWithBruteForce",
            (void (*)(
                    OutImageT &, InImageT const &, lsst::afw::math::Kernel const &,
                    lsst::afw::math::ConvolutionControl const &))convolveWithBruteForce<OutImageT, InImageT>);
}
template <typename PixelType1, typename PixelType2>
void declareAll(py::module &mod) {
    using M1 = image::MaskedImage<PixelType1, image::MaskPixel, image::VariancePixel>;
    using M2 = image::MaskedImage<PixelType2, image::MaskPixel, image::VariancePixel>;

    declareByType<image::Image<PixelType1>, image::Image<PixelType2>>(mod);
    declareByType<M1, M2>(mod);
}
}

WRAP(Spline) {
    /* Module level */
    py::class_<Spline> clsSpline(mod, "Spline");
    py::class_<TautSpline, Spline> clsTautSpline(mod, "TautSpline");

    /* Member types and enums */
    py::enum_<TautSpline::Symmetry>(clsTautSpline, "Symmetry")
            .value("Unknown", TautSpline::Symmetry::Unknown)
            .value("Odd", TautSpline::Symmetry::Odd)
            .value("Even", TautSpline::Symmetry::Even)
            .export_values();

    /* Constructors */

    /* Operators */

    /* Members */
    clsSpline.def("interpolate", &Spline::interpolate);
    clsSpline.def("derivative", &Spline::derivative);

    clsTautSpline.def(py::init<std::vector<double> const&, std::vector<double> const&, double const,
                               TautSpline::Symmetry>(),
                      "x"_a, "y"_a, "gamma"_a = 0, "type"_a = lsst::afw::math::detail::TautSpline::Unknown);
    clsTautSpline.def("roots", &TautSpline::roots);

}

WRAP(Convolve) {
    declareAll<double, double>(mod);
    declareAll<double, float>(mod);
    declareAll<double, int>(mod);
    declareAll<double, std::uint16_t>(mod);
    declareAll<float, float>(mod);
    declareAll<float, int>(mod);
    declareAll<float, std::uint16_t>(mod);
    declareAll<int, int>(mod);
    declareAll<std::uint16_t, std::uint16_t>(mod);

    py::class_<KernelImagesForRegion, std::shared_ptr<KernelImagesForRegion>> clsKernelImagesForRegion(
            mod, "KernelImagesForRegion");

    py::enum_<KernelImagesForRegion::Location>(clsKernelImagesForRegion, "Location")
            .value("BOTTOM_LEFT", KernelImagesForRegion::Location::BOTTOM_LEFT)
            .value("BOTTOM_RIGHT", KernelImagesForRegion::Location::BOTTOM_RIGHT)
            .value("TOP_LEFT", KernelImagesForRegion::Location::TOP_LEFT)
            .value("TOP_RIGHT", KernelImagesForRegion::Location::TOP_RIGHT)
            .export_values();

    clsKernelImagesForRegion.def(
            py::init<KernelImagesForRegion::KernelConstPtr, lsst::geom::Box2I const &,
                     lsst::geom::Point2I const &, bool>(),
            "kernelPtr"_a, "bbox"_a, "xy0"_a, "doNormalize"_a);
    clsKernelImagesForRegion.def(
            py::init<KernelImagesForRegion::KernelConstPtr, lsst::geom::Box2I const &,
                     lsst::geom::Point2I const &, bool, KernelImagesForRegion::ImagePtr,
                     KernelImagesForRegion::ImagePtr, KernelImagesForRegion::ImagePtr,
                     KernelImagesForRegion::ImagePtr>(),
            "kernelPtr"_a, "bbox"_a, "xy0"_a, "doNormalize"_a, "bottomLeftImagePtr"_a,
            "bottomRightImagePtr"_a, "topLeftImagePtr"_a, "topRightImagePtr"_a);

    clsKernelImagesForRegion.def("getBBox", &KernelImagesForRegion::getBBox);
    clsKernelImagesForRegion.def("getXY0", &KernelImagesForRegion::getXY0);
    clsKernelImagesForRegion.def("getDoNormalize", &KernelImagesForRegion::getDoNormalize);
    clsKernelImagesForRegion.def("getImage", &KernelImagesForRegion::getImage);
    clsKernelImagesForRegion.def("getKernel", &KernelImagesForRegion::getKernel);
    clsKernelImagesForRegion.def("getPixelIndex", &KernelImagesForRegion::getPixelIndex);
    clsKernelImagesForRegion.def("computeNextRow", &KernelImagesForRegion::computeNextRow);
    clsKernelImagesForRegion.def_static("getMinInterpolationSize",
                                        KernelImagesForRegion::getMinInterpolationSize);

    py::class_<RowOfKernelImagesForRegion, std::shared_ptr<RowOfKernelImagesForRegion>>
            clsRowOfKernelImagesForRegion(mod, "RowOfKernelImagesForRegion");

    clsRowOfKernelImagesForRegion.def(py::init<int, int>(), "nx"_a, "ny"_a);

    clsRowOfKernelImagesForRegion.def("front", &RowOfKernelImagesForRegion::front);
    clsRowOfKernelImagesForRegion.def("back", &RowOfKernelImagesForRegion::back);
    clsRowOfKernelImagesForRegion.def("getNX", &RowOfKernelImagesForRegion::getNX);
    clsRowOfKernelImagesForRegion.def("getNY", &RowOfKernelImagesForRegion::getNY);
    clsRowOfKernelImagesForRegion.def("getYInd", &RowOfKernelImagesForRegion::getYInd);
    clsRowOfKernelImagesForRegion.def("getRegion", &RowOfKernelImagesForRegion::getRegion);
    clsRowOfKernelImagesForRegion.def("hasData", &RowOfKernelImagesForRegion::hasData);
    clsRowOfKernelImagesForRegion.def("isLastRow", &RowOfKernelImagesForRegion::isLastRow);
    clsRowOfKernelImagesForRegion.def("incrYInd", &RowOfKernelImagesForRegion::incrYInd);

}
} // namespace detail

WRAP(Math) {
    wrapSpatialCellCandidate(mod);
    wrapSpatialCellCandidateIterator(mod);
    wrapSpatialCell(mod);
    wrapSpatialCellSet(mod);
    wrapCandidateVisitor(mod);
    wrapSpatialCellImageCandidate(mod);

    wrapBoundedField(mod);
    wrapTransformBoundedField(mod);
    wrapPixelAreaBound(mod);
    wrapApproximate(mod);
    wrapChebyshexBoundedField(mod);
    wrapMinimize(mod);
    wrapGaussian(mod);
    wrapFunction(mod);
    wrapFunctionLibrary(mod);
    wrapLeastSquares(mod);
    wrapInterpolate(mod);
    wrapStatistics(mod);
    wrapOffsetImage(mod);
    wrapTestClasses(mod);
    //wrapSpatial(mod);
    wrapKernel(mod);
    wrapConvolveImage(mod);
    wrapBackground(mod);
    wrapWarpExposure(mod);
    wrapStack(mod);
    wrapRandom(mod);
    wrapProductBoundField(mod);
    detail::wrapSpline(mod);
    detail::wrapConvolve(mod);

}

}
}
}  // namespace lsst::afw::math
