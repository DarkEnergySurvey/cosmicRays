#include "pybind/meas_bind.h"
#include "lsst/meas/algorithms/PsfCandidate.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
#include "lsst/pex/config/python.h"  // for LSST_DECLARE_CONTROL_FIELD
#include "lsst/meas/algorithms/WarpedPsf.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/SingleGaussianPsf.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/PcaPsf.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/KernelPsf.h"
#include "lsst/geom/Box.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Interp.h"
#include "lsst/utils/python/PySharedPtr.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/meas/algorithms/python.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/DoubleGaussianPsf.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/CR.h"
#include "lsst/meas/algorithms/CoaddTransmissionCurve.h"
#include "lsst/geom/Box.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/CoaddBoundedField.h"



namespace lsst {
namespace meas {
namespace algorithms {
namespace {
template <typename PixelT>
void declarePsfCandidate(py::module& mod, std::string const& suffix) {
    using Class = PsfCandidate<PixelT>;

    py::class_<Class, std::shared_ptr<Class>, afw::math::SpatialCellImageCandidate> cls(
            mod, ("PsfCandidate" + suffix).c_str());

    cls.def(py::init<std::shared_ptr<afw::table::SourceRecord> const&,
                     std::shared_ptr<afw::image::Exposure<PixelT> const>>(),
            "source"_a, "parentExposure"_a);
    cls.def(py::init<std::shared_ptr<afw::table::SourceRecord> const&,
                     std::shared_ptr<afw::image::Exposure<PixelT> const>, double, double>(),
            "source"_a, "parentExposure"_a, "xCenter"_a, "yCenter"_a);

    /* SpatialCellCandidate.getCandidateRating is defined in Python.
     * Therefore we cannot override it from the C++ wrapper.
     * Instead we give it a temporary name here, and assign it to the
     * class from Python. */
    cls.def("_getCandidateRating", &Class::getCandidateRating);
    cls.def("getSource", &Class::getSource);
    cls.def("getAmplitude", &Class::getAmplitude);
    cls.def("setAmplitude", &Class::setAmplitude);
    cls.def("getVar", &Class::getVar);
    cls.def("setVar", &Class::setVar);
    cls.def("getMaskedImage", (std::shared_ptr<afw::image::MaskedImage<PixelT> const>(Class::*)() const) &
                                      Class::getMaskedImage);
    cls.def("getMaskedImage",
            (std::shared_ptr<afw::image::MaskedImage<PixelT> const>(Class::*)(int, int) const) &
                    Class::getMaskedImage,
            "width"_a, "height"_a);
    cls.def("getOffsetImage", &Class::getOffsetImage);
    cls.def_static("getBorderWidth", &Class::getBorderWidth);
    cls.def_static("setBorderWidth", &Class::setBorderWidth);
    cls.def_static("setPixelThreshold", &Class::setPixelThreshold);
    cls.def_static("getPixelThreshold", &Class::getPixelThreshold);
    cls.def_static("setMaskBlends", &Class::setMaskBlends);
    cls.def_static("getMaskBlends", &Class::getMaskBlends);

    mod.def("makePsfCandidate", makePsfCandidate<PixelT>, "source"_a, "image"_a);
}


template <typename PixelT>
static void declareFunctions(py::module &mod) {
    using MaskedImageT = afw::image::MaskedImage<PixelT, afw::image::MaskPixel, afw::image::VariancePixel>;

    mod.def("createKernelFromPsfCandidates", createKernelFromPsfCandidates<PixelT>, "psfCells"_a, "dims"_a,
            "xy0"_a, "nEigenComponents"_a, "spatialOrder"_a, "ksize"_a, "nStarPerCell"_a = -1,
            "constantWeight"_a = true, "border"_a = 3);
    mod.def("countPsfCandidates", countPsfCandidates<PixelT>, "psfCells"_a, "nStarPerCell"_a = -1);
    mod.def("fitSpatialKernelFromPsfCandidates",
            (std::pair<bool, double>(*)(afw::math::Kernel *, afw::math::SpatialCellSet const &, int const,
                                        double const, double const))fitSpatialKernelFromPsfCandidates<PixelT>,
            "kernel"_a, "psfCells"_a, "nStarPerCell"_a = -1, "tolerance"_a = 1e-5, "lambda"_a = 0.0);
    mod.def("fitSpatialKernelFromPsfCandidates",
            (std::pair<bool, double>(*)(afw::math::Kernel *, afw::math::SpatialCellSet const &, bool const,
                                        int const, double const,
                                        double const))fitSpatialKernelFromPsfCandidates<PixelT>,
            "kernel"_a, "psfCells"_a, "doNonLinearFit"_a, "nStarPerCell"_a = -1, "tolerance"_a = 1e-5,
            "lambda"_a = 0.0);
    mod.def("subtractPsf", subtractPsf<MaskedImageT>, "psf"_a, "data"_a, "x"_a, "y"_a,
            "psfFlux"_a = std::numeric_limits<double>::quiet_NaN());
    mod.def("fitKernelParamsToImage", fitKernelParamsToImage<MaskedImageT>, "kernel"_a, "image"_a, "pos"_a);
    mod.def("fitKernelToImage", fitKernelToImage<MaskedImageT>, "kernel"_a, "image"_a, "pos"_a);
}

template <typename PixelT>
void declareInterpolateOverDefects(py::module& mod) {
    mod.def("interpolateOverDefects",
            interpolateOverDefects<
                    afw::image::MaskedImage<PixelT, afw::image::MaskPixel, afw::image::VariancePixel>>,
            "image"_a, "psf"_a, "badList"_a, "fallBackValue"_a = 0.0, "useFallbackValueAtEdge"_a = false);
}

template <typename PixelT>
void declareFindCosmicRays(py::module& mod) {
    mod.def("findCosmicRays", &findCosmicRays<PixelT>, "image"_a, "psf"_a, "bkgd"_a,
            "policy"_a, "keep"_a = false);
}

WRAP(Cr) {
    declareFindCosmicRays<float>(mod);
}

WRAP(Interp) {
    py::class_<Defect, std::shared_ptr<Defect>, afw::image::DefectBase> clsDefect(mod, "Defect");

    py::enum_<Defect::DefectPosition>(clsDefect, "DefectPosition")
            .value("LEFT", Defect::DefectPosition::LEFT)
            .value("NEAR_LEFT", Defect::DefectPosition::NEAR_LEFT)
            .value("WIDE_LEFT", Defect::DefectPosition::WIDE_LEFT)
            .value("MIDDLE", Defect::DefectPosition::MIDDLE)
            .value("WIDE_NEAR_LEFT", Defect::DefectPosition::WIDE_NEAR_LEFT)
            .value("WIDE", Defect::DefectPosition::WIDE)
            .value("WIDE_NEAR_RIGHT", Defect::DefectPosition::WIDE_NEAR_RIGHT)
            .value("NEAR_RIGHT", Defect::DefectPosition::NEAR_RIGHT)
            .value("WIDE_RIGHT", Defect::DefectPosition::WIDE_RIGHT)
            .value("RIGHT", Defect::DefectPosition::RIGHT)
            .export_values();

    clsDefect.def(py::init<const geom::BoxI&>(), "bbox"_a = geom::BoxI());

    clsDefect.def("classify", &Defect::classify);
    clsDefect.def("getType", &Defect::getType);
    clsDefect.def("getPos", &Defect::getPos);

    declareInterpolateOverDefects<float>(mod);
}
WRAP(SpatialModelPsf) {
    declareFunctions<float>(mod);
}
WRAP(PsfCandidate) {
    declarePsfCandidate<float>(mod, "F");
}

WRAP(CoaddPsf) {
    /* CoaddPsfControl */
    py::class_<CoaddPsfControl, std::shared_ptr<CoaddPsfControl>> clsControl(mod, "CoaddPsfControl");
    clsControl.def(py::init<std::string, int>(), "warpingKernelName"_a = "lanczos3", "cacheSize"_a = 10000);
    LSST_DECLARE_CONTROL_FIELD(clsControl, CoaddPsfControl, warpingKernelName);
    LSST_DECLARE_CONTROL_FIELD(clsControl, CoaddPsfControl, cacheSize);

    /* CoaddPsf */
    py::class_<CoaddPsf, std::shared_ptr<CoaddPsf>, ImagePsf> clsCoaddPsf(mod, "CoaddPsf");
    afw::table::io::python::addPersistableMethods<CoaddPsf>(clsCoaddPsf);

    /* Constructors */
    clsCoaddPsf.def(py::init<afw::table::ExposureCatalog const &, afw::geom::SkyWcs const &,
                             std::string const &, std::string const &, int>(),
                    "catalog"_a, "coaddWcs"_a, "weightFieldName"_a = "weight",
                    "warpingKernelName"_a = "lanczos3", "cacheSize"_a = 10000);
    clsCoaddPsf.def(py::init<afw::table::ExposureCatalog const &, afw::geom::SkyWcs const &,
                             CoaddPsfControl const &, std::string const &>(),
                    "catalog"_a, "coaddWcs"_a, "ctrl"_a, "weightFieldName"_a = "weight");

    /* Members */
    clsCoaddPsf.def("clone", &CoaddPsf::clone);
    clsCoaddPsf.def("getAveragePosition", &CoaddPsf::getAveragePosition);
    clsCoaddPsf.def("getCoaddWcs", &CoaddPsf::getCoaddWcs);
    clsCoaddPsf.def("getComponentCount", &CoaddPsf::getComponentCount);
    clsCoaddPsf.def("getPsf", &CoaddPsf::getPsf);
    clsCoaddPsf.def("getWcs", &CoaddPsf::getWcs);
    clsCoaddPsf.def("getWeight", &CoaddPsf::getWeight);
    clsCoaddPsf.def("getId", &CoaddPsf::getId);
    clsCoaddPsf.def("getBBox", &CoaddPsf::getBBox);
    clsCoaddPsf.def("getValidPolygon", &CoaddPsf::getValidPolygon);
    clsCoaddPsf.def("isPersistable", &CoaddPsf::isPersistable);
}

WRAP(WarpedPsf) {
    py::class_<WarpedPsf, std::shared_ptr<WarpedPsf>, ImagePsf> clsWarpedPsf(mod, "WarpedPsf");

    /* Constructors */
    clsWarpedPsf.def(py::init<std::shared_ptr<afw::detection::Psf const>,
                              std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const>,
                              std::shared_ptr<afw::math::WarpingControl const>>(),
                     "undistortedPsf"_a, "distortion"_a, "control"_a);
    clsWarpedPsf.def(py::init<std::shared_ptr<afw::detection::Psf const>,
                              std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const>, std::string const &,
                              unsigned int>(),
                     "undistortedPsf"_a, "distortion"_a, "kernelName"_a = "lanczos3", "cache"_a = 10000);

    /* Members */
    clsWarpedPsf.def("getAveragePosition", &WarpedPsf::getAveragePosition);
    clsWarpedPsf.def("clone", &WarpedPsf::clone);
}

WRAP(SingleGaussianPsf) {
    py::class_<SingleGaussianPsf, std::shared_ptr<SingleGaussianPsf>, KernelPsf> clsSingleGaussianPsf(
            mod, "SingleGaussianPsf");
    afw::table::io::python::addPersistableMethods<SingleGaussianPsf>(clsSingleGaussianPsf);

    clsSingleGaussianPsf.def(py::init<int, int, double>(), "width"_a, "height"_a, "sigma"_a);

    clsSingleGaussianPsf.def("clone", &SingleGaussianPsf::clone);
    clsSingleGaussianPsf.def("resized", &SingleGaussianPsf::resized, "width"_a, "height"_a);
    clsSingleGaussianPsf.def("getSigma", &SingleGaussianPsf::getSigma);
}

WRAP(PcaPsf) {
    py::class_<PcaPsf, std::shared_ptr<PcaPsf>, KernelPsf> clsPcaPsf(mod, "PcaPsf");
    afw::table::io::python::addPersistableMethods<PcaPsf>(clsPcaPsf);

    clsPcaPsf.def(py::init<std::shared_ptr<afw::math::LinearCombinationKernel>, geom::Point2D const &>(),
                  "kernel"_a, "averagePosition"_a = geom::Point2D());

    clsPcaPsf.def("clone", &PcaPsf::clone);
    clsPcaPsf.def("getKernel", &PcaPsf::getKernel);
}

WRAP(KernelPsf) {
    py::class_<KernelPsf, std::shared_ptr<KernelPsf>, ImagePsf> clsKernelPsf(mod, "KernelPsf");
    afw::table::io::python::addPersistableMethods<KernelPsf>(clsKernelPsf);

    clsKernelPsf.def(py::init<afw::math::Kernel const &, geom::Point2D const &>(), "kernel"_a,
                     "averagePosition"_a = geom::Point2D());

    clsKernelPsf.def("getKernel", &KernelPsf::getKernel);
    clsKernelPsf.def("getAveragePosition", &KernelPsf::getAveragePosition);
    clsKernelPsf.def("clone", &KernelPsf::clone);
}
using lsst::utils::python::PySharedPtr;

WRAP(ImagePsf) {

    py::class_<ImagePsf, PySharedPtr<ImagePsf>, afw::detection::Psf, ImagePsfTrampoline<>> clsImagePsf(
            mod, "ImagePsf");
    clsImagePsf.def(py::init<bool>(), "init", "isFixed"_a = false);  // Ctor for pure python subclasses
    afw::table::io::python::addPersistableMethods<ImagePsf>(clsImagePsf);
}

WRAP(DoubleGaussianPsf) {
    py::class_<DoubleGaussianPsf, std::shared_ptr<DoubleGaussianPsf>, KernelPsf> clsDoubleGaussianPsf(
            mod, "DoubleGaussianPsf");
    afw::table::io::python::addPersistableMethods<DoubleGaussianPsf>(clsDoubleGaussianPsf);

    clsDoubleGaussianPsf.def(py::init<int, int, double, double, double>(), "width"_a, "height"_a, "sigma1"_a,
                             "sigma2"_a = 0.0, "b"_a = 0.0);

    clsDoubleGaussianPsf.def("clone", &DoubleGaussianPsf::clone);
    clsDoubleGaussianPsf.def("resized", &DoubleGaussianPsf::resized, "width"_a, "height"_a);
    clsDoubleGaussianPsf.def("getSigma1", &DoubleGaussianPsf::getSigma1);
    clsDoubleGaussianPsf.def("getSigma2", &DoubleGaussianPsf::getSigma2);
    clsDoubleGaussianPsf.def("getB", &DoubleGaussianPsf::getB);
}

WRAP(CoaddTransmissionCurve) {
    mod.def("makeCoaddTransmissionCurve", &makeCoaddTransmissionCurve, "coaddWcs"_a, "inputSensors"_a);
}

WRAP(CoaddBoundedField) {
    py::class_<CoaddBoundedFieldElement> clsCoaddBoundedFieldElement(mod, "CoaddBoundedFieldElement");

    clsCoaddBoundedFieldElement.def(
            py::init([](std::shared_ptr<afw::math::BoundedField> field,
                        std::shared_ptr<afw::geom::SkyWcs const> wcs, py::object polygon, double weight) {
                if (polygon.is(py::none())) {
                    return new CoaddBoundedFieldElement(field, wcs, nullptr, weight);
                } else {
                    auto pgon = py::cast<std::shared_ptr<afw::geom::polygon::Polygon const>>(polygon);
                    return new CoaddBoundedFieldElement(field, wcs, pgon, weight);
                }
            }),
            "field"_a, "wcs"_a, "validPolygon"_a, "weight"_a = 1.0);

    clsCoaddBoundedFieldElement.def_readwrite("field", &CoaddBoundedFieldElement::field);
    clsCoaddBoundedFieldElement.def_readwrite("wcs", &CoaddBoundedFieldElement::wcs);
    clsCoaddBoundedFieldElement.def_readwrite("validPolygon", &CoaddBoundedFieldElement::validPolygon);
    clsCoaddBoundedFieldElement.def_readwrite("weight", &CoaddBoundedFieldElement::weight);

    clsCoaddBoundedFieldElement.def("__eq__", &CoaddBoundedFieldElement::operator==, py::is_operator());
    clsCoaddBoundedFieldElement.def("__ne__", &CoaddBoundedFieldElement::operator!=, py::is_operator());

    py::class_<CoaddBoundedField, std::shared_ptr<CoaddBoundedField>, afw::math::BoundedField>
            clsCoaddBoundedField(mod, "CoaddBoundedField");
    afw::table::io::python::addPersistableMethods<CoaddBoundedField>(clsCoaddBoundedField);

    clsCoaddBoundedField.attr("Element") = clsCoaddBoundedFieldElement;

    /* Constructors */
    clsCoaddBoundedField.def(py::init<geom::Box2I const &, std::shared_ptr<afw::geom::SkyWcs const>,
                                      typename CoaddBoundedField::ElementVector const &>(),
                             "bbox"_a, "coaddWcs"_a, "elements"_a);
    clsCoaddBoundedField.def(py::init<geom::Box2I const &, std::shared_ptr<afw::geom::SkyWcs const>,
                                      typename CoaddBoundedField::ElementVector const &, double>(),
                             "bbox"_a, "coaddWcs"_a, "elements"_a, "default"_a);

    /* Operators */
    clsCoaddBoundedField.def("__eq__", &CoaddBoundedField::operator==, py::is_operator());
    clsCoaddBoundedField.def("__ne__", &CoaddBoundedField::operator!=, py::is_operator());
    clsCoaddBoundedField.def("__imul__", &CoaddBoundedField::operator*);

    /* Members */
    clsCoaddBoundedField.def("evaluate", &CoaddBoundedField::evaluate);
    clsCoaddBoundedField.def("getCoaddWcs", &CoaddBoundedField::getCoaddWcs);
    clsCoaddBoundedField.def("getDefault", &CoaddBoundedField::getDefault);
    clsCoaddBoundedField.def("getElements", &CoaddBoundedField::getElements);
}
}

WRAP(Alg) {
    wrapPsfCandidate(mod);
    wrapSpatialModelPsf(mod);
    wrapImagePsf(mod);
    wrapKernelPsf(mod);
    wrapPcaPsf(mod);
    wrapSingleGaussianPsf(mod);
    wrapInterp(mod);
    wrapCoaddPsf(mod);
    wrapWarpedPsf(mod);
    wrapDoubleGaussianPsf(mod);
    wrapCr(mod);
    wrapCoaddTransmissionCurve(mod);
    wrapCoaddBoundedField(mod);
}
}
}
}
