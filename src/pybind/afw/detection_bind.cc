#include "pybind/afw_bind.h"

#include "lsst/afw/detection/Threshold.h"
#include <cstdint>

#include "ndarray/pybind11.h"

#include "lsst/utils/python.h"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/io/python.h"  // for addPersistableMethods
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/FootprintSet.h"
#include "lsst/afw/detection/FootprintMerge.h"
#include "lsst/afw/detection/FootprintCtrl.h"
#include "lsst/afw/table/io/python.h"  // for addPersistableMethods
#include "lsst/afw/detection/GaussianPsf.h"
#include "lsst/afw/detection/HeavyFootprint.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/table/BaseTable.h"
#include "lsst/afw/detection/Peak.h"
#include "lsst/afw/table/python/catalog.h"
#include "lsst/afw/table/python/columnView.h"
#include "lsst/utils/python/PySharedPtr.h"

#include "lsst/geom/Point.h"
#include "lsst/afw/image/Color.h"
#include "lsst/afw/table/io/python.h"  // for addPersistableMethods
#include "lsst/afw/typehandling/Storable.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/python.h"


namespace py = pybind11;
using namespace py::literals;
using lsst::utils::python::PySharedPtr;

namespace lsst {
namespace afw {
namespace detection {
using PyPeakRecord = py::class_<PeakRecord, std::shared_ptr<PeakRecord>, table::BaseRecord>;
using PyPeakTable = py::class_<PeakTable, std::shared_ptr<PeakTable>, table::BaseTable>;

namespace {
WRAP(Threshold) {
    py::class_<Threshold, std::shared_ptr<Threshold>> cls(mod, "Threshold");
    py::enum_<Threshold::ThresholdType>(cls, "ThresholdType")
            .value("VALUE", Threshold::ThresholdType::VALUE)
            .value("BITMASK", Threshold::ThresholdType::BITMASK)
            .value("STDEV", Threshold::ThresholdType::STDEV)
            .value("VARIANCE", Threshold::ThresholdType::VARIANCE)
            .value("PIXEL_STDEV", Threshold::ThresholdType::PIXEL_STDEV)
            .export_values();
    cls.def(py::init<double const, typename Threshold::ThresholdType const, bool const,
            double const>(),
            "value"_a, "type"_a = Threshold::VALUE, "polarity"_a = true,
            "includeMultiplier"_a = 1.0);

    cls.def("getType", &Threshold::getType);
    cls.def_static("parseTypeString", Threshold::parseTypeString);
    cls.def_static("getTypeString", Threshold::getTypeString);
    cls.def("getValue", (double (Threshold::*)(const double) const) & Threshold::getValue,
            "param"_a = -1);
    //
    //    template<typename ImageT>
    //    double getValue(ImageT const& image) const;
    //
    cls.def("getPolarity", &Threshold::getPolarity);
    cls.def("setPolarity", &Threshold::setPolarity);
    cls.def("getIncludeMultiplier", &Threshold::getIncludeMultiplier);
    cls.def("setIncludeMultiplier", &Threshold::setIncludeMultiplier);


    mod.def("createThreshold", createThreshold, "value"_a, "type"_a = "value", "polarity"_a = true);
}

template <typename MaskT>
void declareMaskFromFootprintList(py::module_ &mod) {
    auto maskSetter = [](lsst::afw::image::Mask<MaskT> *mask,
                         std::vector<std::shared_ptr<lsst::afw::detection::Footprint>> const &footprints,
                         MaskT const bitmask, bool doClip) {
        for (auto const &foot : footprints) {
            try {
                if (doClip) {
                    auto tmpSpan = foot->getSpans()->clippedTo(mask->getBBox());
                    tmpSpan->setMask(*mask, bitmask);
                } else {
                    foot->getSpans()->setMask(*mask, bitmask);
                }
            } catch (lsst::pex::exceptions::OutOfRangeError const &e) {
                throw LSST_EXCEPT(lsst::pex::exceptions::OutOfRangeError,
                                  "Bounds of a Footprint fall outside mask set doClip to force");
            }
        }
    };

    mod.def("setMaskFromFootprintList", std::move(maskSetter), "mask"_a, "footprints"_a, "bitmask"_a,
                "doClip"_a = true);

}


WRAP(Footprint) {
    //wrappers.wrapType(
            py::class_<Footprint, std::shared_ptr<Footprint>> cls(mod, "Footprint");
    //        [](auto &mod, auto &cls) {
                cls.def(py::init<std::shared_ptr<geom::SpanSet>, lsst::geom::Box2I const &>(), "inputSpans"_a,
                        "region"_a = lsst::geom::Box2I());

                cls.def(py::init<std::shared_ptr<geom::SpanSet>, afw::table::Schema const &,
                                 lsst::geom::Box2I const &>(),
                        "inputSpans"_a, "peakSchema"_a, "region"_a = lsst::geom::Box2I());
                cls.def(py::init<Footprint const &>());
                cls.def(py::init<>());

                table::io::python::addPersistableMethods<Footprint>(cls);

                cls.def("getSpans", &Footprint::getSpans);
                cls.def("setSpans", &Footprint::setSpans);
                cls.def("getPeaks", (PeakCatalog & (Footprint::*)()) & Footprint::getPeaks,
                        py::return_value_policy::reference_internal);
                cls.def("addPeak", &Footprint::addPeak);
                cls.def("sortPeaks", &Footprint::sortPeaks, "key"_a = afw::table::Key<float>());
                cls.def("setPeakSchema", &Footprint::setPeakSchema);
                cls.def("setPeakCatalog", &Footprint::setPeakCatalog, "otherPeaks"_a);
                cls.def("getArea", &Footprint::getArea);
                cls.def("getCentroid", &Footprint::getCentroid);
                cls.def("getShape", &Footprint::getShape);
                cls.def("shift", (void (Footprint::*)(int, int)) & Footprint::shift);
                cls.def("shift", (void (Footprint::*)(lsst::geom::ExtentI const &)) & Footprint::shift);
                cls.def("getBBox", &Footprint::getBBox);
                cls.def("getRegion", &Footprint::getRegion);
                cls.def("setRegion", &Footprint::setRegion);
                cls.def("clipTo", &Footprint::clipTo);
                cls.def("contains", &Footprint::contains);
                cls.def("transform",
                        (std::shared_ptr<Footprint>(Footprint::*)(std::shared_ptr<geom::SkyWcs>,
                                                                  std::shared_ptr<geom::SkyWcs>,
                                                                  lsst::geom::Box2I const &, bool) const) &
                                Footprint::transform,
                        "source"_a, "target"_a, "region"_a, "doClip"_a = true);
                cls.def("transform",
                        (std::shared_ptr<Footprint>(Footprint::*)(lsst::geom::LinearTransform const &,
                                                                  lsst::geom::Box2I const &, bool) const) &
                                Footprint::transform);
                cls.def("transform",
                        (std::shared_ptr<Footprint>(Footprint::*)(lsst::geom::AffineTransform const &,
                                                                  lsst::geom::Box2I const &, bool) const) &
                                Footprint::transform);
                cls.def("transform",
                        (std::shared_ptr<Footprint>(Footprint::*)(geom::TransformPoint2ToPoint2 const &,
                                                                  lsst::geom::Box2I const &, bool) const) &
                                Footprint::transform);
                cls.def("dilate", (void (Footprint::*)(int, geom::Stencil)) & Footprint::dilate, "r"_a,
                        "stencil"_a = geom::Stencil::CIRCLE);
                cls.def("dilate", (void (Footprint::*)(geom::SpanSet const &)) & Footprint::dilate);
                cls.def("erode", (void (Footprint::*)(int, geom::Stencil)) & Footprint::erode, "r"_a,
                        "stencil"_a = geom::Stencil::CIRCLE);
                cls.def("erode", (void (Footprint::*)(geom::SpanSet const &)) & Footprint::erode);
                cls.def("removeOrphanPeaks", &Footprint::removeOrphanPeaks);
                cls.def("isContiguous", &Footprint::isContiguous);
                cls.def("isHeavy", &Footprint::isHeavy);
                cls.def("assign", (Footprint & (Footprint::*)(Footprint const &)) & Footprint::operator=);

                cls.def("split", [](Footprint const &self) -> py::list {
                    /* This is a work around for pybind not properly
                     * handling converting a vector of unique pointers
                     * to python lists of shared pointers */
                    py::list l;
                    for (auto &ptr : self.split()) {
                        l.append(py::cast(std::shared_ptr<Footprint>(std::move(ptr))));
                    }
                    return l;
                });

                cls.def_property("spans", &Footprint::getSpans, &Footprint::setSpans);
                cls.def_property_readonly("peaks", (PeakCatalog & (Footprint::*)()) & Footprint::getPeaks,
                                          py::return_value_policy::reference_internal);

                cls.def("__contains__", [](Footprint const &self, lsst::geom::Point2I const &point) -> bool {
                    return self.contains(point);
                });
                cls.def("__eq__",
                        [](Footprint const &self, Footprint const &other) -> bool { return self == other; },
                        py::is_operator());
            //});

    declareMaskFromFootprintList<lsst::afw::image::MaskPixel>(mod);

    //wrappers.wrap([](auto &mod) {
    mod.def("mergeFootprints", &mergeFootprints);
    mod.def("footprintToBBoxList", &footprintToBBoxList);
    //});
}

template <typename PixelT, typename PyClass>
void declareMakeHeavy(PyClass &cls) {
    //    cls.def("makeHeavy", [](FootprintSet & self, image::MaskedImage<PixelT, image::MaskPixel> const&
    //    mimg) {
    //            return self.makeHeavy(mimg);
    //            });
    //    cls.def("makeHeavy", [](FootprintSet & self, image::MaskedImage<PixelT, image::MaskPixel> const&
    //    mimg, HeavyFootprintCtrl const* ctrl) {
    //            return self.makeHeavy(mimg, ctrl);
    //            });
    cls.def("makeHeavy",
            (void (FootprintSet::*)(image::MaskedImage<PixelT, image::MaskPixel> const &,
                                    HeavyFootprintCtrl const *)) &
                    FootprintSet::makeHeavy<PixelT, image::MaskPixel>,
            "mimg"_a, "ctrl"_a = nullptr);
}

template <typename PixelT, typename PyClass>
void declareSetMask(PyClass &cls) {
    cls.def("setMask",
            (void (FootprintSet::*)(image::Mask<PixelT> *, std::string const &)) &
                    FootprintSet::setMask<PixelT>,
            "mask"_a, "planeName"_a);
}

template <typename PixelT, typename PyClass>
void declareTemplatedMembers(PyClass &cls) {
    /* Constructors */
    cls.def(py::init<image::Image<PixelT> const &, Threshold const &, int const, bool const>(), "img"_a,
            "threshold"_a, "npixMin"_a = 1, "setPeaks"_a = true);
    cls.def(py::init<image::MaskedImage<PixelT, image::MaskPixel> const &, Threshold const &,
                     std::string const &, int const, bool const>(),
            "img"_a, "threshold"_a, "planeName"_a = "", "npixMin"_a = 1, "setPeaks"_a = true);

    /* Members */
    declareMakeHeavy<int>(cls);
    declareMakeHeavy<float>(cls);
    declareMakeHeavy<double>(cls);
    declareMakeHeavy<std::uint16_t>(cls);
    //    declareMakeHeavy<std::uint64_t>(cls);
    declareSetMask<image::MaskPixel>(cls);
}

WRAP(FootprintSet) {
    py::class_<FootprintSet, std::shared_ptr<FootprintSet>> cls(mod, "FootprintSet");
    declareTemplatedMembers<std::uint16_t>(cls);
    declareTemplatedMembers<int>(cls);
    declareTemplatedMembers<float>(cls);
    declareTemplatedMembers<double>(cls);

    cls.def(py::init<image::Mask<image::MaskPixel> const &, Threshold const &, int const>(),
            "img"_a, "threshold"_a, "npixMin"_a = 1);

    cls.def(py::init<lsst::geom::Box2I>(), "region"_a);
    cls.def(py::init<FootprintSet const &>(), "set"_a);
    cls.def(py::init<FootprintSet const &, int, FootprintControl const &>(), "set"_a, "rGrow"_a,
            "ctrl"_a);
    cls.def(py::init<FootprintSet const &, int, bool>(), "set"_a, "rGrow"_a, "isotropic"_a);
    cls.def(py::init<FootprintSet const &, FootprintSet const &, bool>(), "footprints1"_a,
            "footprints2"_a, "includePeaks"_a);

    cls.def("swap", &FootprintSet::swap);
    // setFootprints takes shared_ptr<FootprintList> and getFootprints returns it,
    // but pybind11 can't handle that type, so use a custom getter and setter
    cls.def("setFootprints", [](FootprintSet &self, FootprintSet::FootprintList footList) {
        self.setFootprints(std::make_shared<FootprintSet::FootprintList>(std::move(footList)));
    });
    cls.def("getFootprints", [](FootprintSet &self) { return *(self.getFootprints()); });
    cls.def("makeSources", &FootprintSet::makeSources);
    cls.def("setRegion", &FootprintSet::setRegion);
    cls.def("getRegion", &FootprintSet::getRegion);
    cls.def("insertIntoImage", &FootprintSet::insertIntoImage);
    cls.def("setMask", (void (FootprintSet::*)(image::Mask<lsst::afw::image::MaskPixel> *,
                                               std::string const &)) &
            FootprintSet::setMask<lsst::afw::image::MaskPixel>);
    cls.def("setMask",
            (void (FootprintSet::*)(std::shared_ptr<image::Mask<lsst::afw::image::MaskPixel>>,
                                    std::string const &)) &
            FootprintSet::setMask<lsst::afw::image::MaskPixel>);
    cls.def("merge", &FootprintSet::merge, "rhs"_a, "tGrow"_a = 0, "rGrow"_a = 0,
            "isotropic"_a = true);
}

WRAP(FootprintMerge) {
    py::class_<FootprintMergeList>(mod, "FootprintMergeList")
            .def(py::init<afw::table::Schema &, std::vector<std::string> const &,
                 afw::table::Schema const &>(),
                 "sourceSchema"_a, "filterList"_a, "initialPeakSchema"_a)
            .def(py::init<afw::table::Schema &, std::vector<std::string> const &>(), "sourceSchema"_a,
                 "filterList"_a)

            .def("getPeakSchema", &FootprintMergeList::getPeakSchema)
            .def("addCatalog", &FootprintMergeList::addCatalog, "sourceTable"_a, "inputCat"_a,
                 "filter"_a, "minNewPeakDist"_a = -1., "doMerge"_a = true, "maxSamePeakDist"_a = -1.)
            .def("clearCatalog", &FootprintMergeList::clearCatalog)
            .def("getFinalSources", &FootprintMergeList::getFinalSources, "outputCat"_a);

}

WRAP(FootprintCtrl) {
    py::class_<FootprintControl>(mod, "FootprintControl")
            .def(py::init<>())
            .def(py::init<bool, bool>(), "circular"_a, "isotropic"_a = false)
            .def(py::init<bool, bool, bool, bool>(), "left"_a, "right"_a, "up"_a, "down"_a)

            .def("growCircular", &FootprintControl::growCircular)
            .def("growIsotropic", &FootprintControl::growIsotropic)
            .def("growLeft", &FootprintControl::growLeft)
            .def("growRight", &FootprintControl::growRight)
            .def("growUp", &FootprintControl::growUp)
            .def("growDown", &FootprintControl::growDown)

            .def("isCircular", &FootprintControl::isCircular)
            .def("isIsotropic", &FootprintControl::isIsotropic)
            .def("isLeft", &FootprintControl::isLeft)
            .def("isRight", &FootprintControl::isRight)
            .def("isUp", &FootprintControl::isUp)
            .def("isDown", &FootprintControl::isDown);

    py::class_<HeavyFootprintCtrl> cls(mod, "HeavyFootprintCtrl");
    py::enum_<HeavyFootprintCtrl::ModifySource>(cls, "ModifySource")
            .value("NONE", HeavyFootprintCtrl::ModifySource::NONE)
            .value("SET", HeavyFootprintCtrl::ModifySource::SET)
            .export_values();

    cls.def(py::init<HeavyFootprintCtrl::ModifySource>(),
            "modifySource"_a = HeavyFootprintCtrl::ModifySource::NONE);

    cls.def("getModifySource", &HeavyFootprintCtrl::getModifySource);
    cls.def("setModifySource", &HeavyFootprintCtrl::setModifySource);
    cls.def("getImageVal", &HeavyFootprintCtrl::getImageVal);
    cls.def("setImageVal", &HeavyFootprintCtrl::setImageVal);
    cls.def("getMaskVal", &HeavyFootprintCtrl::getMaskVal);
    cls.def("setMaskVal", &HeavyFootprintCtrl::setMaskVal);
    cls.def("getVarianceVal", &HeavyFootprintCtrl::getVarianceVal);
    cls.def("setVarianceVal", &HeavyFootprintCtrl::setVarianceVal);

}

WRAP(GaussianPsf) {
    py::class_<GaussianPsf, std::shared_ptr<GaussianPsf>, Psf> cls(mod, "GaussianPsf");
    table::io::python::addPersistableMethods<GaussianPsf>(cls);

    cls.def(py::init<int, int, double>(), "width"_a, "height"_a, "sigma"_a);
    cls.def(py::init<lsst::geom::Extent2I const&, double>(), "dimensions"_a, "sigma"_a);

    cls.def("clone", &GaussianPsf::clone);
    cls.def("resized", &GaussianPsf::resized, "width"_a, "height"_a);
    cls.def("getDimensions", &GaussianPsf::getDimensions);
    cls.def("getSigma", &GaussianPsf::getSigma);
    cls.def("isPersistable", &GaussianPsf::isPersistable);
}

template <typename ImagePixelT, typename MaskPixelT = lsst::afw::image::MaskPixel,
          typename VariancePixelT = lsst::afw::image::VariancePixel>
void declareHeavyFootprint(py::module_ &mod, std::string const &suffix) {
    using Class = HeavyFootprint<ImagePixelT>;
    py::class_<Class, std::shared_ptr<Class>, Footprint>(mod, ("HeavyFootprint" + suffix).c_str())
            .def(py::init<Footprint const &,
                 lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> const&,
                 HeavyFootprintCtrl const *>(),
                 "foot"_a, "mimage"_a, "ctrl"_a = nullptr)
            .def(py::init<Footprint const &, HeavyFootprintCtrl const *>(), "foot"_a,
                 "ctrl"_a = nullptr)

            .def("isHeavy", &Class::isHeavy)
            .def("insert", (void (Class::*)(lsst::afw::image::MaskedImage<ImagePixelT> &) const) &
                 Class::insert)
            .def("insert",
                 (void (Class::*)(lsst::afw::image::Image<ImagePixelT> &) const) & Class::insert)
            .def("getImageArray",
                 (ndarray::Array<ImagePixelT, 1, 1>(Class::*)()) & Class::getImageArray)
            .def("getMaskArray", (ndarray::Array<MaskPixelT, 1, 1>(Class::*)()) & Class::getMaskArray)
            .def("getVarianceArray",
                 (ndarray::Array<VariancePixelT, 1, 1>(Class::*)()) & Class::getVarianceArray)
            .def("getMaskBitsSet", &Class::getMaskBitsSet)
            .def("dot", &Class::dot);

    mod.def("makeHeavyFootprint",
            (Class(*)(Footprint const &,
                      lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> const &,
                      HeavyFootprintCtrl const *))
            makeHeavyFootprint<ImagePixelT, MaskPixelT, VariancePixelT>,
            "foot"_a, "img"_a, "ctrl"_a = nullptr);

    mod.def("mergeHeavyFootprints", mergeHeavyFootprints<ImagePixelT, MaskPixelT, VariancePixelT>);

}

WRAP(HeavyFootprint) {
    declareHeavyFootprint<int>(mod, "I");
    declareHeavyFootprint<std::uint16_t>(mod, "U");
    declareHeavyFootprint<float>(mod, "F");
    declareHeavyFootprint<double>(mod, "D");
}

/**
@internal Declare constructors and member and static functions for a pybind11 PeakRecord
*/
void declarePeakRecord(PyPeakRecord &cls) {
    cls.def("getTable", &PeakRecord::getTable);
    cls.def_property_readonly("table", &PeakRecord::getTable);
    cls.def("getId", &PeakRecord::getId);
    cls.def("setId", &PeakRecord::setId);
    cls.def("getIx", &PeakRecord::getIx);
    cls.def("getIy", &PeakRecord::getIy);
    cls.def("setIx", &PeakRecord::setIx);
    cls.def("setIy", &PeakRecord::setIy);
    cls.def("getI", &PeakRecord::getI);
    cls.def("getCentroid", (lsst::geom::Point2I(PeakRecord::*)(bool) const) & PeakRecord::getCentroid);
    cls.def("getCentroid", (lsst::geom::Point2D(PeakRecord::*)() const) & PeakRecord::getCentroid);
    cls.def("getFx", &PeakRecord::getFx);
    cls.def("getFy", &PeakRecord::getFy);
    cls.def("setFx", &PeakRecord::setFx);
    cls.def("setFy", &PeakRecord::setFy);
    cls.def("getF", &PeakRecord::getF);
    cls.def("getPeakValue", &PeakRecord::getPeakValue);
    cls.def("setPeakValue", &PeakRecord::setPeakValue);
    utils::python::addOutputOp(cls, "__str__");
    utils::python::addOutputOp(cls, "__repr__");
}

/**
@internal Declare constructors and member and static functions for a pybind11 PeakTable
*/
void declarePeakTable(PyPeakTable &cls) {
    cls.def_static("make", &PeakTable::make, "schema"_a, "forceNew"_a = false);
    cls.def_static("makeMinimalSchema", &PeakTable::makeMinimalSchema);
    cls.def_static("checkSchema", &PeakTable::checkSchema, "schema"_a);
    cls.def("getIdFactory", (std::shared_ptr<table::IdFactory>(PeakTable::*)()) & PeakTable::getIdFactory);
    cls.def("setIdFactory", &PeakTable::setIdFactory, "factory"_a);
    cls.def_static("getIdKey", &PeakTable::getIdKey);
    cls.def_static("getIxKey", &PeakTable::getIxKey);
    cls.def_static("getIyKey", &PeakTable::getIyKey);
    cls.def_static("getFxKey", &PeakTable::getFxKey);
    cls.def_static("getFyKey", &PeakTable::getFyKey);
    cls.def_static("getPeakValueKey", &PeakTable::getPeakValueKey);
    cls.def("clone", &PeakTable::clone);
    cls.def("makeRecord", &PeakTable::makeRecord);
    cls.def("copyRecord", (std::shared_ptr<PeakRecord>(PeakTable::*)(afw::table::BaseRecord const &)) &
                                  PeakTable::copyRecord);
    cls.def("copyRecord", (std::shared_ptr<PeakRecord>(PeakTable::*)(afw::table::BaseRecord const &,
                                                                     afw::table::SchemaMapper const &)) &
                                  PeakTable::copyRecord);
}

WRAP(Peak) {

    PyPeakRecord clsPeakRecord(mod, "PeakRecord");
    declarePeakRecord(clsPeakRecord);
    PyPeakTable clsPeakTable(mod, "PeakTable");
    declarePeakTable(clsPeakTable);

    auto clsPeakColumnView = table::python::declareColumnView<PeakRecord>(mod, "Peak");
    auto clsPeakCatalog = table::python::declareCatalog<PeakRecord>(mod, "Peak");

    clsPeakRecord.attr("Table") = clsPeakTable;
    clsPeakRecord.attr("ColumnView") = clsPeakColumnView;
    clsPeakRecord.attr("Catalog") = clsPeakCatalog;
    clsPeakTable.attr("Record") = clsPeakRecord;
    clsPeakTable.attr("ColumnView") = clsPeakColumnView;
    clsPeakTable.attr("Catalog") = clsPeakCatalog;
    clsPeakCatalog.attr("Record") = clsPeakRecord;
    clsPeakCatalog.attr("Table") = clsPeakTable;
    clsPeakCatalog.attr("ColumnView") = clsPeakColumnView;
}

auto const NullPoint = lsst::geom::Point2D(std::numeric_limits<double>::quiet_NaN());

WRAP(Psf) {
    py::class_<Psf, PySharedPtr<Psf>, typehandling::Storable, PsfTrampoline<> > cls(mod, "Psf");
    py::enum_<Psf::ImageOwnerEnum>(cls, "ImageOwnerEnum")
        .value("COPY", Psf::ImageOwnerEnum::COPY)
        .value("INTERNAL", Psf::ImageOwnerEnum::INTERNAL)
        .export_values();

    table::io::python::addPersistableMethods<Psf>(cls);
    cls.def(py::init<bool, size_t>(), "isFixed"_a=false, "capacity"_a=100);  // Constructor for pure-Python subclasses
    cls.def("clone", &Psf::clone);
    cls.def("resized", &Psf::resized, "width"_a, "height"_a);
    cls.def("computeImage", &Psf::computeImage, "position"_a = NullPoint,
            "color"_a = image::Color(), "owner"_a = Psf::ImageOwnerEnum::COPY);
    cls.def("computeKernelImage", &Psf::computeKernelImage, "position"_a = NullPoint,
            "color"_a = image::Color(), "owner"_a = Psf::ImageOwnerEnum::COPY);
    cls.def("computePeak", &Psf::computePeak, "position"_a = NullPoint,
            "color"_a = image::Color());
    cls.def("computeApertureFlux", &Psf::computeApertureFlux, "radius"_a,
            "position"_a = NullPoint, "color"_a = image::Color());
    cls.def("computeShape", &Psf::computeShape, "position"_a = NullPoint,
            "color"_a = image::Color());
    cls.def("computeBBox", &Psf::computeBBox, "position"_a = NullPoint,
            "color"_a = image::Color());
    cls.def("getLocalKernel", &Psf::getLocalKernel, "position"_a = NullPoint,
            "color"_a = image::Color());
    cls.def("getAverageColor", &Psf::getAverageColor);
    cls.def("getAveragePosition", &Psf::getAveragePosition);
    cls.def_static("recenterKernelImage", &Psf::recenterKernelImage, "im"_a, "position"_a,
                   "warpAlgorithm"_a = "lanczos5", "warpBuffer"_a = 5);
    cls.def("getCacheCapacity", &Psf::getCacheCapacity);
    cls.def("setCacheCapacity", &Psf::setCacheCapacity);

}

} // namespace


WRAP(Detection){
    wrapPsf(mod);
    wrapFootprintCtrl(mod);
    wrapFootprint(mod);
    wrapThreshold(mod);
    wrapFootprintSet(mod);
    wrapFootprintMerge(mod);
    wrapPeak(mod);
    wrapGaussianPsf(mod);
    wrapHeavyFootprint(mod);
}

}
}
}
