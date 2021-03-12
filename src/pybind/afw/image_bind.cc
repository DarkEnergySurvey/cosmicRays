#include "pybind/afw_bind.h"
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "ndarray/pybind11.h"

#include "lsst/afw/cameraGeom/Detector.h"
#include "lsst/afw/coord/Observatory.h"
#include "lsst/afw/coord/Weather.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/fits.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/geom/polygon/Polygon.h"
#include "lsst/afw/image/ApCorrMap.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/CoaddInputs.h"
#include "lsst/afw/image/Color.h"
#include "lsst/afw/image/Defect.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/ExposureFitsReader.h"
#include "lsst/afw/image/ExposureInfo.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/afw/image/FilterLabel.h"
//#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/ImageBaseFitsReader.h"
#include "lsst/afw/image/ImageFitsReader.h"
#include "lsst/afw/image/ImagePca.h"
#include "lsst/afw/image/ImageSlice.h"
#include "lsst/afw/image/ImageUtils.h"
//#include "lsst/afw/image/Mask.h"
#include "lsst/afw/image/MaskFitsReader.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/MaskedImageFitsReader.h"
#include "lsst/afw/image/PhotoCalib.h"
//#include "lsst/afw/image/Pixel.h"
#include "lsst/afw/image/TransmissionCurve.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/afw/image/python/indexing.h"
#include "lsst/afw/math/BoundedField.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/Schema.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/afw/table/misc.h"
#include "lsst/afw/typehandling/Storable.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/SpherePoint.h"
#include "lsst/utils/python.h"
#include "lsst/utils/python/TemplateInvoker.h"


namespace py = pybind11;
using namespace pybind11::literals;
using namespace std::string_literals;

namespace lsst {
namespace afw {
namespace image {

namespace {

template <typename ImagePixelT>  // only the image type varies; mask and variance are fixed
using PyMaskedImage = py::class_<MaskedImage<ImagePixelT>, std::shared_ptr<MaskedImage<ImagePixelT>>>;

/**
@internal Declare a constructor that takes a MaskedImage of FromPixelT and returns a MaskedImage cast to
ToPixelT

The mask and variance must be of the standard types.

@param[in] cls  The pybind11 class to which add the constructor
*/
template <typename FromPixelT, typename ToPixelT>
void declareCastConstructor(PyMaskedImage<ToPixelT> &cls) {
    cls.def(py::init<MaskedImage<FromPixelT> const &, bool const>(), "src"_a, "deep"_a);
}

template <typename ImagePixelT>
PyMaskedImage<ImagePixelT> declareMaskedImage(py::module &mod, const std::string &suffix) {
    using MI = MaskedImage<ImagePixelT>;

    PyMaskedImage<ImagePixelT> cls(mod, ("MaskedImage" + suffix).c_str());

    mod.def("makeMaskedImage", &makeMaskedImage<ImagePixelT, MaskPixel, VariancePixel>, "image"_a,
            "mask"_a = nullptr, "variance"_a = nullptr);

    /* Member types and enums */

    /* Constructors */
    cls.def(py::init<unsigned int, unsigned int, typename MI::MaskPlaneDict const &>(), "width"_a, "height"_a,
            "planeDict"_a = typename MI::MaskPlaneDict());
    cls.def(py::init<lsst::geom::Extent2I, typename MI::MaskPlaneDict const &>(), "dimensions"_a,
            "planeDict"_a = typename MI::MaskPlaneDict());
    cls.def(py::init<typename MI::ImagePtr, typename MI::MaskPtr, typename MI::VariancePtr>(), "image"_a,
            "mask"_a = nullptr, "variance"_a = nullptr);
    cls.def(py::init<lsst::geom::Box2I const &, typename MI::MaskPlaneDict const &>(), "bbox"_a,
            "planeDict"_a = typename MI::MaskPlaneDict());
    cls.def(py::init<std::string const &, std::shared_ptr<daf::base::PropertySet>, lsst::geom::Box2I const &,
                     ImageOrigin, bool, bool, std::shared_ptr<daf::base::PropertySet>,
                     std::shared_ptr<daf::base::PropertySet>, std::shared_ptr<daf::base::PropertySet>,
                     bool>(),
            "fileName"_a, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT,
            "conformMasks"_a = false, "needAllHdus"_a = false, "imageMetadata"_a = nullptr,
            "maskMetadata"_a = nullptr, "varianceMetadata"_a = nullptr, "allowUnsafe"_a = false);
    cls.def(py::init<fits::MemFileManager &, std::shared_ptr<daf::base::PropertySet>,
                     lsst::geom::Box2I const &, ImageOrigin, bool, bool,
                     std::shared_ptr<daf::base::PropertySet>, std::shared_ptr<daf::base::PropertySet>,
                     std::shared_ptr<daf::base::PropertySet>, bool>(),
            "manager"_a, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT,
            "conformMasks"_a = false, "needAllHdus"_a = false, "imageMetadata"_a = nullptr,
            "maskMetadata"_a = nullptr, "varianceMetadata"_a = nullptr, "allowUnsafe"_a = false);
    cls.def(py::init<MI const &, bool>(), "rhs"_a, "deep"_a = false);
    cls.def(py::init<MI const &, lsst::geom::Box2I const &, ImageOrigin, bool>(), "rhs"_a, "bbox"_a,
            "origin"_a = PARENT, "deep"_a = false);

    /* Operators */
    cls.def("swap", &MI::swap);
    cls.def("assign", &MI::assign, "rhs"_a, "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT,
            py::is_operator()  // py::is_operator is a workaround for code in slicing.py
                               // that expects NotImplemented to be returned on failure.
    );

    cls.def("subset", &MI::subset, "bbox"_a, "origin"_a = PARENT);

    cls.def("__iadd__", (MI & (MI::*)(ImagePixelT const)) & MI::operator+=);
    cls.def("__iadd__", (MI & (MI::*)(MI const &)) & MI::operator+=);
    cls.def("__iadd__", (MI & (MI::*)(Image<ImagePixelT> const &)) & MI::operator+=);
    cls.def("__iadd__", (MI & (MI::*)(math::Function2<double> const &)) & MI::operator+=);
    cls.def("scaledPlus", &MI::scaledPlus);
    cls.def("__isub__", (MI & (MI::*)(ImagePixelT const)) & MI::operator-=);
    cls.def("__isub__", (MI & (MI::*)(MI const &)) & MI::operator-=);
    cls.def("__isub__", (MI & (MI::*)(Image<ImagePixelT> const &)) & MI::operator-=);
    cls.def("__isub__", (MI & (MI::*)(math::Function2<double> const &)) & MI::operator-=);
    cls.def("scaledMinus", &MI::scaledMinus);
    cls.def("__imul__", (MI & (MI::*)(ImagePixelT const)) & MI::operator*=);
    cls.def("__imul__", (MI & (MI::*)(MI const &)) & MI::operator*=);
    cls.def("__imul__", (MI & (MI::*)(Image<ImagePixelT> const &)) & MI::operator*=);
    cls.def("scaledMultiplies", &MI::scaledMultiplies);
    cls.def("__itruediv__", (MI & (MI::*)(ImagePixelT const)) & MI::operator/=);
    cls.def("__itruediv__", (MI & (MI::*)(MI const &)) & MI::operator/=);
    cls.def("__itruediv__", (MI & (MI::*)(Image<ImagePixelT> const &)) & MI::operator/=);
    cls.def("scaledDivides", &MI::scaledDivides);

    /* Members */
    cls.def("writeFits",
            (void (MI::*)(std::string const &, std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>) const) &
                    MI::writeFits,
            "fileName"_a, "metadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "imageMetadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "maskMetadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "varianceMetadata"_a = std::shared_ptr<daf::base::PropertySet const>());
    cls.def("writeFits",
            (void (MI::*)(fits::MemFileManager &, std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>) const) &
                    MI::writeFits,
            "manager"_a, "metadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "imageMetadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "maskMetadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "varianceMetadata"_a = std::shared_ptr<daf::base::PropertySet const>());
    cls.def("writeFits",
            (void (MI::*)(fits::Fits &, std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>,
                          std::shared_ptr<daf::base::PropertySet const>) const) &
                    MI::writeFits,
            "fitsfile"_a, "metadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "imageMetadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "maskMetadata"_a = std::shared_ptr<daf::base::PropertySet const>(),
            "varianceMetadata"_a = std::shared_ptr<daf::base::PropertySet const>());

    cls.def("writeFits",
            [](MI &self, std::string const &filename, fits::ImageWriteOptions const &imageOptions,
               fits::ImageWriteOptions const &maskOptions, fits::ImageWriteOptions const &varianceOptions,
               std::shared_ptr<daf::base::PropertySet const> header) {
                self.writeFits(filename, imageOptions, maskOptions, varianceOptions, header);
            },
            "filename"_a, "imageOptions"_a, "maskOptions"_a, "varianceOptions"_a,
            "header"_a = std::shared_ptr<daf::base::PropertyList>());
    cls.def("writeFits",
            [](MI &self, fits::MemFileManager &manager, fits::ImageWriteOptions const &imageOptions,
               fits::ImageWriteOptions const &maskOptions, fits::ImageWriteOptions const &varianceOptions,
               std::shared_ptr<daf::base::PropertySet const> header) {
                self.writeFits(manager, imageOptions, maskOptions, varianceOptions, header);
            },
            "manager"_a, "imageOptions"_a, "maskOptions"_a, "varianceOptions"_a,
            "header"_a = std::shared_ptr<daf::base::PropertyList>());
    cls.def("writeFits",
            [](MI &self, fits::Fits &fits, fits::ImageWriteOptions const &imageOptions,
               fits::ImageWriteOptions const &maskOptions, fits::ImageWriteOptions const &varianceOptions,
               std::shared_ptr<daf::base::PropertySet const> header) {
                self.writeFits(fits, imageOptions, maskOptions, varianceOptions, header);
            },
            "fits"_a, "imageOptions"_a, "maskOptions"_a, "varianceOptions"_a,
            "header"_a = std::shared_ptr<daf::base::PropertyList>());

    cls.def_static("readFits", (MI(*)(std::string const &))MI::readFits, "filename"_a);
    cls.def_static("readFits", (MI(*)(fits::MemFileManager &))MI::readFits, "manager"_a);
    cls.def("getImage", &MI::getImage);
    cls.def("setImage", &MI::setImage);
    cls.def_property("image", &MI::getImage, &MI::setImage);
    cls.def("getMask", &MI::getMask);
    cls.def("setMask", &MI::setMask);
    cls.def_property("mask", &MI::getMask, &MI::setMask);
    cls.def("getVariance", &MI::getVariance);
    cls.def("setVariance", &MI::setVariance);
    cls.def_property("variance", &MI::getVariance, &MI::setVariance);
    cls.def("getWidth", &MI::getWidth);
    cls.def("getHeight", &MI::getHeight);
    cls.def("getDimensions", &MI::getDimensions);
    cls.def("getBBox", &MI::getBBox, "origin"_a = PARENT);
    cls.def("getX0", &MI::getX0);
    cls.def("getY0", &MI::getY0);
    cls.def("getXY0", &MI::getXY0);
    cls.def("setXY0", (void (MI::*)(int const, int const)) & MI::setXY0, "x0"_a, "y0"_a);
    cls.def("setXY0", (void (MI::*)(lsst::geom::Point2I const)) & MI::setXY0, "origin"_a);
    cls.def("indexToPosition", &MI::indexToPosition);
    cls.def("positionToIndex", &MI::positionToIndex);

    return cls;
}

template <typename ImagePixelT>  // addtional template types do not seem to be needed
void declareMakeMaskedImage(py::module &mod) {
    mod.def("makeMaskedImage", makeMaskedImage<ImagePixelT, MaskPixel, VariancePixel>, "image"_a,
            "mask"_a = nullptr, "variance"_a = nullptr);
}

template <typename ImagePixelT1, typename ImagePixelT2>
void declareImagesOverlap(py::module &mod) {
    // wrap both the Image and MaskedImage versions of imagesOverlap here, as wrapping
    // the Image version in the Image wrapper results in it being invisible in lsst.afw.image
    mod.def("imagesOverlap",
            py::overload_cast<ImageBase<ImagePixelT1> const &, ImageBase<ImagePixelT2> const &>(
                    &imagesOverlap<ImagePixelT1, ImagePixelT2>),
            "image1"_a, "image2"_a);

    mod.def("imagesOverlap",
            py::overload_cast<MaskedImage<ImagePixelT1> const &, MaskedImage<ImagePixelT2> const &>(
                    &imagesOverlap<ImagePixelT1, ImagePixelT2>),
            "image1"_a, "image2"_a);
}

template <typename PixelT>
using PyImageBase = py::class_<ImageBase<PixelT>, std::shared_ptr<ImageBase<PixelT>>>;

template <typename PixelT>
using PyImage = py::class_<Image<PixelT>, std::shared_ptr<Image<PixelT>>, ImageBase<PixelT>>;

template <typename PixelT>
using PyDecoratedImage = py::class_<DecoratedImage<PixelT>, std::shared_ptr<DecoratedImage<PixelT>>>;

template <typename MaskPixelT>
using PyMask = py::class_<Mask<MaskPixelT>, std::shared_ptr<Mask<MaskPixelT>>, ImageBase<MaskPixelT>>;

/**
@internal Declare a constructor that takes a MaskedImage of FromPixelT and returns a MaskedImage cast to
ToPixelT

The mask and variance must be of the standard types.

@param[in] cls  The pybind11 class to which add the constructor
*/
template <typename FromPixelT, typename ToPixelT>
static void declareCastConstructor(PyImage<ToPixelT> &cls) {
    cls.def(py::init<Image<FromPixelT> const &, bool const>(), "src"_a, "deep"_a);
}

template <typename PixelT>
static void declareImageBase(py::module &mod, std::string const &suffix) {
    PyImageBase<PixelT> cls(mod, ("ImageBase" + suffix).c_str());

    using Array = typename ImageBase<PixelT>::Array;

    cls.def(py::init<lsst::geom::Extent2I const &>(), "dimensions"_a = lsst::geom::Extent2I());
    cls.def(py::init<ImageBase<PixelT> const &, bool>(), "src"_a, "deep"_a = false);
    cls.def(py::init<ImageBase<PixelT> const &, lsst::geom::Box2I const &, ImageOrigin, bool>(), "src"_a,
            "bbox"_a, "origin"_a = PARENT, "deep"_a = false);
    cls.def(py::init<Array const &, bool, lsst::geom::Point2I const &>(), "array"_a, "deep"_a = false,
            "xy0"_a = lsst::geom::Point2I());

    cls.def("assign", &ImageBase<PixelT>::assign, "rhs"_a, "bbox"_a = lsst::geom::Box2I(),
            "origin"_a = PARENT,
            py::is_operator());  // py::is_operator is a workaround for code in slicing.py
                                 // that expects NotImplemented to be returned on failure.
    cls.def("getWidth", &ImageBase<PixelT>::getWidth);
    cls.def("getHeight", &ImageBase<PixelT>::getHeight);
    cls.def("getX0", &ImageBase<PixelT>::getX0);
    cls.def("getY0", &ImageBase<PixelT>::getY0);
    cls.def("getXY0", &ImageBase<PixelT>::getXY0);
    cls.def("positionToIndex", &ImageBase<PixelT>::positionToIndex, "position"_a, "xOrY"_a);
    cls.def("indexToPosition", &ImageBase<PixelT>::indexToPosition, "index"_a, "xOrY"_a);
    cls.def("getDimensions", &ImageBase<PixelT>::getDimensions);
    cls.def("getArray", (Array(ImageBase<PixelT>::*)()) & ImageBase<PixelT>::getArray);
    cls.def_property("array", (Array(ImageBase<PixelT>::*)()) & ImageBase<PixelT>::getArray,
                     [](ImageBase<PixelT> &self, ndarray::Array<PixelT const, 2, 0> const &array) {
                         // Avoid self-assignment, which is invoked when a Python in-place operator is used.
                         if (array.shallow() != self.getArray().shallow()) {
                             self.getArray().deep() = array;
                         }
                     });
    cls.def("setXY0", (void (ImageBase<PixelT>::*)(lsst::geom::Point2I const)) & ImageBase<PixelT>::setXY0,
            "xy0"_a);
    cls.def("setXY0", (void (ImageBase<PixelT>::*)(int const, int const)) & ImageBase<PixelT>::setXY0, "x0"_a,
            "y0"_a);
    cls.def("getBBox", &ImageBase<PixelT>::getBBox, "origin"_a = PARENT);

    cls.def("set", [](ImageBase<PixelT> &img, PixelT val) { img = val; });
    cls.def("_set",
            [](ImageBase<PixelT> &img, lsst::geom::Point2I const &index, PixelT val, ImageOrigin origin) {
                python::checkBounds(index, img.getBBox(origin));
                img.get(index, origin) = val;
            },
            "index"_a, "value"_a, "origin"_a);
    cls.def("_get",
            [](ImageBase<PixelT> &img, lsst::geom::Point2I const &index, ImageOrigin origin) {
                python::checkBounds(index, img.getBBox(origin));
                return img.get(index, origin);
            },
            "index"_a, "origin"_a);
}

template <typename MaskPixelT>
static void declareMask(py::module &mod, std::string const &suffix) {
    PyMask<MaskPixelT> cls(mod, ("Mask" + suffix).c_str());

    /* Constructors */
    cls.def(py::init<unsigned int, unsigned int, typename Mask<MaskPixelT>::MaskPlaneDict const &>(),
            "width"_a, "height"_a, "planeDefs"_a = typename Mask<MaskPixelT>::MaskPlaneDict());
    cls.def(py::init<unsigned int, unsigned int, MaskPixelT,
                     typename Mask<MaskPixelT>::MaskPlaneDict const &>(),
            "width"_a, "height"_a, "initialValue"_a,
            "planeDefs"_a = typename Mask<MaskPixelT>::MaskPlaneDict());
    cls.def(py::init<lsst::geom::Extent2I const &, typename Mask<MaskPixelT>::MaskPlaneDict const &>(),
            "dimensions"_a = lsst::geom::Extent2I(),
            "planeDefs"_a = typename Mask<MaskPixelT>::MaskPlaneDict());
    cls.def(py::init<lsst::geom::Extent2I const &, MaskPixelT,
                     typename Mask<MaskPixelT>::MaskPlaneDict const &>(),
            "dimensions"_a = lsst::geom::Extent2I(), "initialValue"_a,
            "planeDefs"_a = typename Mask<MaskPixelT>::MaskPlaneDict());
    cls.def(py::init<lsst::geom::Box2I const &, typename Mask<MaskPixelT>::MaskPlaneDict const &>(), "bbox"_a,
            "planeDefs"_a = typename Mask<MaskPixelT>::MaskPlaneDict());
    cls.def(py::init<lsst::geom::Box2I const &, MaskPixelT,
                     typename Mask<MaskPixelT>::MaskPlaneDict const &>(),
            "bbox"_a, "initialValue"_a, "planeDefs"_a = typename Mask<MaskPixelT>::MaskPlaneDict());
    cls.def(py::init<const Mask<MaskPixelT> &, const bool>(), "src"_a, "deep"_a = false);
    cls.def(py::init<const Mask<MaskPixelT> &, const lsst::geom::Box2I &, ImageOrigin const, const bool>(),
            "src"_a, "bbox"_a, "origin"_a = PARENT, "deep"_a = false);
    cls.def(py::init<ndarray::Array<MaskPixelT, 2, 1> const &, bool, lsst::geom::Point2I const &>(),
            "array"_a, "deep"_a = false, "xy0"_a = lsst::geom::Point2I());
    cls.def(py::init<std::string const &, int, std::shared_ptr<lsst::daf::base::PropertySet>,
                     lsst::geom::Box2I const &, ImageOrigin, bool, bool>(),
            "fileName"_a, "hdu"_a = fits::DEFAULT_HDU, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(),
            "origin"_a = PARENT, "conformMasks"_a = false, "allowUnsafe"_a = false);
    cls.def(py::init<fits::MemFileManager &, int, std::shared_ptr<lsst::daf::base::PropertySet>,
                     lsst::geom::Box2I const &, ImageOrigin, bool, bool>(),
            "manager"_a, "hdu"_a = fits::DEFAULT_HDU, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(),
            "origin"_a = PARENT, "conformMasks"_a = false, "allowUnsafe"_a = false);
    cls.def(py::init<fits::Fits &, std::shared_ptr<lsst::daf::base::PropertySet>, lsst::geom::Box2I const &,
                     ImageOrigin, bool, bool>(),
            "fitsFile"_a, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT,
            "conformMasks"_a = false, "allowUnsafe"_a = false);

    /* Operators */
    cls.def("__ior__", [](Mask<MaskPixelT> &self, Mask<MaskPixelT> &other) { return self |= other; });
    cls.def("__ior__", [](Mask<MaskPixelT> &self, MaskPixelT const other) { return self |= other; });
    cls.def("__ior__", [](Mask<MaskPixelT> &self, int other) { return self |= other; });
    cls.def("__iand__", [](Mask<MaskPixelT> &self, Mask<MaskPixelT> &other) { return self &= other; });
    cls.def("__iand__", [](Mask<MaskPixelT> &self, MaskPixelT const other) { return self &= other; });
    cls.def("__iand__", [](Mask<MaskPixelT> &self, int other) { return self &= other; });
    cls.def("__ixor__", [](Mask<MaskPixelT> &self, Mask<MaskPixelT> &other) { return self ^= other; });
    cls.def("__ixor__", [](Mask<MaskPixelT> &self, MaskPixelT const other) { return self ^= other; });
    cls.def("__ixor__", [](Mask<MaskPixelT> &self, int other) { return self ^= other; });

    /* Members */
    cls.def("swap", (void (Mask<MaskPixelT>::*)(Mask<MaskPixelT> &)) & Mask<MaskPixelT>::swap);
    cls.def("writeFits",
            (void (Mask<MaskPixelT>::*)(std::string const &,
                                        std::shared_ptr<lsst::daf::base::PropertySet const>,
                                        std::string const &) const) &
                    Mask<MaskPixelT>::writeFits,
            "fileName"_a, "metadata"_a = std::shared_ptr<lsst::daf::base::PropertySet>(), "mode"_a = "w");
    cls.def("writeFits",
            (void (Mask<MaskPixelT>::*)(fits::MemFileManager &,
                                        std::shared_ptr<lsst::daf::base::PropertySet const>,
                                        std::string const &) const) &
                    Mask<MaskPixelT>::writeFits,
            "manager"_a, "metadata"_a = std::shared_ptr<lsst::daf::base::PropertySet>(), "mode"_a = "w");
    cls.def("writeFits",
            (void (Mask<MaskPixelT>::*)(fits::Fits &, std::shared_ptr<lsst::daf::base::PropertySet const>)
                     const) &
                    Mask<MaskPixelT>::writeFits,
            "fitsfile"_a, "metadata"_a = std::shared_ptr<lsst::daf::base::PropertySet const>());
    cls.def("writeFits",
            (void (Mask<MaskPixelT>::*)(std::string const &, fits::ImageWriteOptions const &,
                                        std::string const &, std::shared_ptr<daf::base::PropertySet const>)
                     const) &
                    Mask<MaskPixelT>::writeFits,
            "filename"_a, "options"_a, "mode"_a = "w",
            "header"_a = std::shared_ptr<daf::base::PropertyList>());
    cls.def("writeFits",
            (void (Mask<MaskPixelT>::*)(fits::MemFileManager &, fits::ImageWriteOptions const &,
                                        std::string const &, std::shared_ptr<daf::base::PropertySet const>)
                     const) &
                    Mask<MaskPixelT>::writeFits,
            "manager"_a, "options"_a, "mode"_a = "w",
            "header"_a = std::shared_ptr<daf::base::PropertyList>());
    cls.def("writeFits",
            (void (Mask<MaskPixelT>::*)(fits::Fits &, fits::ImageWriteOptions const &,
                                        std::shared_ptr<daf::base::PropertySet const>) const) &
                    Mask<MaskPixelT>::writeFits,
            "fits"_a, "options"_a, "header"_a = std::shared_ptr<daf::base::PropertyList>());
    cls.def_static("readFits", (Mask<MaskPixelT>(*)(std::string const &, int))Mask<MaskPixelT>::readFits,
                   "filename"_a, "hdu"_a = fits::DEFAULT_HDU);
    cls.def_static("readFits", (Mask<MaskPixelT>(*)(fits::MemFileManager &, int))Mask<MaskPixelT>::readFits,
                   "manager"_a, "hdu"_a = fits::DEFAULT_HDU);
    cls.def_static("interpret", Mask<MaskPixelT>::interpret);
    cls.def("subset", &Mask<MaskPixelT>::subset, "bbox"_a, "origin"_a = PARENT);
    cls.def("getAsString", &Mask<MaskPixelT>::getAsString);
    cls.def("clearAllMaskPlanes", &Mask<MaskPixelT>::clearAllMaskPlanes);
    cls.def("clearMaskPlane", &Mask<MaskPixelT>::clearMaskPlane);
    cls.def("setMaskPlaneValues", &Mask<MaskPixelT>::setMaskPlaneValues);
    cls.def_static("parseMaskPlaneMetadata", Mask<MaskPixelT>::parseMaskPlaneMetadata);
    cls.def_static("clearMaskPlaneDict", Mask<MaskPixelT>::clearMaskPlaneDict);
    cls.def_static("removeMaskPlane", Mask<MaskPixelT>::removeMaskPlane);
    cls.def("removeAndClearMaskPlane", &Mask<MaskPixelT>::removeAndClearMaskPlane, "name"_a,
            "removeFromDefault"_a = false);
    cls.def_static("getMaskPlane", Mask<MaskPixelT>::getMaskPlane);
    cls.def_static("getPlaneBitMask", (MaskPixelT(*)(const std::string &))Mask<MaskPixelT>::getPlaneBitMask);
    cls.def_static("getPlaneBitMask",
                   (MaskPixelT(*)(const std::vector<std::string> &))Mask<MaskPixelT>::getPlaneBitMask);
    cls.def_static("getNumPlanesMax", Mask<MaskPixelT>::getNumPlanesMax);
    cls.def_static("getNumPlanesUsed", Mask<MaskPixelT>::getNumPlanesUsed);
    cls.def("getMaskPlaneDict", &Mask<MaskPixelT>::getMaskPlaneDict);
    cls.def("printMaskPlanes", &Mask<MaskPixelT>::printMaskPlanes);
    cls.def_static("addMaskPlanesToMetadata", Mask<MaskPixelT>::addMaskPlanesToMetadata);
    cls.def("conformMaskPlanes", &Mask<MaskPixelT>::conformMaskPlanes);
    cls.def_static("addMaskPlane", (int (*)(const std::string &))Mask<MaskPixelT>::addMaskPlane);
}

template <typename PixelT>
static PyImage<PixelT> declareImage(py::module &mod, const std::string &suffix) {
    PyImage<PixelT> cls(mod, ("Image" + suffix).c_str());

    /* Constructors */
    cls.def(py::init<unsigned int, unsigned int, PixelT>(), "width"_a, "height"_a, "intialValue"_a = 0);
    cls.def(py::init<lsst::geom::Extent2I const &, PixelT>(), "dimensions"_a = lsst::geom::Extent2I(),
            "initialValue"_a = 0);
    cls.def(py::init<lsst::geom::Box2I const &, PixelT>(), "bbox"_a, "initialValue"_a = 0);
    cls.def(py::init<Image<PixelT> const &, lsst::geom::Box2I const &, ImageOrigin const, const bool>(),
            "rhs"_a, "bbox"_a, "origin"_a = PARENT, "deep"_a = false);
    cls.def(py::init<ndarray::Array<PixelT, 2, 1> const &, bool, lsst::geom::Point2I const &>(), "array"_a,
            "deep"_a = false, "xy0"_a = lsst::geom::Point2I());
    cls.def(py::init<std::string const &, int, std::shared_ptr<daf::base::PropertySet>,
                     lsst::geom::Box2I const &, ImageOrigin, bool>(),
            "fileName"_a, "hdu"_a = fits::DEFAULT_HDU, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(),
            "origin"_a = PARENT, "allowUnsafe"_a = false);
    cls.def(py::init<fits::MemFileManager &, int, std::shared_ptr<daf::base::PropertySet>,
                     lsst::geom::Box2I const &, ImageOrigin, bool>(),
            "manager"_a, "hdu"_a = fits::DEFAULT_HDU, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(),
            "origin"_a = PARENT, "allowUnsafe"_a = false);
    cls.def(py::init<fits::Fits &, std::shared_ptr<daf::base::PropertySet>, lsst::geom::Box2I const &,
                     ImageOrigin, bool>(),
            "fitsFile"_a, "metadata"_a = nullptr, "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT,
            "allowUnsafe"_a = false);

    /* Operators */
    cls.def("__iadd__", [](Image<PixelT> &self, PixelT const &other) { return self += other; });
    cls.def("__iadd__", [](Image<PixelT> &self, Image<PixelT> const &other) { return self += other; });
    cls.def("__iadd__", [](Image<PixelT> &self, lsst::afw::math::Function2<double> const &other) {
        return self += other;
    });
    cls.def("__isub__", [](Image<PixelT> &self, PixelT const &other) { return self -= other; });
    cls.def("__isub__", [](Image<PixelT> &self, Image<PixelT> const &other) { return self -= other; });
    cls.def("__isub__", [](Image<PixelT> &self, lsst::afw::math::Function2<double> const &other) {
        return self -= other;
    });
    cls.def("__imul__", [](Image<PixelT> &self, PixelT const &other) { return self *= other; });
    cls.def("__imul__", [](Image<PixelT> &self, Image<PixelT> const &other) { return self *= other; });
    cls.def("__itruediv__", [](Image<PixelT> &self, PixelT const &other) { return self /= other; });
    cls.def("__itruediv__", [](Image<PixelT> &self, Image<PixelT> const &other) { return self /= other; });

    /* Members */
    cls.def("scaledPlus", &Image<PixelT>::scaledPlus);
    cls.def("scaledMinus", &Image<PixelT>::scaledMinus);
    cls.def("scaledMultiplies", &Image<PixelT>::scaledMultiplies);
    cls.def("scaledDivides", &Image<PixelT>::scaledDivides);

    cls.def("subset", &Image<PixelT>::subset, "bbox"_a, "origin"_a = PARENT);

    cls.def("writeFits",
            (void (Image<PixelT>::*)(std::string const &, std::shared_ptr<daf::base::PropertySet const>,
                                     std::string const &) const) &
                    Image<PixelT>::writeFits,
            "fileName"_a, "metadata"_a = std::shared_ptr<daf::base::PropertySet const>(), "mode"_a = "w");
    cls.def("writeFits",
            (void (Image<PixelT>::*)(fits::MemFileManager &, std::shared_ptr<daf::base::PropertySet const>,
                                     std::string const &) const) &
                    Image<PixelT>::writeFits,
            "manager"_a, "metadata"_a = std::shared_ptr<daf::base::PropertySet const>(), "mode"_a = "w");
    cls.def("writeFits",
            (void (Image<PixelT>::*)(fits::Fits &, std::shared_ptr<daf::base::PropertySet const>) const) &
                    Image<PixelT>::writeFits,
            "fitsfile"_a, "metadata"_a = std::shared_ptr<daf::base::PropertySet const>());
    cls.def("writeFits",
            (void (Image<PixelT>::*)(std::string const &, fits::ImageWriteOptions const &,
                                     std::string const &, std::shared_ptr<daf::base::PropertySet const>,
                                     std::shared_ptr<image::Mask<image::MaskPixel> const>) const) &
                    Image<PixelT>::writeFits,
            "filename"_a, "options"_a, "mode"_a = "w",
            "header"_a = std::shared_ptr<daf::base::PropertyList>(),
            "mask"_a = std::shared_ptr<image::Mask<image::MaskPixel>>());
    cls.def("writeFits",
            (void (Image<PixelT>::*)(fits::MemFileManager &, fits::ImageWriteOptions const &,
                                     std::string const &, std::shared_ptr<daf::base::PropertySet const>,
                                     std::shared_ptr<image::Mask<image::MaskPixel> const>) const) &
                    Image<PixelT>::writeFits,
            "manager"_a, "options"_a, "mode"_a = "w", "header"_a = std::shared_ptr<daf::base::PropertyList>(),
            "mask"_a = std::shared_ptr<image::Mask<image::MaskPixel>>());
    cls.def("writeFits",
            (void (Image<PixelT>::*)(fits::Fits &, fits::ImageWriteOptions const &,
                                     std::shared_ptr<daf::base::PropertySet const>,
                                     std::shared_ptr<image::Mask<image::MaskPixel> const>) const) &
                    Image<PixelT>::writeFits,
            "fits"_a, "options"_a, "header"_a = std::shared_ptr<daf::base::PropertyList>(),
            "mask"_a = std::shared_ptr<image::Mask<image::MaskPixel>>());

    cls.def_static("readFits", (Image<PixelT>(*)(std::string const &, int))Image<PixelT>::readFits,
                   "filename"_a, "hdu"_a = fits::DEFAULT_HDU);
    cls.def_static("readFits", (Image<PixelT>(*)(fits::MemFileManager &, int))Image<PixelT>::readFits,
                   "manager"_a, "hdu"_a = fits::DEFAULT_HDU);
    cls.def("sqrt", &Image<PixelT>::sqrt);

    return cls;
}

template <typename PixelT>
static void declareDecoratedImage(py::module &mod, std::string const &suffix) {
    PyDecoratedImage<PixelT> cls(mod, ("DecoratedImage" + suffix).c_str());

    cls.def(py::init<const lsst::geom::Extent2I &>(), "dimensions"_a = lsst::geom::Extent2I());
    cls.def(py::init<const lsst::geom::Box2I &>(), "bbox"_a);
    cls.def(py::init<std::shared_ptr<Image<PixelT>>>(), "rhs"_a);
    cls.def(py::init<DecoratedImage<PixelT> const &, const bool>(), "rhs"_a, "deep"_a = false);
    cls.def(py::init<std::string const &, const int, lsst::geom::Box2I const &, ImageOrigin const, bool>(),
            "fileName"_a, "hdu"_a = fits::DEFAULT_HDU, "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT,
            "allowUnsafe"_a = false);

    cls.def("getMetadata", &DecoratedImage<PixelT>::getMetadata);
    cls.def("setMetadata", &DecoratedImage<PixelT>::setMetadata);
    cls.def("getWidth", &DecoratedImage<PixelT>::getWidth);
    cls.def("getHeight", &DecoratedImage<PixelT>::getHeight);
    cls.def("getX0", &DecoratedImage<PixelT>::getX0);
    cls.def("getY0", &DecoratedImage<PixelT>::getY0);
    cls.def("getDimensions", &DecoratedImage<PixelT>::getDimensions);
    cls.def("swap", &DecoratedImage<PixelT>::swap);
    cls.def("writeFits",
            py::overload_cast<std::string const &, std::shared_ptr<daf::base::PropertySet const>,
                              std::string const &>(&DecoratedImage<PixelT>::writeFits, py::const_),
            "filename"_a, "metadata"_a = std::shared_ptr<daf::base::PropertyList>(), "mode"_a = "w");
    cls.def("writeFits",
            py::overload_cast<std::string const &, fits::ImageWriteOptions const &,
                              std::shared_ptr<daf::base::PropertySet const>, std::string const &>(
                    &DecoratedImage<PixelT>::writeFits, py::const_),
            "filename"_a, "options"_a, "metadata"_a = std::shared_ptr<daf::base::PropertyList>(),
            "mode"_a = "w");
    cls.def("getImage", py::overload_cast<>(&DecoratedImage<PixelT>::getImage));
    cls.def_property_readonly("image", py::overload_cast<>(&DecoratedImage<PixelT>::getImage));
    cls.def("getGain", &DecoratedImage<PixelT>::getGain);
    cls.def("setGain", &DecoratedImage<PixelT>::setGain);
}

/* Declare ImageSlice operators separately since they are only instantiated for float double */
template <typename PixelT>
static void addImageSliceOperators(
        py::class_<Image<PixelT>, std::shared_ptr<Image<PixelT>>, ImageBase<PixelT>> &cls) {
    cls.def("__add__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self + other; },
            py::is_operator());
    cls.def("__sub__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self - other; },
            py::is_operator());
    cls.def("__mul__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self * other; },
            py::is_operator());
    cls.def("__truediv__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self / other; },
            py::is_operator());
    cls.def("__iadd__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self += other;
        return self;
    });
    cls.def("__isub__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self -= other;
        return self;
    });
    cls.def("__imul__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self *= other;
        return self;
    });
    cls.def("__itruediv__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self /= other;
        return self;
    });
}

template <typename PixelT, typename PyClass>
static void addGeneralizedCopyConstructors(PyClass &cls) {
    cls.def(py::init<Image<int> const &, const bool>(), "rhs"_a, "deep"_a = false);
    cls.def(py::init<Image<float> const &, const bool>(), "rhs"_a, "deep"_a = false);
    cls.def(py::init<Image<double> const &, const bool>(), "rhs"_a, "deep"_a = false);
    cls.def(py::init<Image<std::uint16_t> const &, const bool>(), "rhs"_a, "deep"_a = false);
    cls.def(py::init<Image<std::uint64_t> const &, const bool>(), "rhs"_a, "deep"_a = false);

    cls.def("convertI", [](Image<PixelT> const &self) { return Image<int>(self, true); });
    cls.def("convertF", [](Image<PixelT> const &self) { return Image<float>(self, true); });
    cls.def("convertD", [](Image<PixelT> const &self) { return Image<double>(self, true); });
    cls.def("convertU", [](Image<PixelT> const &self) { return Image<std::uint16_t>(self, true); });
    cls.def("convertL", [](Image<PixelT> const &self) { return Image<std::uint64_t>(self, true); });

    cls.def("convertFloat", [](Image<PixelT> const &self) { return Image<float>(self, true); });
    cls.def("convertDouble", [](Image<PixelT> const &self) { return Image<double>(self, true); });
}

WRAP(Img) {

    py::enum_<ImageOrigin>(mod, "ImageOrigin")
            .value("PARENT", ImageOrigin::PARENT)
            .value("LOCAL", ImageOrigin::LOCAL)
            .export_values();

    declareImageBase<int>(mod, "I");
    declareImageBase<float>(mod, "F");
    declareImageBase<double>(mod, "D");
    declareImageBase<std::uint16_t>(mod, "U");
    declareImageBase<std::uint64_t>(mod, "L");

    // Mask must be declared before Image because a mask is used as a default value in at least one method
    declareMask<MaskPixel>(mod, "X");

    auto clsImageI = declareImage<int>(mod, "I");
    auto clsImageF = declareImage<float>(mod, "F");
    auto clsImageD = declareImage<double>(mod, "D");
    auto clsImageU = declareImage<std::uint16_t>(mod, "U");
    auto clsImageL = declareImage<std::uint64_t>(mod, "L");

    // Add generalized copy constructors
    addGeneralizedCopyConstructors<int>(clsImageI);
    addGeneralizedCopyConstructors<float>(clsImageF);
    addGeneralizedCopyConstructors<double>(clsImageD);
    addGeneralizedCopyConstructors<std::uint16_t>(clsImageU);
    addGeneralizedCopyConstructors<std::uint64_t>(clsImageL);

    // Add slice operators only for float and double
    addImageSliceOperators<float>(clsImageF);
    addImageSliceOperators<double>(clsImageD);

    declareDecoratedImage<int>(mod, "I");
    declareDecoratedImage<float>(mod, "F");
    declareDecoratedImage<double>(mod, "D");
    declareDecoratedImage<std::uint16_t>(mod, "U");
    declareDecoratedImage<std::uint64_t>(mod, "L");

    // Declare constructors for casting all exposure types to to float and double
    // (the only two types of casts that Python supports)
    declareCastConstructor<int, float>(clsImageF);
    declareCastConstructor<int, double>(clsImageD);

    declareCastConstructor<float, double>(clsImageD);

    declareCastConstructor<double, float>(clsImageF);

    declareCastConstructor<std::uint16_t, float>(clsImageF);
    declareCastConstructor<std::uint16_t, double>(clsImageD);

    declareCastConstructor<std::uint64_t, float>(clsImageF);
    declareCastConstructor<std::uint64_t, double>(clsImageD);

    // Note: wrap both the Image and MaskedImage versions of imagesOverlap in the MaskedImage wrapper,
    // as wrapping the Image version here results in it being invisible in lsst.afw.image

    mod.def("bboxFromMetadata", &bboxFromMetadata);
}

template <typename PixelT>
using PyExposure = py::class_<Exposure<PixelT>, std::shared_ptr<Exposure<PixelT>>>;

/*
Declare a constructor that takes an Exposure of FromPixelT and returns an Exposure cast to ToPixelT

The mask and variance must be of the standard types.

@param[in] cls  The pybind11 class to which add the constructor
*/
template <typename FromPixelT, typename ToPixelT>
void declareCastConstructor(PyExposure<ToPixelT> &cls) {
    cls.def(py::init<Exposure<FromPixelT> const &, bool const>(), "src"_a, "deep"_a);
}

template <typename PixelT>  // only the image type varies; mask and variance are fixed to default tparams
PyExposure<PixelT> declareExposure(py::module &mod, const std::string &suffix) {
    using ExposureT = Exposure<PixelT>;
    using MaskedImageT = typename ExposureT::MaskedImageT;

    PyExposure<PixelT> cls(mod, ("Exposure" + suffix).c_str());

    mod.def("makeExposure", &makeExposure<PixelT, MaskPixel, VariancePixel>, "maskedImage"_a,
            "wcs"_a = std::shared_ptr<geom::SkyWcs const>());

    /* Constructors */
    cls.def(py::init<unsigned int, unsigned int, std::shared_ptr<geom::SkyWcs const>>(), "width"_a,
            "height"_a, "wcs"_a = std::shared_ptr<geom::SkyWcs const>());
    cls.def(py::init<lsst::geom::Extent2I const &, std::shared_ptr<geom::SkyWcs const>>(),
            "dimensions"_a = lsst::geom::Extent2I(), "wcs"_a = std::shared_ptr<geom::SkyWcs const>());
    cls.def(py::init<lsst::geom::Box2I const &, std::shared_ptr<geom::SkyWcs const>>(), "bbox"_a,
            "wcs"_a = std::shared_ptr<geom::SkyWcs const>());
    cls.def(py::init<MaskedImageT &, std::shared_ptr<geom::SkyWcs const>>(), "maskedImage"_a,
            "wcs"_a = std::shared_ptr<geom::SkyWcs const>());
    cls.def(py::init<MaskedImageT &, std::shared_ptr<ExposureInfo>>(), "maskedImage"_a, "exposureInfo"_a);
    cls.def(py::init<std::string const &, lsst::geom::Box2I const &, ImageOrigin, bool, bool>(), "fileName"_a,
            "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT, "conformMasks"_a = false,
            "allowUnsafe"_a = false);
    cls.def(py::init<fits::MemFileManager &, lsst::geom::Box2I const &, ImageOrigin, bool, bool>(),
            "manager"_a, "bbox"_a = lsst::geom::Box2I(), "origin"_a = PARENT, "conformMasks"_a = false,
            "allowUnsafe"_a = false);
    cls.def(py::init<ExposureT const &, bool>(), "other"_a, "deep"_a = false);
    cls.def(py::init<ExposureT const &, lsst::geom::Box2I const &, ImageOrigin, bool>(), "other"_a, "bbox"_a,
            "origin"_a = PARENT, "deep"_a = false);

    /* Members */
    cls.def("getMaskedImage", (MaskedImageT(ExposureT::*)()) & ExposureT::getMaskedImage);
    cls.def("setMaskedImage", &ExposureT::setMaskedImage, "maskedImage"_a);
    cls.def_property("maskedImage", (MaskedImageT(ExposureT::*)()) & ExposureT::getMaskedImage,
                     &ExposureT::setMaskedImage);
    cls.def("getMetadata", &ExposureT::getMetadata);
    cls.def("setMetadata", &ExposureT::setMetadata, "metadata"_a);
    cls.def("getWidth", &ExposureT::getWidth);
    cls.def("getHeight", &ExposureT::getHeight);
    cls.def("getDimensions", &ExposureT::getDimensions);
    cls.def("getX0", &ExposureT::getX0);
    cls.def("getY0", &ExposureT::getY0);
    cls.def("getXY0", &ExposureT::getXY0);
    cls.def("setXY0", &ExposureT::setXY0, "xy0"_a);
    cls.def("getBBox", &ExposureT::getBBox, "origin"_a = PARENT);
    cls.def("getWcs", (std::shared_ptr<geom::SkyWcs>(ExposureT::*)()) & ExposureT::getWcs);
    cls.def("setWcs", &ExposureT::setWcs, "wcs"_a);
    cls.def("hasWcs", &ExposureT::hasWcs);
    cls.def("getDetector", &ExposureT::getDetector);
    cls.def("setDetector", &ExposureT::setDetector, "detector"_a);
    cls.def("getFilter", &ExposureT::getFilter);
    cls.def("setFilter", &ExposureT::setFilter, "filter"_a);

    cls.def("getPhotoCalib", &ExposureT::getPhotoCalib);
    cls.def("setPhotoCalib", &ExposureT::setPhotoCalib, "photoCalib"_a);
    cls.def("getPsf", (std::shared_ptr<detection::Psf>(ExposureT::*)()) & ExposureT::getPsf);
    cls.def("setPsf", &ExposureT::setPsf, "psf"_a);
    cls.def("hasPsf", &ExposureT::hasPsf);
    cls.def("getInfo", (std::shared_ptr<ExposureInfo>(ExposureT::*)()) & ExposureT::getInfo);
    cls.def("setInfo", &ExposureT::setInfo, "exposureInfo"_a);

    cls.def("subset", &ExposureT::subset, "bbox"_a, "origin"_a = PARENT);

    cls.def("writeFits", (void (ExposureT::*)(std::string const &) const) & ExposureT::writeFits);
    cls.def("writeFits", (void (ExposureT::*)(fits::MemFileManager &) const) & ExposureT::writeFits);
    cls.def("writeFits", [](ExposureT &self, fits::Fits &fits) { self.writeFits(fits); });

    cls.def("writeFits",
            [](ExposureT &self, std::string const &filename, fits::ImageWriteOptions const &imageOptions,
               fits::ImageWriteOptions const &maskOptions, fits::ImageWriteOptions const &varianceOptions) {
                self.writeFits(filename, imageOptions, maskOptions, varianceOptions);
            },
            "filename"_a, "imageOptions"_a, "maskOptions"_a, "varianceOptions"_a);
    cls.def("writeFits",
            [](ExposureT &self, fits::MemFileManager &manager, fits::ImageWriteOptions const &imageOptions,
               fits::ImageWriteOptions const &maskOptions, fits::ImageWriteOptions const &varianceOptions) {
                self.writeFits(manager, imageOptions, maskOptions, varianceOptions);
            },
            "manager"_a, "imageOptions"_a, "maskOptions"_a, "varianceOptions"_a);
    cls.def("writeFits",
            [](ExposureT &self, fits::Fits &fits, fits::ImageWriteOptions const &imageOptions,
               fits::ImageWriteOptions const &maskOptions, fits::ImageWriteOptions const &varianceOptions) {
                self.writeFits(fits, imageOptions, maskOptions, varianceOptions);
            },
            "fits"_a, "imageOptions"_a, "maskOptions"_a, "varianceOptions"_a);

    cls.def_static("readFits", (ExposureT(*)(std::string const &))ExposureT::readFits);
    cls.def_static("readFits", (ExposureT(*)(fits::MemFileManager &))ExposureT::readFits);

    cls.def("getCutout", &ExposureT::getCutout, "center"_a, "size"_a);

    return cls;
}

WRAP(Exposure) {

    auto clsExposureF = declareExposure<float>(mod, "F");
    auto clsExposureD = declareExposure<double>(mod, "D");
    declareExposure<int>(mod, "I");
    declareExposure<std::uint16_t>(mod, "U");
    declareExposure<std::uint64_t>(mod, "L");

    // Declare constructors for casting all exposure types to to float and double
    // (the only two types of casts that Python supports)
    declareCastConstructor<int, float>(clsExposureF);
    declareCastConstructor<int, double>(clsExposureD);

    declareCastConstructor<float, double>(clsExposureD);

    declareCastConstructor<double, float>(clsExposureF);

    declareCastConstructor<std::uint16_t, float>(clsExposureF);
    declareCastConstructor<std::uint16_t, double>(clsExposureD);

    declareCastConstructor<std::uint64_t, float>(clsExposureF);
    declareCastConstructor<std::uint64_t, double>(clsExposureD);
}

using PyApCorrMap = py::class_<ApCorrMap, std::shared_ptr<ApCorrMap>, typehandling::Storable>;

WRAP(ApCorrMap) {

    /* Declare CRTP base class. */
    PyApCorrMap cls(mod, "ApCorrMap");

    /* Constructors */
    cls.def(py::init<>());

    table::io::python::addPersistableMethods<ApCorrMap>(cls);

    /* Operators */
    cls.def("__imul__", &ApCorrMap::operator*=);
    cls.def("__itruediv__", &ApCorrMap::operator/=);

    /* Members */
    cls.def("get", &ApCorrMap::get);
    cls.def("set", &ApCorrMap::set);
    cls.def("items", [](ApCorrMap const& self) { return py::make_iterator(self.begin(), self.end()); },
            py::keep_alive<0, 1>());
    // values, keys, and __iter__ defined in apCorrMap.py

    cls.def("__len__", &ApCorrMap::size);
    cls.def("__getitem__", &ApCorrMap::operator[]);
    cls.def("__setitem__", &ApCorrMap::set);
    cls.def("__contains__",
            [](ApCorrMap const& self, std::string name) { return static_cast<bool>(self.get(name)); });
}

static double const nan(std::numeric_limits<double>::quiet_NaN());
static lsst::geom::Angle const nanAngle(nan);

using PyTransmissionCurve =
        py::class_<TransmissionCurve, std::shared_ptr<TransmissionCurve>, typehandling::Storable>;

PyTransmissionCurve declare(py::module & mod) {
    return PyTransmissionCurve(mod, "TransmissionCurve");
}

void define(PyTransmissionCurve & cls) {
    table::io::python::addPersistableMethods(cls);

    cls.def_static("makeIdentity", &TransmissionCurve::makeIdentity);
    cls.def_static("makeSpatiallyConstant", &TransmissionCurve::makeSpatiallyConstant,
                   "throughput"_a, "wavelengths"_a,
                   "throughputAtMin"_a=0.0, "throughputAtMax"_a=0.0);
    cls.def_static("makeRadial", &TransmissionCurve::makeRadial,
                   "throughput"_a, "wavelengths"_a, "radii"_a,
                   "throughputAtMin"_a=0.0, "throughputAtMax"_a=0.0);
    cls.def("__mul__", &TransmissionCurve::multipliedBy, py::is_operator());
    cls.def("multipliedBy", &TransmissionCurve::multipliedBy);
    cls.def("transformedBy", &TransmissionCurve::transformedBy, "transform"_a);
    cls.def("getWavelengthBounds", &TransmissionCurve::getWavelengthBounds);
    cls.def("getThroughputAtBounds", &TransmissionCurve::getThroughputAtBounds);
    cls.def(
        "sampleAt",
        (void (TransmissionCurve::*)(
            lsst::geom::Point2D const &,
            ndarray::Array<double const,1,1> const &,
            ndarray::Array<double,1,1> const &
        ) const) &TransmissionCurve::sampleAt,
        "position"_a, "wavelengths"_a, "out"_a
    );
    cls.def(
        "sampleAt",
        (ndarray::Array<double,1,1> (TransmissionCurve::*)(
            lsst::geom::Point2D const &,
            ndarray::Array<double const,1,1> const &
        ) const) &TransmissionCurve::sampleAt,
        "position"_a, "wavelengths"_a
    );
}

WRAP(TransmissionCurve) {
    // then declare classes
    auto cls = declare(mod);
    // and now we can safely define methods and other attributes
    define(cls);
}

// ImageBaseFitsReader is an implementation detail and is not exposed directly
// to Python, as we have better ways to share wrapper code between classes
// at the pybind11 level (e.g. declareCommon below).
using PyImageFitsReader = py::class_<ImageFitsReader, std::shared_ptr<ImageFitsReader>>;
using PyMaskFitsReader = py::class_<MaskFitsReader, std::shared_ptr<MaskFitsReader>>;
using PyMaskedImageFitsReader = py::class_<MaskedImageFitsReader, std::shared_ptr<MaskedImageFitsReader>>;
using PyExposureFitsReader = py::class_<ExposureFitsReader, std::shared_ptr<ExposureFitsReader>>;

// Declare attributes common to all FitsReaders.  Excludes constructors
// because ExposureFitsReader's don't take an HDU argument.
template <typename Class, typename ...Args>
void declareCommonMethods(py::class_<Class, Args...> & cls) {
    cls.def("readBBox", &Class::readBBox, "origin"_a=PARENT);
    cls.def("readXY0", &Class::readXY0, "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT);
    cls.def("getFileName", &Class::getFileName);
    cls.def_property_readonly("fileName", &Class::getFileName);
}

// Declare attributes common to ImageFitsReader and MaskFitsReader
template <typename Class, typename ...Args>
void declareSinglePlaneMethods(py::class_<Class, Args...> & cls) {
    cls.def(py::init<std::string const &, int>(), "fileName"_a, "hdu"_a=fits::DEFAULT_HDU);
    cls.def(py::init<fits::MemFileManager&, int>(), "manager"_a, "hdu"_a=fits::DEFAULT_HDU);
    cls.def("readMetadata", &Class::readMetadata);
    cls.def("readDType", [](Class & self) { return py::dtype(self.readDType()); });
    cls.def("getHdu", &Class::getHdu);
    cls.def_property_readonly("hdu", &Class::getHdu);
    cls.def(
        "readArray",
        [](Class & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool allowUnsafe,
           py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.template readArray<decltype(t)>(bbox, origin, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<std::uint16_t, int, float, double, std::uint64_t>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
}

// Declare attributes shared by MaskedImageFitsReader and MaskedImageFitsReader.
template <typename Class, typename ...Args>
void declareMultiPlaneMethods(py::class_<Class, Args...> & cls) {
    cls.def("readImageDType", [](Class & self) { return py::dtype(self.readImageDType()); } );
    cls.def("readMaskDType", [](Class & self) { return py::dtype(self.readMaskDType()); });
    cls.def("readVarianceDType", [](Class & self) { return py::dtype(self.readVarianceDType()); });
    cls.def(
        "readImage",
        [](Class & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool allowUnsafe,
           py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readImageDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.template readImage<decltype(t)>(bbox, origin, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<std::uint16_t, int, float, double, std::uint64_t>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
    cls.def(
        "readImageArray",
        [](Class & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool allowUnsafe,
           py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readImageDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.template readImageArray<decltype(t)>(bbox, origin, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<std::uint16_t, int, float, double, std::uint64_t>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
    cls.def(
        "readMask",
        [](Class & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool conformMasks,
           bool allowUnsafe, py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readMaskDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.template readMask<decltype(t)>(bbox, origin, conformMasks, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<MaskPixel>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "conformMasks"_a=false, "allowUnsafe"_a=false,
        "dtype"_a=py::none()
    );
    cls.def(
        "readMaskArray",
        [](Class & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool allowUnsafe,
           py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readMaskDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.template readMaskArray<decltype(t)>(bbox, origin, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<MaskPixel>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
    cls.def(
        "readVariance",
        [](Class & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool allowUnsafe,
           py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readVarianceDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.template readVariance<decltype(t)>(bbox, origin, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<VariancePixel>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
    cls.def(
        "readVarianceArray",
        [](Class & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool allowUnsafe,
           py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readVarianceDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.template readVarianceArray<decltype(t)>(bbox, origin, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<VariancePixel>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
}

void declareImageFitsReader(py::module & mod) {
    PyImageFitsReader cls(mod, "ImageFitsReader");
    declareCommonMethods(cls);
    declareSinglePlaneMethods(cls);
    cls.def(
        "read",
        [](ImageFitsReader & self, lsst::geom::Box2I const & bbox, ImageOrigin origin, bool allowUnsafe,
           py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.read<decltype(t)>(bbox, origin, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<std::uint16_t, int, float, double, std::uint64_t>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
}

void declareMaskFitsReader(py::module & mod) {
    PyMaskFitsReader cls(mod, "MaskFitsReader");
    declareCommonMethods(cls);
    declareSinglePlaneMethods(cls);
    cls.def(
        "read",
        [](MaskFitsReader & self, lsst::geom::Box2I const & bbox, ImageOrigin origin,
           bool conformMasks, bool allowUnsafe, py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.read<decltype(t)>(bbox, origin, conformMasks, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<MaskPixel>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "conformMasks"_a=false, "allowUnsafe"_a=false,
        "dtype"_a=py::none()
    );
    // all other methods provided by base class wrappers
}

void declareMaskedImageFitsReader(py::module & mod) {
    PyMaskedImageFitsReader cls(mod, "MaskedImageFitsReader");
    cls.def(py::init<std::string const &, int>(), "fileName"_a, "hdu"_a=fits::DEFAULT_HDU);
    cls.def(py::init<fits::MemFileManager&, int>(), "manager"_a, "hdu"_a=fits::DEFAULT_HDU);
    declareCommonMethods(cls);
    declareMultiPlaneMethods(cls);
    cls.def("readPrimaryMetadata", &MaskedImageFitsReader::readPrimaryMetadata);
    cls.def("readImageMetadata", &MaskedImageFitsReader::readImageMetadata);
    cls.def("readMaskMetadata", &MaskedImageFitsReader::readMaskMetadata);
    cls.def("readVarianceMetadata", &MaskedImageFitsReader::readVarianceMetadata);
    cls.def(
        "read",
        [](MaskedImageFitsReader & self, lsst::geom::Box2I const & bbox, ImageOrigin origin,
           bool conformMasks, bool needAllHdus, bool allowUnsafe, py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readImageDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.read<decltype(t)>(bbox, origin, conformMasks, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<std::uint16_t, int, float, double, std::uint64_t>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "conformMasks"_a=false, "needAllHdus"_a=false,
        "allowUnsafe"_a=false, "dtype"_a=py::none()
    );
}

void declareExposureFitsReader(py::module & mod) {
    PyExposureFitsReader cls(mod, "ExposureFitsReader");
    cls.def(py::init<std::string const &>(), "fileName"_a);
    cls.def(py::init<fits::MemFileManager&>(), "manager"_a);
    declareCommonMethods(cls);
    declareMultiPlaneMethods(cls);
    cls.def("readMetadata", &ExposureFitsReader::readMetadata);
    cls.def("readWcs", &ExposureFitsReader::readWcs);
    cls.def("readFilter", &ExposureFitsReader::readFilter);
    cls.def("readPhotoCalib", &ExposureFitsReader::readPhotoCalib);
    cls.def("readPsf", &ExposureFitsReader::readPsf);
    cls.def("readValidPolygon", &ExposureFitsReader::readValidPolygon);
    cls.def("readApCorrMap", &ExposureFitsReader::readApCorrMap);
    cls.def("readCoaddInputs", &ExposureFitsReader::readCoaddInputs);
    cls.def("readVisitInfo", &ExposureFitsReader::readVisitInfo);
    cls.def("readTransmissionCurve", &ExposureFitsReader::readTransmissionCurve);
    cls.def("readDetector", &ExposureFitsReader::readDetector);
    cls.def("readExposureInfo", &ExposureFitsReader::readExposureInfo);
    cls.def(
        "readMaskedImage",
        [](ExposureFitsReader & self, lsst::geom::Box2I const & bbox, ImageOrigin origin,
           bool conformMasks, bool allowUnsafe, py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readImageDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.readMaskedImage<decltype(t)>(bbox, origin, conformMasks, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<std::uint16_t, int, float, double, std::uint64_t>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "conformMasks"_a=false, "allowUnsafe"_a=false,
        "dtype"_a=py::none()
    );
    cls.def(
        "read",
        [](ExposureFitsReader & self, lsst::geom::Box2I const & bbox, ImageOrigin origin,
           bool conformMasks, bool allowUnsafe, py::object dtype) {
            if (dtype.is(py::none())) {
                dtype = py::dtype(self.readImageDType());
            }
            return utils::python::TemplateInvoker().apply(
                [&](auto t) {
                    return self.read<decltype(t)>(bbox, origin, conformMasks, allowUnsafe);
                },
                py::dtype(dtype),
                utils::python::TemplateInvoker::Tag<std::uint16_t, int, float, double, std::uint64_t>()
            );
        },
        "bbox"_a=lsst::geom::Box2I(), "origin"_a=PARENT, "conformMasks"_a=false, "allowUnsafe"_a=false,
        "dtype"_a=py::none()
    );
}


WRAP(Readers) {
    declareImageFitsReader(mod);
    declareMaskFitsReader(mod);
    declareMaskedImageFitsReader(mod);
    declareExposureFitsReader(mod);
}


void declareMeasurement(py::module &mod) {
    py::class_<Measurement, std::shared_ptr<Measurement>> cls(mod, "Measurement");

    cls.def(py::init<double, double>(), "value"_a, "error"_a);
    cls.def_readonly("value", &Measurement::value);
    cls.def_readonly("error", &Measurement::error);

    utils::python::addOutputOp(cls, "__str__");
    cls.def("__repr__", [](Measurement const &self) {
        std::ostringstream os;
        os << "Measurement(" << self << ")";
        return os.str();
    });
}

WRAP(Photocalib) {

    declareMeasurement(mod);

    py::class_<PhotoCalib, std::shared_ptr<PhotoCalib>, typehandling::Storable> cls(mod, "PhotoCalib");

    /* Constructors */
    cls.def(py::init<>());
    cls.def(py::init<double, double, lsst::geom::Box2I>(), "calibrationMean"_a, "calibrationErr"_a = 0.0,
            "bbox"_a = lsst::geom::Box2I());
    cls.def(py::init<std::shared_ptr<afw::math::BoundedField>, double>(), "calibration"_a,
            "calibrationErr"_a = 0.0);
    cls.def(py::init<double, double, std::shared_ptr<afw::math::BoundedField>, bool>(), "calibrationMean"_a,
            "calibrationErr"_a, "calibration"_a, "isConstant"_a);

    table::io::python::addPersistableMethods<PhotoCalib>(cls);

    /* Members - nanojansky */
    cls.def("instFluxToNanojansky",
            (double (PhotoCalib::*)(double, lsst::geom::Point<double, 2> const &) const) &
                    PhotoCalib::instFluxToNanojansky,
            "instFlux"_a, "point"_a);
    cls.def("instFluxToNanojansky", (double (PhotoCalib::*)(double) const) & PhotoCalib::instFluxToNanojansky,
            "instFlux"_a);

    cls.def("instFluxToNanojansky",
            (Measurement(PhotoCalib::*)(double, double, lsst::geom::Point<double, 2> const &) const) &
                    PhotoCalib::instFluxToNanojansky,
            "instFlux"_a, "instFluxErr"_a, "point"_a);
    cls.def("instFluxToNanojansky",
            (Measurement(PhotoCalib::*)(double, double) const) & PhotoCalib::instFluxToNanojansky,
            "instFlux"_a, "instFluxErr"_a);

    cls.def("instFluxToNanojansky",
            (Measurement(PhotoCalib::*)(afw::table::SourceRecord const &, std::string const &) const) &
                    PhotoCalib::instFluxToNanojansky,
            "sourceRecord"_a, "instFluxField"_a);

    cls.def("instFluxToNanojansky",
            (ndarray::Array<double, 2, 2>(PhotoCalib::*)(afw::table::SourceCatalog const &,
                                                         std::string const &) const) &
                    PhotoCalib::instFluxToNanojansky,
            "sourceCatalog"_a, "instFluxField"_a);

    cls.def("instFluxToNanojansky",
            (void (PhotoCalib::*)(afw::table::SourceCatalog &, std::string const &, std::string const &)
                     const) &
                    PhotoCalib::instFluxToNanojansky,
            "sourceCatalog"_a, "instFluxField"_a, "outField"_a);

    /* Members - magnitudes */
    cls.def("instFluxToMagnitude",
            (double (PhotoCalib::*)(double, lsst::geom::Point<double, 2> const &) const) &
                    PhotoCalib::instFluxToMagnitude,
            "instFlux"_a, "point"_a);
    cls.def("instFluxToMagnitude", (double (PhotoCalib::*)(double) const) & PhotoCalib::instFluxToMagnitude,
            "instFlux"_a);

    cls.def("instFluxToMagnitude",
            (Measurement(PhotoCalib::*)(double, double, lsst::geom::Point<double, 2> const &) const) &
                    PhotoCalib::instFluxToMagnitude,
            "instFlux"_a, "instFluxErr"_a, "point"_a);
    cls.def("instFluxToMagnitude",
            (Measurement(PhotoCalib::*)(double, double) const) & PhotoCalib::instFluxToMagnitude,
            "instFlux"_a, "instFluxErr"_a);

    cls.def("instFluxToMagnitude",
            (Measurement(PhotoCalib::*)(afw::table::SourceRecord const &, std::string const &) const) &
                    PhotoCalib::instFluxToMagnitude,
            "sourceRecord"_a, "instFluxField"_a);

    cls.def("instFluxToMagnitude",
            (ndarray::Array<double, 2, 2>(PhotoCalib::*)(afw::table::SourceCatalog const &,
                                                         std::string const &) const) &
                    PhotoCalib::instFluxToMagnitude,
            "sourceCatalog"_a, "instFluxField"_a);

    cls.def("instFluxToMagnitude",
            (void (PhotoCalib::*)(afw::table::SourceCatalog &, std::string const &, std::string const &)
                     const) &
                    PhotoCalib::instFluxToMagnitude,
            "sourceCatalog"_a, "instFluxField"_a, "outField"_a);

    /* from magnitude. */
    cls.def("magnitudeToInstFlux",
            py::overload_cast<double, lsst::geom::Point<double, 2> const &>(&PhotoCalib::magnitudeToInstFlux,
                                                                            py::const_),
            "instFlux"_a, "point"_a);
    cls.def("magnitudeToInstFlux", py::overload_cast<double>(&PhotoCalib::magnitudeToInstFlux, py::const_),
            "instFlux"_a);

    /* utilities */
    cls.def("getCalibrationMean", &PhotoCalib::getCalibrationMean);
    cls.def("getCalibrationErr", &PhotoCalib::getCalibrationErr);
    cls.def("getInstFluxAtZeroMagnitude", &PhotoCalib::getInstFluxAtZeroMagnitude);
    cls.def("getLocalCalibration", &PhotoCalib::getLocalCalibration, "point"_a);

    cls.def("computeScaledCalibration", &PhotoCalib::computeScaledCalibration);
    cls.def("computeScalingTo", &PhotoCalib::computeScalingTo);

    cls.def("calibrateImage", &PhotoCalib::calibrateImage, "maskedImage"_a,
            "includeScaleUncertainty"_a = true);

    cls.def("calibrateCatalog",
            py::overload_cast<afw::table::SourceCatalog const &, std::vector<std::string> const &>(
                    &PhotoCalib::calibrateCatalog, py::const_),
            "maskedImage"_a, "fluxFields"_a);
    cls.def("calibrateCatalog",
            py::overload_cast<afw::table::SourceCatalog const &>(&PhotoCalib::calibrateCatalog, py::const_),
            "maskedImage"_a);

    /* Operators */
    cls.def("__eq__", &PhotoCalib::operator==, py::is_operator());
    cls.def("__ne__", &PhotoCalib::operator!=, py::is_operator());

    /* Utility functions */
    mod.def("makePhotoCalibFromMetadata",
            py::overload_cast<daf::base::PropertySet &, bool>(makePhotoCalibFromMetadata), "metadata"_a,
            "strip"_a = false);
    mod.def("makePhotoCalibFromCalibZeroPoint",
            py::overload_cast<double, double>(makePhotoCalibFromCalibZeroPoint), "instFluxMag0"_a,
            "instFluxMag0Err"_a = false);

    utils::python::addOutputOp(cls, "__str__");
    cls.def("__repr__", [](PhotoCalib const &self) {
        std::ostringstream os;
        os << "PhotoCalib(" << self << ")";
        return os.str();
    });
}

template <typename PixelT>
static void declareImageSlice(py::module &mod, std::string const &suffix) {
    using Class = ImageSlice<PixelT>;

    py::class_<Class, std::shared_ptr<Class>, Image<PixelT>> cls(mod, ("ImageSlice" + suffix).c_str());

    cls.def(py::init<Image<PixelT> const &>(), "img"_a);

    py::enum_<typename Class::ImageSliceType>(cls, "ImageSliceType")
            .value("ROW", Class::ImageSliceType::ROW)
            .value("COLUMN", Class::ImageSliceType::COLUMN)
            .export_values();

    cls.def("getImageSliceType", &Class::getImageSliceType);

    cls.def("__add__", [](ImageSlice<PixelT> &self, Image<PixelT> const &other) { return self + other; },
            py::is_operator());
    cls.def("__mul__", [](ImageSlice<PixelT> &self, Image<PixelT> const &other) { return self * other; },
            py::is_operator());
    cls.def("__iadd__", [](ImageSlice<PixelT> &self, Image<PixelT> const &other) {
        self += other;
        return self;
    });
    cls.def("__imul__", [](ImageSlice<PixelT> &self, Image<PixelT> const &other) {
        self *= other;
        return self;
    });

    cls.def("__add__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self + other; },
            py::is_operator());
    cls.def("__sub__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self - other; },
            py::is_operator());
    cls.def("__mul__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self * other; },
            py::is_operator());
    cls.def("__truediv__",
            [](Image<PixelT> const &self, ImageSlice<PixelT> const &other) { return self / other; },
            py::is_operator());
    cls.def("__iadd__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self += other;
        return self;
    });
    cls.def("__isub__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self -= other;
        return self;
    });
    cls.def("__imul__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self *= other;
        return self;
    });
    cls.def("__itruediv__", [](Image<PixelT> &self, ImageSlice<PixelT> const &other) {
        self /= other;
        return self;
    });
}

WRAP(ImageSlice) {

    declareImageSlice<float>(mod, "F");
    declareImageSlice<double>(mod, "D");
}

template <typename ImageT>
static void declareImagePca(py::module &mod, std::string const &suffix) {
    py::class_<ImagePca<ImageT>, std::shared_ptr<ImagePca<ImageT>>> cls(mod, ("ImagePca" + suffix).c_str());

    cls.def(py::init<bool>(), "constantWeight"_a = true);

    cls.def("addImage", &ImagePca<ImageT>::addImage, "img"_a, "flux"_a = 0.0);
    cls.def("getImageList", &ImagePca<ImageT>::getImageList);
    cls.def("getDimensions", &ImagePca<ImageT>::getDimensions);
    cls.def("getMean", &ImagePca<ImageT>::getMean);
    cls.def("analyze", &ImagePca<ImageT>::analyze);
    cls.def("updateBadPixels", &ImagePca<ImageT>::updateBadPixels);
    cls.def("getEigenValues", &ImagePca<ImageT>::getEigenValues);
    cls.def("getEigenImages", &ImagePca<ImageT>::getEigenImages);
}

template <typename Image1T, typename Image2T>
static void declareInnerProduct(py::module &mod) {
    mod.def("innerProduct",
            (double (*)(Image1T const &, Image2T const &, int const))innerProduct<Image1T, Image2T>, "lhs"_a,
            "rhs"_a, "border"_a = 0);
}

using PyFilterLabel = py::class_<FilterLabel, std::shared_ptr<FilterLabel>, typehandling::Storable>;

// Macro to convert an exception without defining a global handler for it
#define _DELEGATE_EXCEPTION(call, cpp_ex, py_ex) \
    try {                                        \
        return call;                             \
    } catch (cpp_ex const& e) {                  \
        std::throw_with_nested(py_ex(e.what())); \
    }

PyFilterLabel declareFilter(py::module& mod) {
    // Include Python-only constructor in class doc, since it's what people should use
    auto initDoc = R"delim(
        Attributes
        ----------
        band : str, optional
            The band associated with this label.
        physical : str, optional
            The physical filter associated with this label.
    )delim";
    return PyFilterLabel(mod, "FilterLabel", initDoc);
}

void define(py::module& mod, PyFilterLabel& cls) {
    table::io::python::addPersistableMethods(cls);

    cls.def_static("fromBandPhysical", &FilterLabel::fromBandPhysical, "band"_a, "physical"_a);
    cls.def_static("fromBand", &FilterLabel::fromBand, "band"_a);
    cls.def_static("fromPhysical", &FilterLabel::fromPhysical, "physical"_a);

    // Keyword constructor
    /* This is messy in C++, but it's hard to write a Python __init__ that delegates to a factory,
     * and the pybind11 docs imply that this way is less prone to multiple-definition errors.
     * In C++17, we should be able to replace py::object with std::optional<string>.
     */
    cls.def(py::init([](py::object band, py::object physical) {
                try {
                    // Expand as we get more combinations of keywords
                    if (!band.is_none() && !physical.is_none()) {
                        return FilterLabel::fromBandPhysical(py::cast<std::string>(band),
                                                             py::cast<std::string>(physical));
                    } else if (!band.is_none()) {
                        return FilterLabel::fromBand(py::cast<std::string>(band));
                    } else if (!physical.is_none()) {
                        return FilterLabel::fromPhysical(py::cast<std::string>(physical));
                    } else {
                        throw py::value_error("Need at least one of band, physical");
                    }
                } catch (py::cast_error const& e) {
                    // By default cast_error is wrapped as RuntimeError
                    std::throw_with_nested(py::type_error(e.what()));
                }
            }),
            // TODO: use py::kw_only() in pybind11 2.6 or later (DM-27247)
            "band"_a = py::none(), "physical"_a = py::none());

    cls.def("hasBandLabel", &FilterLabel::hasBandLabel);
    cls.def_property_readonly("bandLabel", [](FilterLabel const& label) {
        _DELEGATE_EXCEPTION(label.getBandLabel(), pex::exceptions::LogicError, std::runtime_error);
    });
    cls.def("hasPhysicalLabel", &FilterLabel::hasPhysicalLabel);
    cls.def_property_readonly("physicalLabel", [](FilterLabel const& label) {
        _DELEGATE_EXCEPTION(label.getPhysicalLabel(), pex::exceptions::LogicError, std::runtime_error);
    });
    cls.def("__eq__", &FilterLabel::operator==, py::is_operator());
    cls.def("__ne__", &FilterLabel::operator!=, py::is_operator());

    cls.def("__repr__", &FilterLabel::toString);
    // Neither __copy__ nor __deepcopy__ default to each other
    cls.def("__copy__", [](const FilterLabel& obj) { return obj.cloneStorable(); });
    cls.def("__deepcopy__", [](const FilterLabel& obj, py::dict& memo) { return obj.cloneStorable(); });

    // Free functions

    mod.def("getDatabaseFilterLabel", &getDatabaseFilterLabel, "filterLabel"_a);
}

WRAP(FilterLabel) {
    // then declare classes
    auto cls = declareFilter(mod);
    // then import dependencies used in method signatures
    // none
    // and now we can safely define methods and other attributes
    define(mod, cls);
}

WRAP(ImagePca) {
    declareImagePca<Image<int>>(mod, "I");
    declareImagePca<Image<float>>(mod, "F");
    declareImagePca<Image<double>>(mod, "D");
    declareImagePca<Image<std::uint16_t>>(mod, "U");
    declareImagePca<Image<std::uint64_t>>(mod, "L");
    declareImagePca<MaskedImage<int>>(mod, "MI");
    declareImagePca<MaskedImage<float>>(mod, "MF");
    declareImagePca<MaskedImage<double>>(mod, "MD");
    declareImagePca<MaskedImage<std::uint16_t>>(mod, "MU");
    declareImagePca<MaskedImage<std::uint64_t>>(mod, "ML");

    declareInnerProduct<Image<int>, Image<int>>(mod);
    declareInnerProduct<Image<float>, Image<float>>(mod);
    declareInnerProduct<Image<double>, Image<double>>(mod);
    declareInnerProduct<Image<std::uint16_t>, Image<std::uint16_t>>(mod);
    declareInnerProduct<Image<std::uint64_t>, Image<std::uint64_t>>(mod);

    declareInnerProduct<Image<float>, Image<double>>(mod);
    declareInnerProduct<Image<double>, Image<float>>(mod);
}

using PyFilterProperty = py::class_<FilterProperty, std::shared_ptr<FilterProperty>>;

using PyFilter = py::class_<Filter, std::shared_ptr<Filter>, typehandling::Storable>;

WRAP(Filter) {

    mod.def("stripFilterKeywords", &detail::stripFilterKeywords, "metadata"_a);

    PyFilterProperty clsFilterProperty(mod, "FilterProperty");
    clsFilterProperty.def(py::init<std::string const &, double, double, double, bool>(),
                          "name"_a, "lambdaEff"_a, "lambdaMin"_a = NAN, "lambdaMax"_a = NAN,
                          "force"_a = false);
    // note: metadata should be defaulted with "metadata"_a=daf::base::PropertySet()
    // but that causes an error about copying when the Python extension is imported
    clsFilterProperty.def(py::init<std::string const &, daf::base::PropertySet const &, bool>(), "name"_a,
                          "metadata"_a, "force"_a = false);
    clsFilterProperty.def(
            "__eq__", [](FilterProperty const &self, FilterProperty const &other) { return self == other; },
            py::is_operator());
    clsFilterProperty.def(
            "__ne__", [](FilterProperty const &self, FilterProperty const &other) { return self != other; },
            py::is_operator());
    clsFilterProperty.def("getName", &FilterProperty::getName);
    clsFilterProperty.def("getLambdaEff", &FilterProperty::getLambdaEff);
    clsFilterProperty.def("getLambdaMin", &FilterProperty::getLambdaMin);
    clsFilterProperty.def("getLambdaMax", &FilterProperty::getLambdaMax);
    clsFilterProperty.def_static("reset", &FilterProperty::reset);
    clsFilterProperty.def_static("lookup", &FilterProperty::lookup, "name"_a);

    PyFilter clsFilter(mod, "Filter");
    clsFilter.def(py::init<std::string const &, bool const>(), "name"_a, "force"_a = false);
    clsFilter.def(py::init<int>(), "id"_a = Filter::UNKNOWN);
    clsFilter.def(py::init<std::shared_ptr<daf::base::PropertySet const>, bool const>(), "metadata"_a,
                  "force"_a = false);
    clsFilter.def("__eq__", [](Filter const &self, Filter const &other) { return self == other; },
                  py::is_operator());
    clsFilter.def("__ne__", [](Filter const &self, Filter const &other) { return self != other; },
                  py::is_operator());
    clsFilter.def_readonly_static("AUTO", &Filter::AUTO);
    clsFilter.def_readonly_static("UNKNOWN", &Filter::UNKNOWN);
    clsFilter.def("getId", &Filter::getId);
    clsFilter.def("getName", &Filter::getName);
    clsFilter.def("getCanonicalName", &Filter::getCanonicalName);
    clsFilter.def("getAliases", &Filter::getAliases);
    clsFilter.def("getFilterProperty", &Filter::getFilterProperty);
    clsFilter.def_static("reset", &Filter::reset);
    clsFilter.def_static("define", &Filter::define, "filterProperty"_a, "id"_a = Filter::AUTO,
                         "force"_a = false);
    clsFilter.def_static("defineAlias", &Filter::defineAlias, "oldName"_a, "newName"_a, "force"_a = false);
    clsFilter.def_static("getNames", &Filter::getNames);
}

using PyExposureInfo = py::class_<ExposureInfo, std::shared_ptr<ExposureInfo>>;

// Template methods where we can use pybind11's overload resolution (T is input)
template <class T>
void declareGenericMethods(PyExposureInfo &cls) {
    using Class = PyExposureInfo::type;
    cls.def("setComponent",
            [](PyExposureInfo::type &self, std::string const &key, T const &object) {
                self.setComponent(typehandling::makeKey<T>(key), object);
            },
            "key"_a, "object"_a);
}
// Template methods where we need to provide a unified interface (T is not input)
void declareGenericMethodsMerged(PyExposureInfo &cls) {
    using typehandling::Storable;
    using Class = PyExposureInfo::type;
    cls.def("hasComponent",
            [](Class const &self, std::string const &key) {
                return self.hasComponent(typehandling::makeKey<std::shared_ptr<Storable const>>(key));
            },
            "key"_a);
    cls.def("getComponent",
            [](Class const &self, std::string const &key) -> py::object {
                auto sharedKey = typehandling::makeKey<std::shared_ptr<Storable const>>(key);
                // Cascading if-elses to support other types in the future
                if (self.hasComponent(sharedKey)) {
                    return py::cast(self.getComponent(sharedKey));
                } else {
                    return py::none();
                }
            },
            "key"_a);
    cls.def("removeComponent",
            [](Class &self, std::string const &key) {
                self.removeComponent(typehandling::makeKey<std::shared_ptr<Storable const>>(key));
            },
            "key"_a);
}

WRAP(ExposureInfo) {

    /* Module level */
    PyExposureInfo cls(mod, "ExposureInfo");

    /* Member types and enums */

    /* Constructors */
    cls.def(py::init<std::shared_ptr<geom::SkyWcs const> const &,
                     std::shared_ptr<detection::Psf const> const &, std::shared_ptr<PhotoCalib const> const &,
                     std::shared_ptr<cameraGeom::Detector const> const &,
                     std::shared_ptr<geom::polygon::Polygon const> const &, Filter const &,
                     std::shared_ptr<daf::base::PropertySet> const &, std::shared_ptr<CoaddInputs> const &,
                     std::shared_ptr<ApCorrMap> const &, std::shared_ptr<VisitInfo const> const &,
                     std::shared_ptr<TransmissionCurve const> const &>(),
            "wcs"_a = std::shared_ptr<geom::SkyWcs const>(),
            "psf"_a = std::shared_ptr<detection::Psf const>(),
            "photoCalib"_a = std::shared_ptr<PhotoCalib const>(),
            "detector"_a = std::shared_ptr<cameraGeom::Detector const>(),
            "polygon"_a = std::shared_ptr<geom::polygon::Polygon const>(), "filter"_a = Filter(),
            "metadata"_a = std::shared_ptr<daf::base::PropertySet>(),
            "coaddInputs"_a = std::shared_ptr<CoaddInputs>(), "apCorrMap"_a = std::shared_ptr<ApCorrMap>(),
            "visitInfo"_a = std::shared_ptr<VisitInfo const>(), "transmissionCurve"_a = nullptr);
    cls.def(py::init<>());
    cls.def(py::init<ExposureInfo>(), "other"_a);
    cls.def(py::init<ExposureInfo, bool>(), "other"_a, "copyMetadata"_a);

    /* Members */
    cls.attr("KEY_WCS") = ExposureInfo::KEY_WCS.getId();
    cls.def("hasWcs", &ExposureInfo::hasWcs);
    cls.def("getWcs", (std::shared_ptr<geom::SkyWcs>(ExposureInfo::*)()) & ExposureInfo::getWcs);
    cls.def("setWcs", &ExposureInfo::setWcs, "wcs"_a);

    cls.attr("KEY_DETECTOR") = ExposureInfo::KEY_DETECTOR.getId();
    cls.def("hasDetector", &ExposureInfo::hasDetector);
    cls.def("getDetector", &ExposureInfo::getDetector);
    cls.def("setDetector",
            [](ExposureInfo &self, py::object detector) {
                if (detector.is(py::none())) {
                    self.setDetector(nullptr);
                } else {
                    self.setDetector(py::cast<std::shared_ptr<afw::cameraGeom::Detector>>(detector));
                }
            },
            "detector"_a);

    cls.def("getFilter", &ExposureInfo::getFilter);
    cls.def("setFilter", &ExposureInfo::setFilter, "filter"_a);

    declareGenericMethods<std::shared_ptr<typehandling::Storable const>>(cls);
    declareGenericMethodsMerged(cls);

    cls.attr("KEY_PHOTO_CALIB") = ExposureInfo::KEY_PHOTO_CALIB.getId();
    cls.def("hasPhotoCalib", &ExposureInfo::hasPhotoCalib);
    cls.def("getPhotoCalib", &ExposureInfo::getPhotoCalib);
    cls.def("setPhotoCalib", &ExposureInfo::setPhotoCalib, "photoCalib"_a);

    cls.def("getMetadata", &ExposureInfo::getMetadata);
    cls.def("setMetadata", &ExposureInfo::setMetadata, "metadata"_a);

    cls.attr("KEY_PSF") = ExposureInfo::KEY_PSF.getId();
    cls.def("hasPsf", &ExposureInfo::hasPsf);
    cls.def("getPsf", &ExposureInfo::getPsf);
    cls.def("setPsf",
            [](ExposureInfo &self, py::object psf) {
                if (psf.is(py::none())) {
                    self.setPsf(nullptr);
                } else {
                    self.setPsf(py::cast<std::shared_ptr<afw::detection::Psf>>(psf));
                }
            },
            "psf"_a);

    cls.attr("KEY_VALID_POLYGON") = ExposureInfo::KEY_VALID_POLYGON.getId();
    cls.def("hasValidPolygon", &ExposureInfo::hasValidPolygon);
    cls.def("getValidPolygon", &ExposureInfo::getValidPolygon);
    cls.def("setValidPolygon",
            [](ExposureInfo &self, py::object polygon) {
                if (polygon.is(py::none())) {
                    self.setValidPolygon(nullptr);
                } else {
                    self.setValidPolygon(py::cast<std::shared_ptr<afw::geom::polygon::Polygon>>(polygon));
                }
            },
            "polygon"_a);

    cls.attr("KEY_AP_CORR_MAP") = ExposureInfo::KEY_AP_CORR_MAP.getId();
    cls.def("hasApCorrMap", &ExposureInfo::hasApCorrMap);
    cls.def("getApCorrMap", (std::shared_ptr<ApCorrMap>(ExposureInfo::*)()) & ExposureInfo::getApCorrMap);
    cls.def("setApCorrMap", &ExposureInfo::setApCorrMap, "apCorrMap"_a);
    cls.def("initApCorrMap", &ExposureInfo::initApCorrMap);

    cls.attr("KEY_COADD_INPUTS") = ExposureInfo::KEY_COADD_INPUTS.getId();
    cls.def("hasCoaddInputs", &ExposureInfo::hasCoaddInputs);
    cls.def("getCoaddInputs", &ExposureInfo::getCoaddInputs);
    cls.def("setCoaddInputs", &ExposureInfo::setCoaddInputs, "coaddInputs"_a);

    cls.def("hasVisitInfo", &ExposureInfo::hasVisitInfo);
    cls.def("getVisitInfo", &ExposureInfo::getVisitInfo);
    cls.def("setVisitInfo", &ExposureInfo::setVisitInfo, "visitInfo"_a);

    cls.attr("KEY_TRANSMISSION_CURVE") = ExposureInfo::KEY_TRANSMISSION_CURVE.getId();
    cls.def("hasTransmissionCurve", &ExposureInfo::hasTransmissionCurve);
    cls.def("getTransmissionCurve", &ExposureInfo::getTransmissionCurve);
    cls.def("setTransmissionCurve", &ExposureInfo::setTransmissionCurve, "transmissionCurve"_a);
}

using PyDefectBase = py::class_<DefectBase, std::shared_ptr<DefectBase>>;

WRAP(Defect) {

    PyDefectBase cls(mod, "DefectBase");

    /* Constructors */
    cls.def(py::init<const lsst::geom::Box2I &>(), "bbox"_a);

    /* Members */
    cls.def("getBBox", &DefectBase::getBBox);
    cls.def("getX0", &DefectBase::getX0);
    cls.def("getX1", &DefectBase::getX1);
    cls.def("getY0", &DefectBase::getY0);
    cls.def("getY1", &DefectBase::getY1);
    cls.def("clip", &DefectBase::clip);
    cls.def("shift", (void (DefectBase::*)(int, int)) & DefectBase::shift, "dx"_a, "dy"_a);
    cls.def("shift", (void (DefectBase::*)(lsst::geom::Extent2I const &)) & DefectBase::shift, "d"_a);
}

using PyColor = py::class_<Color, std::shared_ptr<Color>>;

WRAP(Color) {
    /* Module level */
    PyColor cls(mod, "Color");

    /* Constructors */
    cls.def(py::init<double>(), "g_r"_a = std::numeric_limits<double>::quiet_NaN());

    /* Operators */
    cls.def("__eq__", [](Color const& self, Color const& other) { return self == other; }, py::is_operator());
    cls.def("__ne__", [](Color const& self, Color const& other) { return self != other; }, py::is_operator());

    /* Members */
    cls.def("isIndeterminate", &Color::isIndeterminate);
    cls.def("getLambdaEff", &Color::getLambdaEff, "filter"_a);
}

using PyCoaddInputs = py::class_<CoaddInputs, std::shared_ptr<CoaddInputs>, typehandling::Storable>;

WRAP(CoaddInputs) {
    /* Module level */

    PyCoaddInputs cls(mod, "CoaddInputs");

    /* Constructors */
    cls.def(py::init<>());
    cls.def(py::init<table::Schema const &, table::Schema const &>(), "visitSchema"_a, "ccdSchema"_a);
    cls.def(py::init<table::ExposureCatalog const &, table::ExposureCatalog const &>(), "visits"_a, "ccds"_a);

    table::io::python::addPersistableMethods<CoaddInputs>(cls);

    /* Members */
    cls.def_readwrite("visits", &CoaddInputs::visits);
    cls.def_readwrite("ccds", &CoaddInputs::ccds);
    cls.def("isPersistable", &CoaddInputs::isPersistable);
}

template <typename T>
void declareVectorOperations(py::module& mod) {
    typedef ndarray::Array<T, 1> Array;
    typedef ndarray::Array<T const, 1> ConstArray;
    mod.def("abMagFromFlux", (Array(*)(ConstArray const&)) & abMagFromFlux<T>, "flux"_a);
    mod.def("abMagErrFromFluxErr", (Array(*)(ConstArray const&, ConstArray const&)) & abMagErrFromFluxErr<T>,
            "fluxErr"_a, "flux"_a);
    mod.def("fluxFromABMag", (Array(*)(ConstArray const&)) & fluxFromABMag<T>, "mag"_a);
    mod.def("fluxErrFromABMagErr", (Array(*)(ConstArray const&, ConstArray const&)) & fluxErrFromABMagErr<T>,
            "magErr"_a, "mag"_a);
}

WRAP(Calib) {
    /* Module level */
    mod.def("abMagFromFlux", (double (*)(double)) & abMagFromFlux, "flux"_a);
    mod.def("abMagErrFromFluxErr", (double (*)(double, double)) & abMagErrFromFluxErr, "fluxErr"_a, "flux"_a);
    mod.def("fluxFromABMag", (double (*)(double)) & fluxFromABMag, "mag"_a);
    mod.def("fluxErrFromABMagErr", (double (*)(double, double)) & fluxErrFromABMagErr, "magErr"_a, "mag"_a);
    declareVectorOperations<float>(mod);
    declareVectorOperations<double>(mod);
}


WRAP(MaskedImage) {

    auto clsMaskedImageF = declareMaskedImage<float>(mod, "F");
    auto clsMaskedImageD = declareMaskedImage<double>(mod, "D");
    auto clsMaskedImageI = declareMaskedImage<int>(mod, "I");
    auto clsMaskedImageU = declareMaskedImage<std::uint16_t>(mod, "U");
    auto clsMaskedImageL = declareMaskedImage<std::uint64_t>(mod, "L");

    // Declare constructors for casting all exposure types to to float and double
    // (the only two types of casts that Python supports)
    declareCastConstructor<int, float>(clsMaskedImageF);
    declareCastConstructor<int, double>(clsMaskedImageD);
    declareCastConstructor<float, double>(clsMaskedImageD);
    declareCastConstructor<double, float>(clsMaskedImageF);
    declareCastConstructor<std::uint16_t, float>(clsMaskedImageF);
    declareCastConstructor<std::uint16_t, double>(clsMaskedImageD);
    declareCastConstructor<std::uint64_t, float>(clsMaskedImageF);
    declareCastConstructor<std::uint64_t, double>(clsMaskedImageD);

    /* Module level */
    declareMakeMaskedImage<int>(mod);
    declareMakeMaskedImage<float>(mod);
    declareMakeMaskedImage<double>(mod);
    declareMakeMaskedImage<std::uint16_t>(mod);
    declareMakeMaskedImage<std::uint64_t>(mod);

    declareImagesOverlap<int, int>(mod);
    declareImagesOverlap<int, float>(mod);
    declareImagesOverlap<int, double>(mod);
    declareImagesOverlap<int, std::uint16_t>(mod);
    declareImagesOverlap<int, std::uint64_t>(mod);

    declareImagesOverlap<float, int>(mod);
    declareImagesOverlap<float, float>(mod);
    declareImagesOverlap<float, double>(mod);
    declareImagesOverlap<float, std::uint16_t>(mod);
    declareImagesOverlap<float, std::uint64_t>(mod);

    declareImagesOverlap<double, int>(mod);
    declareImagesOverlap<double, float>(mod);
    declareImagesOverlap<double, double>(mod);
    declareImagesOverlap<double, std::uint16_t>(mod);
    declareImagesOverlap<double, std::uint64_t>(mod);

    declareImagesOverlap<std::uint16_t, int>(mod);
    declareImagesOverlap<std::uint16_t, float>(mod);
    declareImagesOverlap<std::uint16_t, double>(mod);
    declareImagesOverlap<std::uint16_t, std::uint16_t>(mod);
    declareImagesOverlap<std::uint16_t, std::uint64_t>(mod);

    declareImagesOverlap<std::uint64_t, int>(mod);
    declareImagesOverlap<std::uint64_t, float>(mod);
    declareImagesOverlap<std::uint64_t, double>(mod);
    declareImagesOverlap<std::uint64_t, std::uint16_t>(mod);
    declareImagesOverlap<std::uint64_t, std::uint64_t>(mod);
}

WRAP(VisitInfo) {

    /* Module level */
    py::class_<VisitInfo, std::shared_ptr<VisitInfo>, typehandling::Storable> cls(mod, "VisitInfo");

    /* Member types and enums */
    py::enum_<RotType>(mod, "RotType")
            .value("UNKNOWN", RotType::UNKNOWN)
            .value("SKY", RotType::SKY)
            .value("HORIZON", RotType::HORIZON)
            .value("MOUNT", RotType::MOUNT)
            .export_values();

    /* Constructors */
    cls.def(py::init<table::RecordId, double, double, daf::base::DateTime const &, double,
            lsst::geom::Angle const &, lsst::geom::SpherePoint const &, lsst::geom::SpherePoint const &, double,
            lsst::geom::Angle const &, RotType const &, coord::Observatory const &,
            coord::Weather const &>(),
            "exposureId"_a = 0, "exposureTime"_a = nan, "darkTime"_a = nan, "date"_a = daf::base::DateTime(),
            "ut1"_a = nan, "era"_a = nanAngle, "boresightRaDec"_a = lsst::geom::SpherePoint(nanAngle, nanAngle),
            "boresightAzAlt"_a = lsst::geom::SpherePoint(nanAngle, nanAngle), "boresightAirmass"_a = nan,
            "boresightRotAngle"_a = nanAngle, "rotType"_a = RotType::UNKNOWN,
            "observatory"_a = coord::Observatory(nanAngle, nanAngle, nan),
            "weather"_a = coord::Weather(nan, nan, nan));
    cls.def(py::init<daf::base::PropertySet const &>(), "metadata"_a);
    cls.def(py::init<VisitInfo const &>(), "visitInfo"_a);

    table::io::python::addPersistableMethods<VisitInfo>(cls);

    /* Operators */
    cls.def("__eq__", [](VisitInfo const &self, VisitInfo const &other) { return self == other; },
    py::is_operator());
    cls.def("__ne__", [](VisitInfo const &self, VisitInfo const &other) { return self != other; },
    py::is_operator());

    /* Members */
    cls.def("getExposureId", &VisitInfo::getExposureId);
    cls.def("getExposureTime", &VisitInfo::getExposureTime);
    cls.def("getDarkTime", &VisitInfo::getDarkTime);
    cls.def("getDate", &VisitInfo::getDate);
    cls.def("getUt1", &VisitInfo::getUt1);
    cls.def("getEra", &VisitInfo::getEra);
    cls.def("getBoresightRaDec", &VisitInfo::getBoresightRaDec);
    cls.def("getBoresightAzAlt", &VisitInfo::getBoresightAzAlt);
    cls.def("getBoresightAirmass", &VisitInfo::getBoresightAirmass);
    cls.def("getBoresightParAngle", &VisitInfo::getBoresightParAngle);
    cls.def("getBoresightRotAngle", &VisitInfo::getBoresightRotAngle);
    cls.def("getRotType", &VisitInfo::getRotType);
    cls.def("getObservatory", &VisitInfo::getObservatory);
    cls.def("getWeather", &VisitInfo::getWeather);
    cls.def("isPersistable", &VisitInfo::isPersistable);
    cls.def("getLocalEra", &VisitInfo::getLocalEra);
    cls.def("getBoresightHourAngle", &VisitInfo::getBoresightHourAngle);

    utils::python::addOutputOp(cls, "__str__");

    /* Free Functions */
    mod.def("setVisitInfoMetadata", &detail::setVisitInfoMetadata, "metadata"_a, "visitInfo"_a);
    mod.def("stripVisitInfoKeywords", &detail::stripVisitInfoKeywords, "metadata"_a);
}

WRAP(ImageUtils) {
    mod.def("indexToPosition", indexToPosition);
    mod.def("positionToIndex", (int (*)(double))positionToIndex);
    mod.def("positionToIndex", (std::pair<int, double>(*)(double const, bool))positionToIndex);
}

}  // namespace

namespace pixel {
namespace {

/**
@internal Declare a SinglePixel for a MaskedImage

(Note that SinglePixel for Image is just the pixel type)

@tparam  Image plane type
@param mod  pybind11 module
@param[in] name  Name of Python class, e.g. "SinglePixelI" if PixelT is `int`.
*/
template <typename PixelT>
void declareSinglePixel(py::module& mod, std::string const& name) {
    //mod.def("makeSinglePixel", &makeSinglePixel<PixelT, MaskPixel, VariancePixel>, "x"_a, "m"_a, "v"_a);

    py::class_<SinglePixel<PixelT, MaskPixel, VariancePixel>> cls(mod, name.c_str());

    cls.def(py::init<PixelT, MaskPixel, VariancePixel>(), "image"_a, "mask"_a = 0, "variance"_a = 0);
}

static double const nan(std::numeric_limits<double>::quiet_NaN());
static lsst::geom::Angle const nanAngle(nan);

}  // anonymous

WRAP(Pixel) {
    declareSinglePixel<float>(mod, "SinglePixelF");
    declareSinglePixel<double>(mod, "SinglePixelD");
    declareSinglePixel<int>(mod, "SinglePixelI");
    declareSinglePixel<std::uint16_t>(mod, "SinglePixelU");
    declareSinglePixel<std::uint64_t>(mod, "SinglePixelL");
}
} // pixel

WRAP(Image) {
    pixel::wrapPixel(mod);
    wrapImg(mod);
    wrapMaskedImage(mod);
    wrapApCorrMap(mod);
    wrapVisitInfo(mod);
    wrapTransmissionCurve(mod);
    wrapPhotocalib(mod);
    wrapImageSlice(mod);
    wrapFilterLabel(mod);
    wrapFilter(mod);
    wrapCoaddInputs(mod);
    wrapExposureInfo(mod);
    wrapExposure(mod);
    wrapReaders(mod);
    wrapDefect(mod);
    wrapColor(mod);
    wrapCalib(mod);
}
}
}
}  // namespace lsst::afw::image
