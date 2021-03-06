/*
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/geom/polygon/Polygon.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/cameraGeom/Detector.h"
#include "lsst/afw/image/ApCorrMap.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/TransmissionCurve.h"
#include "lsst/afw/image/ExposureFitsReader.h"

namespace lsst {
namespace afw {
namespace image {

namespace {


template <typename T, std::size_t N>
bool _contains(std::array<T, N> const& array, T const& value) {
    for (T const& element : array) {
        if (element == value) {
            return true;
        }
    }
    return false;
}

}  // namespace

class ExposureFitsReader::MetadataReader {
public:
    MetadataReader(std::shared_ptr<daf::base::PropertyList> primaryMetadata,
                   std::shared_ptr<daf::base::PropertyList> imageMetadata, lsst::geom::Point2I const& xy0) {
        auto versionName = ExposureInfo::getFitsSerializationVersionName();
        if (primaryMetadata->exists(versionName)) {
            version = primaryMetadata->getAsInt(versionName);
            primaryMetadata->remove(versionName);
        } else {
            version = 0;  // unversioned files are implicitly version 0
        }
        if (version > ExposureInfo::getFitsSerializationVersion()) {
            throw LSST_EXCEPT(pex::exceptions::TypeError,
                              str(boost::format("Cannot read Exposure FITS version >= %i") %
                                  ExposureInfo::getFitsSerializationVersion()));
        }

        // Try to read WCS from image metadata, and if found, strip the keywords used
        try {
            wcs = afw::geom::makeSkyWcs(*imageMetadata, true);
        } catch (lsst::pex::exceptions::TypeError const&) {
            //LOGLS_DEBUG(_log, "No WCS found in FITS metadata");
        }
        if (wcs && any(xy0.ne(lsst::geom::Point2I(0, 0)))) {
            wcs = wcs->copyAtShiftedPixelOrigin(lsst::geom::Extent2D(xy0));
        }

        // Strip LTV1, LTV2 from imageMetadata, because we don't use it internally
        imageMetadata->remove("LTV1");
        imageMetadata->remove("LTV2");

        if (!imageMetadata->exists("INHERIT")) {
            // New-style exposures put everything but the Wcs in the primary HDU, use
            // INHERIT keyword in the others.  For backwards compatibility, if we don't
            // find the INHERIT keyword, we ignore the primary HDU metadata and expect
            // everything to be in the image HDU metadata.  Note that we can't merge them,
            // because they're probably duplicates.
            metadata = imageMetadata;
        } else {
            metadata = primaryMetadata;
        }

        filter = Filter(metadata, true);
        detail::stripFilterKeywords(metadata);

        visitInfo = std::make_shared<VisitInfo>(*metadata);
        detail::stripVisitInfoKeywords(*metadata);

        // Version 0 persisted Calib FLUXMAG0 in the metadata, >=1 persisted PhotoCalib as a binary table.
        if (version == 0) {
            photoCalib = makePhotoCalibFromMetadata(*metadata, true);
        }

        // Strip MJD-OBS and DATE-OBS from metadata; those may be read by
        // either SkyWcs or VisitInfo or both, so neither can strip them.
        metadata->remove("MJD-OBS");
        metadata->remove("DATE-OBS");

        // Strip DETSER, DETNAME; these are added when writing an Exposure
        // with a Detector
        metadata->remove("DETNAME");
        metadata->remove("DETSER");
    }

    int version;
    std::shared_ptr<daf::base::PropertyList> metadata;
    Filter filter;
    std::shared_ptr<afw::geom::SkyWcs> wcs;
    std::shared_ptr<PhotoCalib> photoCalib;
    std::shared_ptr<VisitInfo> visitInfo;
};

class ExposureFitsReader::ArchiveReader {
public:
    enum Component {
        PSF = 0,
        WCS,
        COADD_INPUTS,
        AP_CORR_MAP,
        VALID_POLYGON,
        TRANSMISSION_CURVE,
        DETECTOR,
        PHOTOCALIB,
        N_ARCHIVE_COMPONENTS
    };

    explicit ArchiveReader(daf::base::PropertyList& metadata) {
        auto popInt = [&metadata](std::string const& name) {
            // The default of zero will cause archive.get to return a
            // null/empty pointer, just as if a null/empty pointer was
            // originally written to the archive.
            int r = 0;
            if (metadata.exists(name)) {
                r = metadata.get<int>(name);
                // We remove metadata entries to maintaing our practice
                // of stripped metadata entries that have been used to
                // construct more structured components.
                metadata.remove(name);
            }
            return r;
        };
        _hdu = popInt("AR_HDU");
        if (_hdu == 0) {
            _state = ArchiveState::MISSING;
        } else {
            --_hdu;  // Switch from FITS 1-indexed convention to LSST 0-indexed convention.
            _state = ArchiveState::PRESENT;
        }
        // Read in traditional components using old-style IDs, for backwards compatibility
        _ids[PSF] = popInt("PSF_ID");
        _ids[WCS] = popInt("SKYWCS_ID");
        _ids[COADD_INPUTS] = popInt("COADD_INPUTS_ID");
        _ids[AP_CORR_MAP] = popInt("AP_CORR_MAP_ID");
        _ids[VALID_POLYGON] = popInt("VALID_POLYGON_ID");
        _ids[TRANSMISSION_CURVE] = popInt("TRANSMISSION_CURVE_ID");
        _ids[DETECTOR] = popInt("DETECTOR_ID");
        _ids[PHOTOCALIB] = popInt("PHOTOCALIB_ID");

        // "Extra" components use a different keyword convention to avoid collisions with non-persistence IDs
        std::vector<std::string> toStrip;
        for (std::string const& headerKey : metadata) {
            static std::string const PREFIX = "ARCHIVE_ID_";
            if (headerKey.substr(0, PREFIX.size()) == PREFIX) {
                std::string componentName = headerKey.substr(PREFIX.size());
                int archiveId = metadata.get<int>(headerKey);
                if (!_contains(_ids, archiveId)) {
                    _extraIds.emplace(componentName, archiveId);
                }
                toStrip.push_back(headerKey);
                toStrip.push_back(componentName + "_ID");  // strip corresponding old-style ID, if it exists
            }
        }
        for (std::string const& key : toStrip) {
            metadata.remove(key);
        }
    }

    /**
     * Read a known component, if available.
     *
     * @param fitsFile The file from which to read the component. Must match
     *                 the metadata used to construct this object.
     * @param c The component to read. Must be convertible to ``T``.
     *
     * @return The desired component, or ``nullptr`` if the file could not be read.
     */
    template <typename T>
    std::shared_ptr<T> readComponent(afw::fits::Fits* fitsFile, Component c) {
        if (!_ensureLoaded(fitsFile)) {
            return nullptr;
        }
        return _archive.get<T>(_ids[c]);
    }

    /**
     * Read the components that are stored using arbitrary-component support.
     *
     * @param fitsFile The file from which to read the components. Must match
     *                 the metadata used to construct this object.
     *
     * @return a map from string IDs to components, or an empty map if the
     *         file could not be read.
     */
    std::map<std::string, std::shared_ptr<table::io::Persistable>> readExtraComponents(
            afw::fits::Fits* fitsFile) {
        std::map<std::string, std::shared_ptr<table::io::Persistable>> result;

        if (!_ensureLoaded(fitsFile)) {
            return result;
        }

        // Not safe to call getAll if a component cannot be unpersisted
        // Instead, look for the archives registered in the metadata
        for (auto const& keyValue : _extraIds) {
            std::string componentName = keyValue.first;
            int archiveId = keyValue.second;

            try {
                result.emplace(componentName, _archive.get(archiveId));
            } catch (pex::exceptions::NotFoundError const& err) {
                //LOGLS_WARN(_log,
                //           "Could not read component " << componentName << "; skipping: " << err.what());
            }
        }
        return result;
    }

private:
    bool _ensureLoaded(afw::fits::Fits* fitsFile) {
        if (_state == ArchiveState::MISSING) {
            return false;
        }
        if (_state == ArchiveState::PRESENT) {
            afw::fits::HduMoveGuard guard(*fitsFile, _hdu);
            _archive = table::io::InputArchive::readFits(*fitsFile);
            _state = ArchiveState::LOADED;
        }
        assert(_state == ArchiveState::LOADED);  // constructor body should guarantee it's not UNKNOWN
        return true;
    }

    enum class ArchiveState { UNKNOWN, MISSING, PRESENT, LOADED };

    int _hdu = 0;
    ArchiveState _state = ArchiveState::UNKNOWN;
    table::io::InputArchive _archive;
    std::array<int, N_ARCHIVE_COMPONENTS> _ids = {0};
    std::map<std::string, int> _extraIds;
};

ExposureFitsReader::ExposureFitsReader(std::string const& fileName) : _maskedImageReader(fileName) {}

ExposureFitsReader::ExposureFitsReader(fits::MemFileManager& manager) : _maskedImageReader(manager) {}

ExposureFitsReader::ExposureFitsReader(fits::Fits* fitsFile) : _maskedImageReader(fitsFile) {}

ExposureFitsReader::~ExposureFitsReader() noexcept = default;

lsst::geom::Box2I ExposureFitsReader::readBBox(ImageOrigin origin) {
    return _maskedImageReader.readBBox(origin);
}

lsst::geom::Point2I ExposureFitsReader::readXY0(lsst::geom::Box2I const& bbox, ImageOrigin origin) {
    return _maskedImageReader.readXY0(bbox, origin);
}

std::string ExposureFitsReader::readImageDType() const { return _maskedImageReader.readImageDType(); }

std::string ExposureFitsReader::readMaskDType() const { return _maskedImageReader.readMaskDType(); }

std::string ExposureFitsReader::readVarianceDType() const { return _maskedImageReader.readVarianceDType(); }

std::shared_ptr<daf::base::PropertyList> ExposureFitsReader::readMetadata() {
    _ensureReaders();
    return _metadataReader->metadata;
}

std::shared_ptr<afw::geom::SkyWcs> ExposureFitsReader::readWcs() {
    _ensureReaders();
    auto r = _archiveReader->readComponent<afw::geom::SkyWcs>(_getFitsFile(), ArchiveReader::WCS);
    if (!r) {
        r = _metadataReader->wcs;
    }
    return r;
}

Filter ExposureFitsReader::readFilter() {
    _ensureReaders();
    return _metadataReader->filter;
}

std::shared_ptr<PhotoCalib> ExposureFitsReader::readPhotoCalib() {
    _ensureReaders();
    if (_metadataReader->version == 0) {
        return _metadataReader->photoCalib;
    } else {
        return _archiveReader->readComponent<image::PhotoCalib>(_getFitsFile(), ArchiveReader::PHOTOCALIB);
    }
}

std::shared_ptr<detection::Psf> ExposureFitsReader::readPsf() {
    _ensureReaders();
    return _archiveReader->readComponent<detection::Psf>(_getFitsFile(), ArchiveReader::PSF);
}

std::shared_ptr<afw::geom::polygon::Polygon> ExposureFitsReader::readValidPolygon() {
    _ensureReaders();
    return _archiveReader->readComponent<afw::geom::polygon::Polygon>(_getFitsFile(),
                                                                      ArchiveReader::VALID_POLYGON);
}

std::shared_ptr<ApCorrMap> ExposureFitsReader::readApCorrMap() {
    _ensureReaders();
    return _archiveReader->readComponent<ApCorrMap>(_getFitsFile(), ArchiveReader::AP_CORR_MAP);
}

std::shared_ptr<CoaddInputs> ExposureFitsReader::readCoaddInputs() {
    _ensureReaders();
    return _archiveReader->readComponent<CoaddInputs>(_getFitsFile(), ArchiveReader::COADD_INPUTS);
}

std::shared_ptr<VisitInfo> ExposureFitsReader::readVisitInfo() {
    _ensureReaders();
    return _metadataReader->visitInfo;
}

std::shared_ptr<TransmissionCurve> ExposureFitsReader::readTransmissionCurve() {
    _ensureReaders();
    return _archiveReader->readComponent<TransmissionCurve>(_getFitsFile(),
                                                            ArchiveReader::TRANSMISSION_CURVE);
}

std::shared_ptr<cameraGeom::Detector> ExposureFitsReader::readDetector() {
    _ensureReaders();
    return _archiveReader->readComponent<cameraGeom::Detector>(_getFitsFile(), ArchiveReader::DETECTOR);
}

std::map<std::string, std::shared_ptr<table::io::Persistable>> ExposureFitsReader::readExtraComponents() {
    _ensureReaders();
    return _archiveReader->readExtraComponents(_getFitsFile());
}

std::shared_ptr<ExposureInfo> ExposureFitsReader::readExposureInfo() {
    auto result = std::make_shared<ExposureInfo>();
    result->setMetadata(readMetadata());
    result->setFilter(readFilter());
    result->setPhotoCalib(readPhotoCalib());
    result->setVisitInfo(readVisitInfo());
    // When reading an ExposureInfo (as opposed to reading individual
    // components), we warn and try to proceed when a component is present
    // but can't be read due its serialization factory not being set up
    // (that's what throws the NotFoundErrors caught below).
    try {
        result->setPsf(readPsf());
    } catch (pex::exceptions::NotFoundError& err) {
        //LOGLS_WARN(_log, "Could not read PSF; setting to null: " << err.what());
    }
    try {
        result->setCoaddInputs(readCoaddInputs());
    } catch (pex::exceptions::NotFoundError& err) {
        //LOGLS_WARN(_log, "Could not read CoaddInputs; setting to null: " << err.what());
    }
    try {
        result->setApCorrMap(readApCorrMap());
    } catch (pex::exceptions::NotFoundError& err) {
        //LOGLS_WARN(_log, "Could not read ApCorrMap; setting to null: " << err.what());
    }
    try {
        result->setValidPolygon(readValidPolygon());
    } catch (pex::exceptions::NotFoundError& err) {
        //LOGLS_WARN(_log, "Could not read ValidPolygon; setting to null: " << err.what());
    }
    try {
        result->setTransmissionCurve(readTransmissionCurve());
    } catch (pex::exceptions::NotFoundError& err) {
        //LOGLS_WARN(_log, "Could not read TransmissionCurve; setting to null: " << err.what());
    }
    try {
        result->setDetector(readDetector());
    } catch (pex::exceptions::NotFoundError& err) {
        //LOGLS_WARN(_log, "Could not read Detector; setting to null: " << err.what());
    }
    // In the case of WCS, we fall back to the metadata WCS if the one from
    // the archive can't be read.
    _ensureReaders();
    result->setWcs(_metadataReader->wcs);
    try {
        auto wcs = _archiveReader->readComponent<afw::geom::SkyWcs>(_getFitsFile(), ArchiveReader::WCS);
        if (!wcs) {
            //LOGLS_DEBUG(_log, "No WCS found in binary table");
        } else {
            result->setWcs(wcs);
        }
    } catch (pex::exceptions::NotFoundError& err) {
        auto msg = str(boost::format("Could not read WCS extension; setting to null: %s") % err.what());
        if (result->hasWcs()) {
            msg += " ; using WCS from FITS header";
        }
        //LOGLS_WARN(_log, msg);
    }
    for (auto keyValue : readExtraComponents()) {
        using StorablePtr = std::shared_ptr<typehandling::Storable const>;
        std::string key = keyValue.first;
        StorablePtr object = std::dynamic_pointer_cast<StorablePtr::element_type>(keyValue.second);

        if (object.use_count() > 0) {  // Failed cast guarantees empty pointer, but not a null one
            result->setComponent(typehandling::makeKey<StorablePtr>(key), object);
        } else {
            //LOGLS_WARN(_log, "Data corruption: generic component " << key << " is not a Storable; skipping.");
        }
    }
    return result;
}  // namespace image

template <typename ImagePixelT>
Image<ImagePixelT> ExposureFitsReader::readImage(lsst::geom::Box2I const& bbox, ImageOrigin origin,
                                                 bool allowUnsafe) {
    return _maskedImageReader.readImage<ImagePixelT>(bbox, origin, allowUnsafe);
}

template <typename ImagePixelT>
ndarray::Array<ImagePixelT, 2, 2> ExposureFitsReader::readImageArray(lsst::geom::Box2I const& bbox,
                                                                     ImageOrigin origin, bool allowUnsafe) {
    return _maskedImageReader.readImageArray<ImagePixelT>(bbox, origin, allowUnsafe);
}

template <typename MaskPixelT>
Mask<MaskPixelT> ExposureFitsReader::readMask(lsst::geom::Box2I const& bbox, ImageOrigin origin,
                                              bool conformMasks, bool allowUnsafe) {
    return _maskedImageReader.readMask<MaskPixelT>(bbox, origin, conformMasks, allowUnsafe);
}

template <typename MaskPixelT>
ndarray::Array<MaskPixelT, 2, 2> ExposureFitsReader::readMaskArray(lsst::geom::Box2I const& bbox,
                                                                   ImageOrigin origin, bool allowUnsafe) {
    return _maskedImageReader.readMaskArray<MaskPixelT>(bbox, origin, allowUnsafe);
}

template <typename VariancePixelT>
Image<VariancePixelT> ExposureFitsReader::readVariance(lsst::geom::Box2I const& bbox, ImageOrigin origin,
                                                       bool allowUnsafe) {
    return _maskedImageReader.readVariance<VariancePixelT>(bbox, origin, allowUnsafe);
}

template <typename VariancePixelT>
ndarray::Array<VariancePixelT, 2, 2> ExposureFitsReader::readVarianceArray(lsst::geom::Box2I const& bbox,
                                                                           ImageOrigin origin,
                                                                           bool allowUnsafe) {
    return _maskedImageReader.readVarianceArray<VariancePixelT>(bbox, origin, allowUnsafe);
}

template <typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> ExposureFitsReader::readMaskedImage(
        lsst::geom::Box2I const& bbox, ImageOrigin origin, bool conformMasks, bool allowUnsafe) {
    return _maskedImageReader.read<ImagePixelT, MaskPixelT, VariancePixelT>(bbox, origin, conformMasks,
                                                                            /* needAllHdus= */ false,
                                                                            allowUnsafe);
}

template <typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
Exposure<ImagePixelT, MaskPixelT, VariancePixelT> ExposureFitsReader::read(lsst::geom::Box2I const& bbox,
                                                                           ImageOrigin origin,
                                                                           bool conformMasks,
                                                                           bool allowUnsafe) {
    auto mi =
            readMaskedImage<ImagePixelT, MaskPixelT, VariancePixelT>(bbox, origin, conformMasks, allowUnsafe);
    return Exposure<ImagePixelT, MaskPixelT, VariancePixelT>(mi, readExposureInfo());
}

void ExposureFitsReader::_ensureReaders() {
    if (!_metadataReader) {
        auto metadataReader = std::make_unique<MetadataReader>(_maskedImageReader.readPrimaryMetadata(),
                                                               _maskedImageReader.readImageMetadata(),
                                                               _maskedImageReader.readXY0());
        _archiveReader = std::make_unique<ArchiveReader>(*metadataReader->metadata);
        _metadataReader = std::move(metadataReader);  // deferred for exception safety
    }
    assert(_archiveReader);  // should always be initialized with _metadataReader.
}

#define INSTANTIATE(ImagePixelT)                                                                            \
    template Exposure<ImagePixelT, MaskPixel, VariancePixel> ExposureFitsReader::read(                      \
            lsst::geom::Box2I const&, ImageOrigin, bool, bool);                                             \
    template Image<ImagePixelT> ExposureFitsReader::readImage(lsst::geom::Box2I const&, ImageOrigin, bool); \
    template ndarray::Array<ImagePixelT, 2, 2> ExposureFitsReader::readImageArray(lsst::geom::Box2I const&, \
                                                                                  ImageOrigin, bool);       \
    template MaskedImage<ImagePixelT, MaskPixel, VariancePixel> ExposureFitsReader::readMaskedImage(        \
            lsst::geom::Box2I const&, ImageOrigin, bool, bool)

INSTANTIATE(std::uint16_t);
INSTANTIATE(int);
INSTANTIATE(float);
INSTANTIATE(double);
INSTANTIATE(std::uint64_t);

template Mask<MaskPixel> ExposureFitsReader::readMask(lsst::geom::Box2I const&, ImageOrigin, bool, bool);
template ndarray::Array<MaskPixel, 2, 2> ExposureFitsReader::readMaskArray(lsst::geom::Box2I const&,
                                                                           ImageOrigin, bool);

template Image<VariancePixel> ExposureFitsReader::readVariance(lsst::geom::Box2I const&, ImageOrigin, bool);
template ndarray::Array<VariancePixel, 2, 2> ExposureFitsReader::readVarianceArray(lsst::geom::Box2I const&,
                                                                                   ImageOrigin, bool);

}  // namespace image
}  // namespace afw
}  // namespace lsst
