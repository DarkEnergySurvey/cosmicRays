
//namespace posix {  // here so no-one includes them first outside namespace posix {}
#include "pybind/afw_bind.h"
#include <unistd.h>
#include <fcntl.h>
//}

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cctype>

#include "boost/format.hpp"
#include "lsst/afw/detection.h"
#include "lsst/pex/exceptions.h"
#include "boost/any.hpp"
#include "lsst/afw/fits.h"

#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/geom/SkyWcs.h"

#define FITS_SIZE 2880

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace afw {
namespace display {

namespace  {
template <typename ImageT>
void replaceSaturatedPixels(
        ImageT& rim,                       //< R image (e.g. i)
        ImageT& gim,                       //< G image (e.g. r)
        ImageT& bim,                       //< B image (e.g. g)
        int borderWidth = 2,               //< width of border used to estimate colour of saturated regions
        float saturatedPixelValue = 65535  //< the brightness of a saturated pixel, once fixed
        );

/**
 * Calculate an IRAF/ds9-style zscaling.
 *
 * To quote Frank Valdes (http://iraf.net/forum/viewtopic.php?showtopic=134139)
 * <blockquote>
ZSCALE ALGORITHM

    The zscale algorithm is designed to display the  image  values  near
    the  median  image  value  without  the  time  consuming  process of
    computing a full image histogram.  This is particularly  useful  for
    astronomical  images  which  generally  have a very peaked histogram
    corresponding to  the  background  sky  in  direct  imaging  or  the
    continuum in a two dimensional spectrum.

    The  sample  of pixels, specified by values greater than zero in the
    sample mask zmask or by an  image  section,  is  selected  up  to  a
    maximum  of nsample pixels.  If a bad pixel mask is specified by the
    bpmask parameter then any pixels with mask values which are  greater
    than  zero  are not counted in the sample.  Only the first pixels up
    to the limit are selected where the order is by line beginning  from
    the  first line.  If no mask is specified then a grid of pixels with
    even spacing along lines and columns that  make  up  a  number  less
    than or equal to the maximum sample size is used.

    If  a  contrast of zero is specified (or the zrange flag is used and
    the image does not have a  valid  minimum/maximum  value)  then  the
    minimum  and maximum of the sample is used for the intensity mapping
    range.

    If the contrast  is  not  zero  the  sample  pixels  are  ranked  in
    brightness  to  form  the  function  I(i) where i is the rank of the
    pixel and I is its value.  Generally the midpoint of  this  function
    (the  median) is very near the peak of the image histogram and there
    is a well defined slope about the midpoint which is related  to  the
    width  of the histogram.  At the ends of the I(i) function there are
    a few very bright and dark pixels due to objects and defects in  the
    field.   To  determine  the  slope  a  linear  function  is fit with
    iterative rejection;

    <code>
            I(i) = intercept + slope * (i - midpoint)
    </code>

    If more than half of the points are rejected then there is  no  well
    defined  slope  and  the full range of the sample defines z1 and z2.
    Otherwise the endpoints of the linear function  are  used  (provided
    they are within the original range of the sample):

    <code>
            z1 = I(midpoint) + (slope / contrast) * (1 - midpoint)
            z2 = I(midpoint) + (slope / contrast) * (npoints - midpoint)
    </code>

    As  can  be  seen,  the parameter contrast may be used to adjust the
    contrast produced by this algorithm.
 * </blockquote>
 */
template <class T>
std::pair<double, double> getZScale(image::Image<T> const& image,  ///< The image we wish to stretch
                                    int const nSamples = 1000,     ///< Number of samples to use
                                    double const contrast = 0.25   ///< Stretch parameter; see description
                                    );

WRAP(Rgb) {
    /* Module level */
    mod.def("replaceSaturatedPixels", replaceSaturatedPixels<lsst::afw::image::MaskedImage<float>>, "rim"_a,
            "gim"_a, "bim"_a, "borderWidth"_a = 2, "saturatedPixelValue"_a = 65535);
    mod.def("getZScale", getZScale<std::uint16_t>, "image"_a, "nsamples"_a = 1000, "contrast"_a = 0.25);
    mod.def("getZScale", getZScale<float>, "image"_a, "nsamples"_a = 1000, "contrast"_a = 0.25);

}

template <typename ImageT>
class SetPixels {
public:
    explicit SetPixels() : _value(0) {}

    void setValue(float value) { _value = value; }

    void operator()(lsst::geom::Point2I const& point, typename ImageT::Pixel& arrayInput) { arrayInput = _value; }

private:
    float _value;
};

template <typename ImageT>
void replaceSaturatedPixels(ImageT& rim,      // R image (e.g. i)
                            ImageT& gim,      // G image (e.g. r)
                            ImageT& bim,      // B image (e.g. g)
                            int borderWidth,  // width of border used to estimate colour of saturated regions
                            float saturatedPixelValue  // the brightness of a saturated pixel, once fixed
                            ) {
    int const width = rim.getWidth(), height = rim.getHeight();
    int const x0 = rim.getX0(), y0 = rim.getY0();

    if (width != gim.getWidth() || height != gim.getHeight() || x0 != gim.getX0() || y0 != gim.getY0()) {
        throw LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                str(boost::format("R image has different size/origin from G image "
                                  "(%dx%d+%d+%d v. %dx%d+%d+%d") %
                    width % height % x0 % y0 % gim.getWidth() % gim.getHeight() % gim.getX0() % gim.getY0()));
    }
    if (width != bim.getWidth() || height != bim.getHeight() || x0 != bim.getX0() || y0 != bim.getY0()) {
        throw LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                str(boost::format("R image has different size/origin from B image "
                                  "(%dx%d+%d+%d v. %dx%d+%d+%d") %
                    width % height % x0 % y0 % bim.getWidth() % bim.getHeight() % bim.getX0() % bim.getY0()));
    }

    bool const useMaxPixel = !std::isfinite(saturatedPixelValue);

    SetPixels<typename ImageT::Image> setR, setG, setB;  // functors used to set pixel values

    // Find all the saturated pixels in any of the three image
    int const npixMin = 1;  // minimum number of pixels in an object
    afw::image::MaskPixel const SAT = rim.getMask()->getPlaneBitMask("SAT");
    detection::Threshold const satThresh(SAT, detection::Threshold::BITMASK);

    detection::FootprintSet sat(*rim.getMask(), satThresh, npixMin);
    sat.merge(detection::FootprintSet(*gim.getMask(), satThresh, npixMin));
    sat.merge(detection::FootprintSet(*bim.getMask(), satThresh, npixMin));
    // go through the list of saturated regions, determining the mean colour of the surrounding pixels
    typedef detection::FootprintSet::FootprintList FootprintList;
    std::shared_ptr<FootprintList> feet = sat.getFootprints();
    for (FootprintList::const_iterator ptr = feet->begin(), end = feet->end(); ptr != end; ++ptr) {
        std::shared_ptr<detection::Footprint> const foot = *ptr;
        auto const bigFoot = std::make_shared<detection::Footprint>(foot->getSpans()->dilated(borderWidth),
                                                                    foot->getRegion());

        double sumR = 0, sumG = 0, sumB = 0;  // sum of all non-saturated adjoining pixels
        double maxR = 0, maxG = 0, maxB = 0;  // maximum of non-saturated adjoining pixels

        for (auto span = bigFoot->getSpans()->begin(), send = bigFoot->getSpans()->end(); span != send;
             ++span) {
            int const y = span->getY() - y0;
            if (y < 0 || y >= height) {
                continue;
            }
            int sx0 = span->getX0() - x0;
            if (sx0 < 0) {
                sx0 = 0;
            }
            int sx1 = span->getX1() - x0;
            if (sx1 >= width) {
                sx1 = width - 1;
            }

            for (typename ImageT::iterator rptr = rim.at(sx0, y), rend = rim.at(sx1 + 1, y),
                                           gptr = gim.at(sx0, y), bptr = bim.at(sx0, y);
                 rptr != rend; ++rptr, ++gptr, ++bptr) {
                if (!((rptr.mask() | gptr.mask() | bptr.mask()) & SAT)) {
                    float val = rptr.image();
                    sumR += val;
                    if (val > maxR) {
                        maxR = val;
                    }

                    val = gptr.image();
                    sumG += val;
                    if (val > maxG) {
                        maxG = val;
                    }

                    val = bptr.image();
                    sumB += val;
                    if (val > maxB) {
                        maxB = val;
                    }
                }
            }
        }
        // OK, we have the mean fluxes for the pixels surrounding this set of saturated pixels
        // so we can figure out the proper values to use for the saturated ones
        float R = 0, G = 0, B = 0;  // mean intensities
        if (sumR + sumB + sumG > 0) {
            if (sumR > sumG) {
                if (sumR > sumB) {
                    R = useMaxPixel ? maxR : saturatedPixelValue;

                    G = (R * sumG) / sumR;
                    B = (R * sumB) / sumR;
                } else {
                    B = useMaxPixel ? maxB : saturatedPixelValue;
                    R = (B * sumR) / sumB;
                    G = (B * sumG) / sumB;
                }
            } else {
                if (sumG > sumB) {
                    G = useMaxPixel ? maxG : saturatedPixelValue;
                    R = (G * sumR) / sumG;
                    B = (G * sumB) / sumG;
                } else {
                    B = useMaxPixel ? maxB : saturatedPixelValue;
                    R = (B * sumR) / sumB;
                    G = (B * sumG) / sumB;
                }
            }
        }
        // Now that we know R, G, and B we can fix the values
        setR.setValue(R);
        foot->getSpans()->applyFunctor(setR, *rim.getImage());
        setG.setValue(G);
        foot->getSpans()->applyFunctor(setG, *gim.getImage());
        setB.setValue(B);
        foot->getSpans()->applyFunctor(setB, *bim.getImage());
    }
}

template void replaceSaturatedPixels(image::MaskedImage<float>& rim, image::MaskedImage<float>& gim,
                                     image::MaskedImage<float>& bim, int borderWidth,
                                     float saturatedPixelValue);

template <class T>
static void getSample(image::Image<T> const& image, std::size_t const nSamples, std::vector<T>& vSample) {
    int const width = image.getWidth();
    int const height = image.getHeight();

    // extract, from image, about nSample samples
    // such that they form a grid.
    vSample.reserve(nSamples);

    int const initialStride =
            std::max(1, static_cast<int>(std::sqrt(width * height / static_cast<double>(nSamples))));

    for (int stride = initialStride; stride >= 1; --stride) {
        vSample.clear();

        for (int y = 0; y < height; y += stride) {
            for (int x = 0; x < width; x += stride) {
                T const elem = image(x, y);
                if (std::isfinite(elem)) {
                    vSample.push_back(elem);
                }
            }
        }

        // if more than 80% of nSamples were sampled, OK.
        if (5 * vSample.size() > 4 * nSamples) {
            break;
        }
    }
}

static inline double computeSigma(std::vector<double> const& vFlat, std::vector<int> const& vBadPix,
                                  std::size_t const nGoodPix) {
    if (nGoodPix <= 1) {
        return 0;
    }

    double sumz = 0, sumsq = 0;
    for (unsigned i = 0; i < vFlat.size(); ++i) {
        if (!vBadPix[i]) {
            double const z = vFlat[i];

            sumz += z;
            sumsq += z * z;
        }
    }

    double const goodPix = nGoodPix;
    double const tmp = sumsq / (goodPix - 1) - sumz * sumz / (goodPix * (goodPix - 1));

    return (tmp > 0) ? std::sqrt(tmp) : 0;
}

template <class T>
static std::pair<double, double> fitLine(
        int* nGoodPixOut,  // returned; it'd be nice to use std::tuple from C++11
        std::vector<T> const& vSample, double const nSigmaClip, int const nGrow, int const minpix,
        int const nIter) {
    // map the indices of vSample to [-1.0, 1.0]
    double const xscale = 2.0 / (vSample.size() - 1);
    std::vector<double> xnorm;
    xnorm.reserve(vSample.size());
    for (std::size_t i = 0; i < vSample.size(); ++i) {
        xnorm.push_back(i * xscale - 1.0);
    }

    // Mask that is used in k-sigma clipping
    std::vector<int> vBadPix(vSample.size(), 0);

    int nGoodPix = vSample.size();
    int nGoodPixOld = nGoodPix + 1;

    // values to be obtained
    double intercept = 0;
    double slope = 0;

    for (int iteration = 0; iteration < nIter; ++iteration) {
        if (nGoodPix < minpix || nGoodPix >= nGoodPixOld) {
            break;
        }

        double sum = nGoodPix;
        double sumx = 0, sumy = 0, sumxx = 0, sumxy = 0;
        for (std::size_t i = 0; i < vSample.size(); ++i) {
            if (!vBadPix[i]) {
                double const x = xnorm[i];
                double const y = vSample[i];

                sumx += x;
                sumy += y;
                sumxx += x * x;
                sumxy += x * y;
            }
        }

        double delta = sum * sumxx - sumx * sumx;

        // slope and intercept
        intercept = (sumxx * sumy - sumx * sumxy) / delta;
        slope = (sum * sumxy - sumx * sumy) / delta;

        // residue
        std::vector<double> vFlat;
        vFlat.reserve(vSample.size());
        for (unsigned i = 0; i < vSample.size(); ++i) {
            vFlat.push_back(vSample[i] - (xnorm[i] * slope + intercept));
        }

        // Threshold of k-sigma clipping
        double const sigma = computeSigma(vFlat, vBadPix, nGoodPix);
        double const hcut = sigma * nSigmaClip;
        double const lcut = -hcut;

        // revise vBadPix
        nGoodPixOld = nGoodPix;
        nGoodPix = 0;

        for (unsigned i = 0; i < vSample.size(); ++i) {
            double val = vFlat[i];
            if (val < lcut || hcut < val) {
                vBadPix[i] = 1;
            }
        }

        // blurr vBadPix
        std::vector<int> vBadPix_new;
        vBadPix_new.reserve(vSample.size());
        for (unsigned x = 0; x < vSample.size(); ++x) {
            int imin = (static_cast<int>(x) > nGrow) ? x - nGrow : -1;
            int val = 0;
            for (int i = x; i > imin; --i) {
                val += vBadPix[i];
            }
            vBadPix_new.push_back(val ? 1 : 0);
            if (!val) {
                ++nGoodPix;
            }
        }
        vBadPix = vBadPix_new;
    }

    // return the scale of x-axis
    *nGoodPixOut = nGoodPix;

    return std::make_pair(intercept - slope, slope * xscale);
}

template <class T>
std::pair<double, double> getZScale(image::Image<T> const& image, int const nSamples, double const contrast) {
    // extract samples
    std::vector<T> vSample;
    getSample(image, nSamples, vSample);
    int nPix = vSample.size();

    if (vSample.empty()) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, "ZScale: No pixel in image is finite");
    }

    std::sort(vSample.begin(), vSample.end());

    // max, min, median
    // N.b. you can get a median in linear time, but we need the sorted array for fitLine()
    // If we wanted to speed this up, the best option would be to quantize
    // the pixel values and build a histogram
    double const zmin = vSample.front();
    double const zmax = vSample.back();
    int const iCenter = nPix / 2;
    T median = (nPix & 1) ? vSample[iCenter] : (vSample[iCenter] + vSample[iCenter + 1]) / 2;

    // fit a line to the sorted sample
    const int maxRejectionRatio = 2;
    const int npixelsMin = 5;

    int minpix = std::max(npixelsMin, nPix / maxRejectionRatio);
    int nGrow = std::max(1, nPix / 100);

    const double nSigmaClip = 2.5;
    const int nIterations = 5;

    int nGoodPix = 0;
    std::pair<double, double> ret = fitLine(&nGoodPix, vSample, nSigmaClip, nGrow, minpix, nIterations);
#if 0  // unused, but calculated and potentially useful
    double const zstart = ret.first;
#endif
    double const zslope = ret.second;

    double z1, z2;
    if (nGoodPix < minpix) {
        z1 = zmin;
        z2 = zmax;
    } else {
        double const slope = zslope / contrast;

        z1 = std::max(zmin, median - iCenter * slope);
        z2 = std::min(zmax, median + (nPix - iCenter - 1) * slope);
    }

    return std::make_pair(z1, z2);
}
//
// Explicit instantiations
#define INSTANTIATE_GETZSCALE(T)                                                                   \
    template std::pair<double, double> getZScale(image::Image<T> const& image, int const nSamples, \
                                                 double const contrast)

template <typename ImageT>
void writeBasicFits(int fd, ImageT const& data, lsst::afw::geom::SkyWcs const* Wcs = NULL,
                    char const* title = NULL);

template <typename ImageT>
void writeBasicFits(std::string const& filename, ImageT const& data, lsst::afw::geom::SkyWcs const* Wcs = NULL,
                    const char* title = NULL);

/// @cond
class Card {
public:
    Card(const std::string &name, bool val, const char *commnt = "")
            : keyword(name), value(val), comment(commnt) {}
    Card(const std::string &name, int val, const char *commnt = "")
            : keyword(name), value(val), comment(commnt) {}
    Card(const std::string &name, double val, const char *commnt = "")
            : keyword(name), value(val), comment(commnt) {}
    Card(const std::string &name, float val, const char *commnt = "")
            : keyword(name), value(val), comment(commnt) {}
    Card(const std::string &name, const std::string &val, const char *commnt = "")
            : keyword(name), value(val), comment(commnt) {}
    Card(const std::string &name, const char *val, const char *commnt = "")
            : keyword(name), value(std::string(val)), comment(commnt) {}

    ~Card() {}

    int write(int fd, int ncard, char *record) const;

    std::string keyword;
    boost::any value;
    std::string comment;
};

/*
 * Write a Card
 */
int Card::write(int fd, int ncard, char *record) const {
    char *card = &record[80 * ncard];

    if (value.type() == typeid(std::string)) {
        const char *str = boost::any_cast<std::string>(value).c_str();
        if (keyword == "" || keyword == "COMMENT" || keyword == "END" || keyword == "HISTORY") {
            sprintf(card, "%-8.8s%-72s", keyword.c_str(), str);
        } else {
            sprintf(card, "%-8.8s= '%s' %c%-*s", keyword.c_str(), str, (comment == "" ? ' ' : '/'),
                    (int)(80 - 14 - strlen(str)), comment.c_str());
        }
    } else {
        sprintf(card, "%-8.8s= ", keyword.c_str());
        card += 10;
        if (value.type() == typeid(bool)) {
            sprintf(card, "%20s", boost::any_cast<bool>(value) ? "T" : "F");
        } else if (value.type() == typeid(int)) {
            sprintf(card, "%20d", boost::any_cast<int>(value));
        } else if (value.type() == typeid(double)) {
            sprintf(card, "%20.10f", boost::any_cast<double>(value));
        } else if (value.type() == typeid(float)) {
            sprintf(card, "%20.7f", boost::any_cast<float>(value));
        }
        card += 20;
        sprintf(card, " %c%-48s", (comment == "" ? ' ' : '/'), comment.c_str());
    }
    /*
     * Write record if full
     */
    if (++ncard == 36) {
        if (::write(fd, record, FITS_SIZE) != FITS_SIZE) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Cannot write header record");
        }
        ncard = 0;
    }

    return ncard;
}
/// @endcond

/*
 * Utilities
 *
 * Flip high-order bit so as to write unsigned short to FITS.  Grrr.
 */
void flip_high_bit(char *arr,      // array that needs bits swapped
                   const int n) {  // number of bytes in arr
    if (n % 2 != 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                          (boost::format("Attempt to bit flip odd number of bytes: %d") % n).str());
    }

    unsigned short *uarr = reinterpret_cast<unsigned short *>(arr);
    for (unsigned short *end = uarr + n / 2; uarr < end; ++uarr) {
        *uarr ^= 0x8000;
    }
}

/*
 * Byte swap ABABAB -> BABABAB in place
 */
void swap_2(char *arr,      // array to swap
            const int n) {  // number of bytes
    if (n % 2 != 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                          (boost::format("Attempt to byte swap odd number of bytes: %d") % n).str());
    }

    for (char *end = arr + n; arr < end; arr += 2) {
        char t = arr[0];
        arr[0] = arr[1];
        arr[1] = t;
    }
}

/*
 * Byte swap ABCDABCD -> DCBADCBA in place (e.g. sun <--> vax)
 */
void swap_4(char *arr,      // array to swap
            const int n) {  // number of bytes
    if (n % 4 != 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                          (boost::format("Attempt to byte swap non-multiple of 4 bytes: %d") % n).str());
    }

    for (char *end = arr + n; arr < end; arr += 4) {
        char t = arr[0];
        arr[0] = arr[3];
        arr[3] = t;
        t = arr[1];
        arr[1] = arr[2];
        arr[2] = t;
    }
}

/*
 * Byte swap ABCDEFGH -> HGFEDCBA in place (e.g. sun <--> vax)
 */
void swap_8(char *arr,      // array to swap
            const int n) {  // number of bytes
    if (n % 8 != 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                          (boost::format("Attempt to byte swap non-multiple of 8 bytes: %d") % n).str());
    }

    for (char *end = arr + n; arr < end; arr += 8) {
        char t = arr[0];
        arr[0] = arr[7];
        arr[7] = t;
        t = arr[1];
        arr[1] = arr[6];
        arr[6] = t;
        t = arr[2];
        arr[2] = arr[5];
        arr[5] = t;
        t = arr[3];
        arr[3] = arr[4];
        arr[4] = t;
    }
}

int write_fits_hdr(int fd, int bitpix, int naxis, int *naxes, std::list<Card> &cards, /* extra header cards */
                   int primary) /* is this the primary HDU? */
{
    int i;
    char record[FITS_SIZE + 1]; /* write buffer */

    int ncard = 0;
    if (primary) {
        Card card("SIMPLE", true);
        ncard = card.write(fd, ncard, record);
    } else {
        Card card("XTENSION", "IMAGE");
        ncard = card.write(fd, ncard, record);
    }

    {
        Card card("BITPIX", bitpix);
        ncard = card.write(fd, ncard, record);
    }
    {
        Card card("NAXIS", naxis);
        ncard = card.write(fd, ncard, record);
    }
    for (i = 0; i < naxis; i++) {
        char key[] = "NAXIS.";
        sprintf(key, "NAXIS%d", i + 1);
        Card card(key, naxes[i]);
        ncard = card.write(fd, ncard, record);
    }
    if (primary) {
        Card card("EXTEND", true, "There may be extensions");
        ncard = card.write(fd, ncard, record);
    }
    /*
     * Write extra header cards
     */
    for (std::list<Card>::const_iterator card = cards.begin(); card != cards.end(); card++) {
        ncard = card->write(fd, ncard, record);
    }

    {
        Card card("END", "");
        ncard = card.write(fd, ncard, record);
    }
    while (ncard != 0) {
        Card card("", "");
        ncard = card.write(fd, ncard, record);
    }

    return 0;
}

/*
 * Pad out to a FITS record boundary
 */
void pad_to_fits_record(int fd,      // output file descriptor
                        int npixel,  // number of pixels already written to HDU
                        int bitpix   // bitpix for this datatype
                        ) {
    const int bytes_per_pixel = (bitpix > 0 ? bitpix : -bitpix) / 8;
    int nbyte = npixel * bytes_per_pixel;

    if (nbyte % FITS_SIZE != 0) {
        char record[FITS_SIZE + 1]; /* write buffer */

        nbyte = FITS_SIZE - nbyte % FITS_SIZE;
        memset(record, ' ', nbyte);
        if (::write(fd, record, nbyte) != nbyte) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                              "error padding file to multiple of fits block size");
        }
    }
}

int write_fits_data(int fd, int bitpix, char *begin, char *end) {
    const int bytes_per_pixel = (bitpix > 0 ? bitpix : -bitpix) / 8;
    int swap_bytes = 0;          // the default
#if defined(LSST_LITTLE_ENDIAN)  // we'll need to byte swap FITS
    if (bytes_per_pixel > 1) {
        swap_bytes = 1;
    }
#endif

    char *buff = NULL;       // I/O buffer
    bool allocated = false;  // do I need to free it?
    if (swap_bytes || bitpix == 16) {
        buff = new char[FITS_SIZE * bytes_per_pixel];
        allocated = true;
    }

    int nbyte = end - begin;
    int nwrite = (nbyte > FITS_SIZE) ? FITS_SIZE : nbyte;
    for (char *ptr = begin; ptr != end; nbyte -= nwrite, ptr += nwrite) {
        if (end - ptr < nwrite) {
            nwrite = end - ptr;
        }

        if (swap_bytes) {
            memcpy(buff, ptr, nwrite);
            if (bitpix == 16) {  // flip high-order bit
                flip_high_bit(buff, nwrite);
            }

            if (bytes_per_pixel == 2) {
                swap_2(buff, nwrite);
            } else if (bytes_per_pixel == 4) {
                swap_4(buff, nwrite);
            } else if (bytes_per_pixel == 8) {
                swap_8(buff, nwrite);
            } else {
                fprintf(stderr, "You cannot get here\n");
                abort();
            }
        } else {
            if (bitpix == 16) {  // flip high-order bit
                memcpy(buff, ptr, nwrite);
                flip_high_bit(buff, nwrite);
            } else {
                buff = ptr;
            }
        }

        if (::write(fd, buff, nwrite) != nwrite) {
            perror("Error writing image: ");
            break;
        }
    }

    if (allocated) {
        delete buff;
    }


    return (nbyte == 0 ? 0 : -1);
}

void
addWcs(std::string const & wcsName, std::list<Card> & cards, int x0=0.0, int y0=0.0)
{
    cards.push_back(Card(str(boost::format("CRVAL1%s") % wcsName),
                         x0, "(output) Column pixel of Reference Pixel"));
    cards.push_back(Card(str(boost::format("CRVAL2%s") % wcsName),
                         y0, "(output) Row pixel of Reference Pixel"));
    cards.push_back(Card(str(boost::format("CRPIX1%s") % wcsName), 1.0,
                         "Column Pixel Coordinate of Reference"));
    cards.push_back(Card(str(boost::format("CRPIX2%s") % wcsName), 1.0,
                         "Row Pixel Coordinate of Reference"));
    cards.push_back(Card(str(boost::format("CTYPE1%s") % wcsName), "LINEAR", "Type of projection"));
    cards.push_back(Card(str(boost::format("CTYPE1%s") % wcsName), "LINEAR", "Type of projection"));
    cards.push_back(Card(str(boost::format("CUNIT1%s") % wcsName), "PIXEL", "Column unit"));
    cards.push_back(Card(str(boost::format("CUNIT2%s") % wcsName), "PIXEL", "Row unit"));
}

template <typename ImageT>
void writeBasicFits(int fd,                 // file descriptor to write to
                    ImageT const &data,     // The data to write
                    geom::SkyWcs const *Wcs,  // which Wcs to use for pixel
                    char const *title       // title to write to DS9
                    ) {
    /*
     * Allocate cards for FITS headers
     */
    std::list<Card> cards;
    /*
     * What sort if image is it?
     */
    int bitpix = lsst::afw::fits::getBitPix<typename ImageT::Pixel>();
    if (bitpix == 20) {  // cfitsio for "Unsigned short"
        cards.push_back(Card("BZERO", 32768.0, ""));
        cards.push_back(Card("BSCALE", 1.0, ""));
        bitpix = 16;
    } else if (bitpix == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Unsupported image type");
    }
    /*
     * Generate WcsA, pixel coordinates, allowing for X0 and Y0
     */
    addWcs("A", cards, data.getX0(), data.getY0());
    /*
     * Now WcsB, so that pixel (0,0) is correctly labelled (but ignoring XY0)
     */
    addWcs("B", cards);

    if (title) {
        cards.push_back(Card("OBJECT", title, "Image being displayed"));
    }
    /*
     * Was there something else?
     */
    if (Wcs == NULL) {
        addWcs("", cards);              // works around a ds9 bug that WCSA/B is ignored if no Wcs is present
    } else {
        typedef std::vector<std::string> NameList;

        auto shift = lsst::geom::Extent2D(-data.getX0(), -data.getY0());
        auto newWcs = Wcs->copyAtShiftedPixelOrigin(shift);

        std::shared_ptr<lsst::daf::base::PropertySet> metadata = newWcs->getFitsMetadata();

        NameList paramNames = metadata->paramNames();

        for (NameList::const_iterator i = paramNames.begin(), end = paramNames.end(); i != end; ++i) {
            if (*i == "SIMPLE" || *i == "BITPIX" || *i == "NAXIS" || *i == "NAXIS1" || *i == "NAXIS2" ||
                *i == "XTENSION" || *i == "PCOUNT" || *i == "GCOUNT") {
                continue;
            }
            std::type_info const &type = metadata->typeOf(*i);
            if (type == typeid(bool)) {
                cards.push_back(Card(*i, metadata->get<bool>(*i)));
            } else if (type == typeid(int)) {
                cards.push_back(Card(*i, metadata->get<int>(*i)));
            } else if (type == typeid(float)) {
                cards.push_back(Card(*i, metadata->get<float>(*i)));
            } else if (type == typeid(double)) {
                cards.push_back(Card(*i, metadata->get<double>(*i)));
            } else {
                cards.push_back(Card(*i, metadata->get<std::string>(*i)));
            }
        }
    }
    /*
     * Basic FITS stuff
     */
    const int naxis = 2;  // == NAXIS
    int naxes[naxis];     /* values of NAXIS1 etc */
    naxes[0] = data.getWidth();
    naxes[1] = data.getHeight();

    write_fits_hdr(fd, bitpix, naxis, naxes, cards, 1);
    for (int y = 0; y != data.getHeight(); ++y) {
        if (write_fits_data(fd, bitpix, (char *)(data.row_begin(y)), (char *)(data.row_end(y))) < 0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                              (boost::format("Error writing data for row %d") % y).str());
        }
    }

    pad_to_fits_record(fd, data.getWidth() * data.getHeight(), bitpix);
}

template <typename ImageT>
void writeBasicFits(std::string const &filename,  // file to write, or "| cmd"
                    ImageT const &data,           // The data to write
                    geom::SkyWcs const *Wcs,        // which Wcs to use for pixel
                    char const *title             // title to write to DS9
                    ) {
    int fd;
    if ((filename.c_str())[0] == '|') {  // a command
        const char *cmd = filename.c_str() + 1;
        while (isspace(*cmd)) {
            cmd++;
        }

        fd = fileno(popen(cmd, "w"));
    } else {
        fd = ::creat(filename.c_str(), 777);
    }

    if (fd < 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                          (boost::format("Cannot open \"%s\"") % filename).str());
    }

    try {
        writeBasicFits(fd, data, Wcs, title);
    } catch (lsst::pex::exceptions::Exception &) {
        (void)::close(fd);
        throw;
    }

    (void)::close(fd);
}


/// @cond
#define INSTANTIATE(IMAGET)                                                              \
    template void writeBasicFits(int, IMAGET const &, geom::SkyWcs const *, char const *); \
    template void writeBasicFits(std::string const &, IMAGET const &, geom::SkyWcs const *, char const *)

#define INSTANTIATE_IMAGE(T) INSTANTIATE(lsst::afw::image::Image<T>)
#define INSTANTIATE_MASK(T) INSTANTIATE(lsst::afw::image::Mask<T>)

INSTANTIATE_GETZSCALE(std::uint16_t);
INSTANTIATE_GETZSCALE(float);


INSTANTIATE_IMAGE(std::uint16_t);
INSTANTIATE_IMAGE(int);
INSTANTIATE_IMAGE(float);
INSTANTIATE_IMAGE(double);
INSTANTIATE_IMAGE(std::uint64_t);

INSTANTIATE_MASK(std::uint16_t);
INSTANTIATE_MASK(image::MaskPixel);
/// @endcond

template <typename ImageT>
void declareAll(py::module &mod) {
    mod.def("writeFitsImage",
            (void (*)(int, ImageT const &, geom::SkyWcs const *, char const *)) & writeBasicFits<ImageT>,
            "fd"_a, "data"_a, "wcs"_a = NULL, "title"_a = NULL);

    mod.def("writeFitsImage",
            (void (*)(std::string const &, ImageT const &, geom::SkyWcs const *, char const *)) &
                    writeBasicFits<ImageT>,
            "filename"_a, "data"_a, "wcs"_a = NULL, "title"_a = NULL);
}

} // namespace

WRAP(Display) {
    auto dispmod = mod.def_submodule("display");
    declareAll<image::Image<std::uint16_t>>(dispmod);
    declareAll<image::Image<std::uint64_t>>(dispmod);
    declareAll<image::Image<int>>(dispmod);
    declareAll<image::Image<float>>(dispmod);
    declareAll<image::Image<double>>(dispmod);
    declareAll<image::Mask<std::uint16_t>>(dispmod);
    declareAll<image::Mask<image::MaskPixel>>(dispmod);
}

}
}
}
