#include <cstddef>
#include <cstdlib>

#include <sstream>

#include "lsst/base/library.h"
#include "lsst/base/versions.h"

namespace lsst {
namespace base {

namespace {

typedef std::string (*VersionGetter)(void);

// List of packages and how to determine the version
std::map<std::string, VersionGetter> const packages {
    {"cfitsio", getCfitsioVersion},
    {"fftw", getFftwVersion},
    {"wcslib", getWcslibVersion},
    {"gsl", getGslVersion},
};

} // anonymous namespace

std::string getCfitsioVersion()
{
    typedef float (GetVersion)(float*);
    float version;
    loadSymbol<GetVersion>("libcfitsio", "ffvers")(&version);

    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    ss << version;
    return ss.str();
}

std::string getFftwVersion()
{
    return std::string(loadSymbol<char const>("libfftw3", "fftw_version"));
}

std::string getWcslibVersion()
{
    typedef char const* (GetVersion)(int[]);
    return std::string(loadSymbol<GetVersion>("libwcs", "wcslib_version")(NULL));
}

std::string getGslVersion()
{
    return std::string(*loadSymbol<char const*>("libgsl", "gsl_version"));
}


std::map<std::string, std::string> getRuntimeVersions()
{
    std::map<std::string, std::string> versions;
    for (auto&& pkg : packages) {
        try {
            versions[pkg.first] = pkg.second();
        } catch (LibraryException const&) {
            // Can't find the module, so ignore it
        }
    }
    return versions;
}


}} // namespace lsst::base
