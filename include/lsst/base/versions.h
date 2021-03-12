#ifndef LSST_BASE_VERSIONS_H
#define LSST_BASE_VERSIONS_H

#include <map>
#include <string>

namespace lsst {
namespace base {

/// Return version strings for dependencies
///
/// It is not clever, and only returns versions of packages declared in
/// an internal list.
///
/// Returns a map of product:version.
std::map<std::string, std::string> getRuntimeVersions();

std::string getCfitsioVersion();
std::string getFftwVersion();
std::string getWcslibVersion();
std::string getGslVersion();

}} // namespace lsst::base


#endif // include guard
