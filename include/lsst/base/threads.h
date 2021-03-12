#ifndef LSST_BASE_THREADS_H
#define LSST_BASE_THREADS_H

#include <stdexcept>
#include <string>

namespace lsst {
namespace base {

extern bool const haveOpenBlas;         ///< Is OpenBLAS available?
extern bool const haveMkl;              ///< Is MKL available?

/// Environment variable to allow implicit threading
///
/// Used by disableImplicitThreading.
std::string const allowEnvvar = "LSST_ALLOW_IMPLICIT_THREADS";


/// No threading library is available
class NoThreadsException : public std::runtime_error {
public:
    NoThreadsException() : std::runtime_error("No threading library is available") {};
};

/// Are threaded packages available?
bool haveThreads() { return haveOpenBlas || haveMkl; }

/// Set number of threads to use
///
/// @throw NoThreadsException if no threading library is available
void setNumThreads(unsigned int numThreads);

/// Get maximum number of threads we might use
///
/// Returns the maximum value of the number of threads being used by
/// the threading libraries that are available.
unsigned int getNumThreads();

/// Disable threading that has not been set explicitly
///
/// Some threaded packages implicitly use multiple threads if the user doesn't
/// explicitly state the number of desired threads. However, this can interfere
/// with operations that are parallelised at a higher level. This function will
/// disable threading unless the user has explicitly specified the number of
/// desired threads through environment variables.
///
/// This behavior may be disabled by setting the environment variable specified
/// by allowEnvvar.
///
/// This is principally intended for Linux machines (we explicitly load .so
/// dynamic libraries); MacOS has its own way of doing threading (Grand Central
/// Dispatch) that throttles threads to avoid overwhelming the machine.
///
/// @ return whether we disabled threading
bool disableImplicitThreading();

}} // namespace lsst::base


#endif // include guard
