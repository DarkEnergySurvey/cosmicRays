#include "pybind/utils_bind.h"
#include "lsst/utils/Backtrace.h"
#include "lsst/utils/packaging.h"
#include "lsst/utils/Demangle.h"

namespace lsst {
namespace utils {

namespace {
WRAP(Backtrace) {
    Backtrace &backtrace = Backtrace::get();
    // Trick to tell the compiler backtrace is used and should not be
    // optimized away, as well as convenient way to check if backtrace
    // is enabled.
    mod.def("isEnabled", [&backtrace]() -> bool { return backtrace.isEnabled(); });
}

WRAP(Packaging) {
    mod.def("getPackageDir", getPackageDir);
}

WRAP(Demangle) {
    mod.def("demangleType", demangleType);
}

}

WRAP(Utils) {
    auto bmod = mod.def_submodule("backtrace");
    wrapBacktrace(bmod);

    wrapPackaging(mod);
    wrapDemangle(mod);
}

}
}
