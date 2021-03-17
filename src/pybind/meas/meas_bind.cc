#include "pybind/meas_bind.h"

namespace lsst {
namespace meas {
WRAP(Meas) {\
    py::print("CALL A");

    auto bmod = mod.def_submodule("base");
    base::wrapBase(bmod);
    py::print("CALL B");
    auto amod = mod.def_submodule("algorithms");
    algorithms::wrapAlg(amod);
}
}
}
