#include "pybind/meas_bind.h"

namespace lsst {
namespace meas {
WRAP(Meas) {\
    auto bmod = mod.def_submodule("base");
    base::wrapBase(bmod);

    auto amod = mod.def_submodule("algorithms");
    algorithms::wrapAlg(amod);
}
}
}
