#pragma once
#include "pybind/bind.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {
WRAP(Base);
}

namespace algorithms {
WRAP(Alg);
}

WRAP(Meas);
}
}
