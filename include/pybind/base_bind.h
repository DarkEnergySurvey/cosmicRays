#pragma once
#include "pybind/bind.h"

namespace py = pybind11;

namespace lsst {
namespace base {
WRAP(Base);
}
}
