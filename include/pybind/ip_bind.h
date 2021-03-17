#pragma once
#include "pybind/bind.h"
namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace isr {
WRAP(IP);
}
}
}
