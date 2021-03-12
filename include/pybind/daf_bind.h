#pragma once
#include "lsst/daf/base.h"

#include "pybind/bind.h"

namespace py = pybind11;
using namespace pybind11::literals;


namespace lsst {
namespace daf {
WRAP(Daf);
}
}
