#include "pybind/afw_bind.h"

#include "lsst/afw/formatters/Utils.h"

namespace py = pybind11;
using namespace py::literals;
namespace lsst {
namespace afw {
namespace formatters {

WRAP(Formatters) {
    mod.def("stringToBytes", stringToBytes);
    mod.def("bytesToString", bytesToString);
}

}
}
}
