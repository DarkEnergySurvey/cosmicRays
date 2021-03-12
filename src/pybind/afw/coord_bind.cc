#include "pybind/afw_bind.h"

#include "lsst/utils/python.h"

#include "lsst/afw/coord/Weather.h"
#include "lsst/afw/coord/Observatory.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace afw {
namespace coord {
namespace {
WRAP(Weather) {
    py::class_<lsst::afw::coord::Weather> cls(mod, "Weather");

    /* Constructors */
    cls.def(py::init<double, double, double>(), "airTemperature"_a, "airPressure"_a, "humidity"_a);
    cls.def(py::init<Weather const &>(), "weather"_a);

    /* Operators */
    cls.def("__eq__", [](Weather const &self, Weather const &other) { return self == other; },
            py::is_operator());
    cls.def("__ne__", [](Weather const &self, Weather const &other) { return self != other; },
            py::is_operator());

    /* Members */
    cls.def("getAirPressure", &lsst::afw::coord::Weather::getAirPressure);
    cls.def("getAirTemperature", &lsst::afw::coord::Weather::getAirTemperature);
    cls.def("getHumidity", &lsst::afw::coord::Weather::getHumidity);
    utils::python::addOutputOp(cls, "__str__");
    utils::python::addOutputOp(cls, "__repr__");
}

WRAP(Observatory) {
    py::class_<Observatory, std::shared_ptr<Observatory>> cls(mod, "Observatory");

    /* Constructors */
    cls.def(py::init<lsst::geom::Angle const, lsst::geom::Angle const, double const>());

    /* Operators */
    cls.def("__eq__", [](Observatory const& self, Observatory const& other) { return self == other; },
            py::is_operator());
    cls.def("__ne__", [](Observatory const& self, Observatory const& other) { return self != other; },
            py::is_operator());
    cls.def("__str__", &Observatory::toString);
    cls.def("__repr__", &Observatory::toString);

    /* Members */
    cls.def("getLongitude", &Observatory::getLongitude);
    cls.def("getLatitude", &Observatory::getLatitude);
    cls.def("getElevation", &Observatory::getElevation);
    cls.def("setLongitude", &Observatory::setLongitude, "longitude"_a);
    cls.def("setLatitude", &Observatory::setLatitude, "latitude"_a);
    cls.def("setElevation", &Observatory::setElevation, "elevation"_a);

}
}

WRAP(Coord) {
    auto coordmod = mod.def_submodule("coord");
    wrapWeather(coordmod);
    wrapObservatory(coordmod);
}
}
}
}
