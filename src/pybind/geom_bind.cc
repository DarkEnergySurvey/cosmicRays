#include "pybind/geom_bind.h"
#include "lsst/utils/python.h"
#include "lsst/geom.h"
#include "pybind11/numpy.h"
namespace python = lsst::utils::python;

namespace lsst {
namespace geom {

namespace {
template <typename Derived, typename T, int N>
using PyCoordinateBase = py::class_<geom::CoordinateBase<Derived, T, N>>;

template <int N>
using PyCoordinateExpr = py::class_<geom::CoordinateExpr<N>, geom::CoordinateBase<geom::CoordinateExpr<N>, bool, N>>;

template <typename T, int N>
using PyExtentBase = py::class_<geom::ExtentBase<T, N>, geom::CoordinateBase<geom::Extent<T, N>, T, N>>;

template <typename T, int N>
using PyExtent = py::class_<geom::Extent<T, N>, geom::ExtentBase<T, N>>;

template <typename T, int N>
using PyPointBase = py::class_<geom::PointBase<T, N>, geom::CoordinateBase<geom::Point<T, N>, T, N>>;

template <typename T, int N>
using PyPoint = py::class_<geom::Point<T, N>, geom::PointBase<T, N>>;

template <typename PyClass>
void declareCommonIntervalInterface(PyClass &cls) {
    using T = typename PyClass::type;
    using Element = typename T::Element;
    cls.def(py::init<>());
    // It's not clear why py::overload cast doesn't work for the next two
    // declarations - maybe the templated context matters?
    cls.def_static("fromSpannedPoints", [](ndarray::Array<Element const, 1> const &elements) {
        return T::fromSpannedPoints(elements);
    });
    cls.def_static("fromSpannedPoints",
                   [](std::vector<Element> const &elements) { return T::fromSpannedPoints(elements); });
    cls.def(py::init([](py::kwargs kw) -> T {
        if (kw.size() == 2) {
            if (kw.contains("min")) {
                if (kw.contains("max")) {
                    return T::fromMinMax(py::cast<Element>(kw["min"]), py::cast<Element>(kw["max"]));
                }
                if (kw.contains("size")) {
                    return T::fromMinSize(py::cast<Element>(kw["min"]), py::cast<Element>(kw["size"]));
                }
            }
            if (kw.contains("max") && kw.contains("size")) {
                return T::fromMaxSize(py::cast<Element>(kw["max"]), py::cast<Element>(kw["size"]));
            }
            if (kw.contains("center") && kw.contains("size")) {
                return T::fromCenterSize(py::cast<Element>(kw["center"]), py::cast<Element>(kw["size"]));
            }
        }
        PyErr_SetString(PyExc_TypeError,
                        "General constructor requires exactly 2 of the following keyword-only "
                        "arguments: (min, max, center, size).");
        throw py::error_already_set();
    }));
    cls.def(py::init<T const &>());
    cls.def("__eq__", [](T const &self, T const &other) { return self == other; }, py::is_operator());
    cls.def("__ne__", [](T const &self, T const &other) { return self != other; }, py::is_operator());
    cls.def("getMin", &T::getMin);
    cls.def_property_readonly("min", &T::getMin);
    cls.def("getMax", &T::getMax);
    cls.def_property_readonly("max", &T::getMax);
    cls.def("getSize", &T::getSize);
    cls.def_property_readonly("size", &T::getSize);
    cls.def("isEmpty", &T::isEmpty);
    cls.def("contains", py::overload_cast<T const &>(&T::contains, py::const_));
    cls.def("contains", py::vectorize(static_cast<bool (T::*)(Element) const noexcept>(&T::contains)));
    cls.def("__contains__", py::overload_cast<Element>(&T::contains, py::const_));
    cls.def("__contains__", py::overload_cast<T const &>(&T::contains, py::const_));
    cls.def("overlaps", &T::overlaps);
    cls.def("intersects", &T::intersects);
    cls.def("isDisjointFrom", &T::isDisjointFrom);
    cls.def("dilatedBy", &T::dilatedBy);
    cls.def("erodedBy", &T::erodedBy);
    cls.def("shiftedBy", &T::shiftedBy);
    cls.def("reflectedAbout", &T::reflectedAbout);
    cls.def("expandedTo", py::overload_cast<Element>(&T::expandedTo, py::const_));
    cls.def("expandedTo", py::overload_cast<T const &>(&T::expandedTo, py::const_));
    cls.def("clippedTo", &T::clippedTo);
    cls.def("__str__", &T::toString);
    utils::python::addOutputOp(cls, "__repr__");
    cls.def("__reduce__", [cls](geom::IntervalD const &self) {
        return py::make_tuple(cls, make_tuple(py::cast(self.getMin()), py::cast(self.getMax())));
    });
}


template <typename Derived, typename T, int N>
void declareCoordinateBase(py::module_ &mod, std::string const &suffix) {
    static std::string const name = "CoordinateBase" + suffix;
    PyCoordinateBase<Derived, T, N>(mod, name.c_str())
            .def("__getitem__", [](geom::CoordinateBase<Derived, T, N> &self, int i) -> T {
        return self[utils::python::cppIndex(N, i)];
    })
    .def("__setitem__", [](geom::CoordinateBase<Derived, T, N> &self, int i, T value) {
        self[utils::python::cppIndex(N, i)] = value;
    })
    .def("__len__", [](geom::CoordinateBase<Derived, T, N> &c) -> int { return N; });
}

// Declare mixed-type and type-overloaded operators for Extent with dimension
// N for both int and double. Because pybind11 tries operators (like any
// overload) `in order', int has to come before double in any overloaded
// operators that dispatch on a scalar, and hence they have to be defined here
// instead of declareExtent.
template <int N>
void declareExtentOperators(py::module &mod, PyExtent<int, N> &clsI, PyExtent<double, N> &clsD) {
    // Python's integer division works differently than C++'s for negative numbers - Python
    // uses floor (rounds towards more negative), while C++ truncates (rounds towards zero).
    // Therefore one needs to be careful in the definition of division operators.
    clsI.def("__floordiv__",
             [](geom::Extent<int, N> const &self, int other) -> geom::Extent<int, N> {
        return floor(self / static_cast<double>(other));
    },
    py::is_operator());

    clsI.def("__truediv__", [](geom::Extent<int, N> const &self, double other) { return self / other; },
    py::is_operator());
    clsD.def("__truediv__", [](geom::Extent<double, N> const &self, double other) { return self / other; },
    py::is_operator());

    clsI.def("__ifloordiv__", [](geom::Extent<int, N> &self, int other) -> geom::Extent<int, N> & {
        self = floor(self / static_cast<double>(other));
        return self;
    });

    clsI.def("__itruediv__", [](geom::Extent<int, N> &self, double other) {
        PyErr_SetString(PyExc_TypeError, "In-place true division not supported for Extent<int,N>.");
        throw py::error_already_set();
    });
    clsD.def("__itruediv__", [](geom::Extent<double, N> &self, double other) -> geom::Extent<double, N> & {
        self /= other;
        return self;
    });

    clsI.def("__iadd__", [](geom::Extent<int, N> &self, geom::Extent<int, N> const &other) -> geom::Extent<int, N> & {
        self += other;
        return self;
    });
    clsD.def("__iadd__",
             [](geom::Extent<double, N> &self, geom::Extent<double, N> const &other) -> geom::Extent<double, N> & {
        self += other;
        return self;
    }
    );
    clsD.def("__iadd__",
             [](geom::Extent<double, N> &self, geom::Extent<int, N> const &other) -> geom::Extent<double, N> & {
        self += other;
        return self;
    }
    );

    clsI.def("__isub__", [](geom::Extent<int, N> &self, geom::Extent<int, N> const &other) -> geom::Extent<int, N> & {
        self -= other;
        return self;
    });
    clsD.def("__isub__", [](geom::Extent<double, N> &self, geom::Extent<double, N> const &other) -> geom::Extent<double, N> & {
        self -= other;
        return self;
    });
    clsD.def("__isub__", [](geom::Extent<double, N> &self, geom::Extent<int, N> const &other) -> geom::Extent<double, N> & {
        self -= other;
        return self;
    });

    clsI.def("__imul__", [](geom::Extent<int, N> &self, int other) -> geom::Extent<int, N> & {
        self *= other;
        return self;
    });
    clsD.def("__imul__", [](geom::Extent<double, N> &self, int other) -> geom::Extent<double, N> & {
        self *= other;
        return self;
    });
    clsD.def("__imul__", [](geom::Extent<double, N> &self, double other) -> geom::Extent<double, N> & {
        self *= other;
        return self;
    });

    // Operator-like free functions
    mod.def("truncate", geom::truncate<N>);
    mod.def("floor", geom::floor<N>);
    mod.def("ceil", geom::ceil<N>);
    // And method versions, since Python doens't have ADL
    clsD.def("truncate", geom::truncate<N>);
    clsD.def("floor", geom::floor<N>);
    clsD.def("ceil", geom::ceil<N>);
}

template <typename T, int N>
void declareExtentBase(py::module &mod, std::string const &suffix) {
    static std::string const name = "ExtentBase" + suffix;
    declareCoordinateBase<geom::Extent<T, N>, T, N>(mod, "Extent" + suffix);
    PyExtentBase<T, N>(mod, name.c_str())
            // These are not the usual Python double-underscore operators - they do elementwise comparisons,
            // returning a CoordinateExpr object with boolean x and y values.  NumPy examples to the contrary
            // notwithstanding, true Python comparison operators are expected to return scalar bools.
            .def("eq", [](geom::ExtentBase<T, N> const &self, geom::Extent<T, N> other) { return self.eq(other); })
            .def("ne", [](geom::ExtentBase<T, N> const &self, geom::Extent<T, N> other) { return self.ne(other); })
            .def("lt", [](geom::ExtentBase<T, N> const &self, geom::Extent<T, N> other) { return self.lt(other); })
            .def("le", [](geom::ExtentBase<T, N> const &self, geom::Extent<T, N> other) { return self.le(other); })
            .def("gt", [](geom::ExtentBase<T, N> const &self, geom::Extent<T, N> other) { return self.gt(other); })
            .def("ge", [](geom::ExtentBase<T, N> const &self, geom::Extent<T, N> other) { return self.ge(other); })
            .def("eq", [](geom::ExtentBase<T, N> const &self, T other) { return self.eq(other); })
            .def("ne", [](geom::ExtentBase<T, N> const &self, T other) { return self.ne(other); })
            .def("lt", [](geom::ExtentBase<T, N> const &self, T other) { return self.lt(other); })
            .def("le", [](geom::ExtentBase<T, N> const &self, T other) { return self.le(other); })
            .def("gt", [](geom::ExtentBase<T, N> const &self, T other) { return self.gt(other); })
            .def("ge", [](geom::ExtentBase<T, N> const &self, T other) { return self.ge(other); })

            .def("asPoint", &geom::ExtentBase<T, N>::asPoint)
            .def("computeNorm", &geom::ExtentBase<T, N>::computeNorm)
            .def("computeSquaredNorm", &geom::ExtentBase<T, N>::computeSquaredNorm);

}

// Common functionality for all Extents, including declaring base classes.
template <typename T, int N>
PyExtent<T, N> declareExtent(py::module &mod, std::string const &suffix) {
    static std::string const name = "Extent" + suffix;
    declareExtentBase<T, N>(mod, suffix);
    return PyExtent<T, N>(mod, name.c_str())
            .def(py::init<T>(), "value"_a = static_cast<T>(0))
            .def(py::init<geom::Point<int, N> const &>())
            .def(py::init<geom::Point<T, N> const &>())
            .def(py::init<geom::Extent<int, N> const &>())
            .def(py::init<geom::Extent<T, N> const &>())
            .def(py::init<typename geom::Extent<T, N>::EigenVector>())
            .def("__neg__", [](geom::Extent<T, N> const &self) { return -self; })
            .def("__pos__", [](geom::Extent<T, N> const &self) { return self; })
            .def("__mul__", [](geom::Extent<T, N> const &self, int other) { return self * other; },
                    py::is_operator())
            .def("__mul__", [](geom::Extent<T, N> const &self, double other) { return self * other; },
                  py::is_operator())
            .def("__rmul__", [](geom::Extent<T, N> const &self, int other) { return self * other; },
                    py::is_operator())
            .def("__rmul__", [](geom::Extent<T, N> const &self, double other) { return self * other; },
                    py::is_operator())
            .def("__add__",
                    [](geom::Extent<T, N> const &self, geom::Extent<int, N> const &other) { return self + other; },
                    py::is_operator())
            .def("__add__",
                    [](geom::Extent<T, N> const &self, geom::Extent<double, N> const &other) { return self + other; },
                    py::is_operator())
            .def("__add__",
                    [](geom::Extent<T, N> const &self, geom::Point<int, N> const &other) {
                        return self + geom::Point<T, N>(other);
                    },
                    py::is_operator())
            .def("__add__",
                    [](geom::Extent<T, N> const &self, geom::Point<double, N> const &other) { return self + other; },
                    py::is_operator())
            .def("__sub__",
                    [](geom::Extent<T, N> const &self, geom::Extent<int, N> const &other) {
                        return self - geom::Extent<T, N>(other);
                    },
                    py::is_operator())
            .def("__sub__",
                    [](geom::Extent<T, N> const &self, geom::Extent<double, N> const &other) { return self - other; },
                    py::is_operator())
            .def("__eq__",
                    [](geom::Extent<T, N> const &self, geom::Extent<T, N> const &other) { return self == other; },
                    py::is_operator())
            .def("__ne__",
                    [](geom::Extent<T, N> const &self, geom::Extent<T, N> const &other) { return self != other; },
                    py::is_operator())
            .def("clone", [](geom::Extent<T, N> const &self) { return geom::Extent<T, N>{self}; });
}

// Add functionality only found in N=2 Extents (and delgate to declareExtent for the rest)
template <typename T>
PyExtent<T, 2> declareExtent2(py::module &mod, std::string const &suffix) {
    auto cls = declareExtent<T, 2>(mod, std::string("2") + suffix);
    /* Members types and enums */
    cls.def_property_readonly_static("dimensions", [](py::object /* cls */) { return 2; });
    /* Constructors */
    cls.def(py::init<int, int>(), "x"_a, "y"_a);
    cls.def(py::init<double, double>(), "x"_a, "y"_a);
    /* Members */
    auto getX = py::overload_cast<>(&geom::Extent<T, 2>::getX, py::const_);
    auto getY = py::overload_cast<>(&geom::Extent<T, 2>::getY, py::const_);
    cls.def("getX", getX);
    cls.def("getY", getY);
    cls.def("setX", &geom::Extent<T, 2>::setX);
    cls.def("setY", &geom::Extent<T, 2>::setY);
    cls.def_property("x", getX, &geom::Extent<T, 2>::setX);
    cls.def_property("y", getY, &geom::Extent<T, 2>::setY);
    return cls;
}

// Add functionality only found in N=3 Extents (and delgate to declareExtent for the rest)
template <typename T>
PyExtent<T, 3> declareExtent3(py::module &mod, const std::string &suffix) {
    auto cls = declareExtent<T, 3>(mod, std::string("3") + suffix);
    /* Member types and enums */
    cls.def_property_readonly_static("dimensions", [](py::object /* cls */) { return 3; });
    /* Constructors */
    cls.def(py::init<int, int, int>(), "x"_a, "y"_a, "z"_a);
    cls.def(py::init<double, double, double>(), "x"_a, "y"_a, "z"_a);
    /* Members */
    auto getX = py::overload_cast<>(&geom::Extent<T, 3>::getX, py::const_);
    auto getY = py::overload_cast<>(&geom::Extent<T, 3>::getY, py::const_);
    auto getZ = py::overload_cast<>(&geom::Extent<T, 3>::getZ, py::const_);
    cls.def("getX", getX);
    cls.def("getY", getY);
    cls.def("getZ", getZ);
    cls.def("setX", &geom::Extent<T, 3>::setX);
    cls.def("setY", &geom::Extent<T, 3>::setY);
    cls.def("setZ", &geom::Extent<T, 3>::setZ);
    cls.def_property("x", getX, &geom::Extent<T, 3>::setX);
    cls.def_property("y", getY, &geom::Extent<T, 3>::setY);
    cls.def_property("z", getZ, &geom::Extent<T, 3>::setZ);
    return cls;
}




template <int N>
void declareCoordinateExpr(py::module &mod, std::string const &suffix) {
    static std::string const name = "CoordinateExpr" + suffix;
    declareCoordinateBase<geom::CoordinateExpr<N>, bool, N>(mod, name);
    PyCoordinateExpr<N>(mod, name.c_str())
            .def(py::init<bool>(), "val"_a = false)
            .def("and_", &geom::CoordinateExpr<N>::and_)
            .def("or_", &geom::CoordinateExpr<N>::or_)
            .def("not_", &geom::CoordinateExpr<N>::not_);
    mod.def("all", geom::all<N>);
    mod.def("any", geom::any<N>);
}

template <typename OtherT>
void declareAngleComparisonOperators(py::class_<geom::Angle>& cls) {
    cls.def("__eq__", [](geom::Angle const& self, OtherT const& other) { return self == other; },
            py::is_operator());
    cls.def("__ne__", [](geom::Angle const& self, OtherT const& other) { return self != other; },
            py::is_operator());
    cls.def("__le__", [](geom::Angle const& self, OtherT const& other) { return self <= other; },
            py::is_operator());
    cls.def("__ge__", [](geom::Angle const& self, OtherT const& other) { return self >= other; },
            py::is_operator());
    cls.def("__lt__", [](geom::Angle const& self, OtherT const& other) { return self < other; }, py::is_operator());
    cls.def("__gt__", [](geom::Angle const& self, OtherT const& other) { return self > other; }, py::is_operator());
}



void wrapAngle(py::module &mod) {
    py::class_<AngleUnit>(mod, "AngleUnit3")
            .def("__eq__", [](AngleUnit const& self, AngleUnit const& other) { return self == other; },
    py::is_operator())
            .def("__ne__", [](AngleUnit const& self, AngleUnit const& other) { return !(self == other); },
    py::is_operator())
            .def("_mul", [](AngleUnit const& self, double other) { return other * self; },
    py::is_operator())
            .def("_rmul", [](AngleUnit const& self, double other) { return other * self; },
    py::is_operator());
    mod.attr("radians") = py::cast(radians);
    mod.attr("degrees") = py::cast(degrees);
    mod.attr("hours") = py::cast(hours);
    mod.attr("arcminutes") = py::cast(arcminutes);
    mod.attr("arcseconds") = py::cast(arcseconds);
    mod.attr("milliarcseconds") = py::cast(milliarcseconds);

    py::class_<Angle> cls(mod, "Angle");
    cls.def(py::init<double, AngleUnit>(), py::arg("val"), py::arg("units") = radians);
    cls.def(py::init<>());
    declareAngleComparisonOperators<Angle>(cls);
    declareAngleComparisonOperators<double>(cls);
    declareAngleComparisonOperators<int>(cls);
    cls.def("__mul__", [](Angle const& self, double other) { return self * other; },
    py::is_operator());
    cls.def("__mul__", [](Angle const& self, int other) { return self * other; },
    py::is_operator());
    cls.def("__rmul__", [](Angle const& self, double other) { return self * other; },
    py::is_operator());
    cls.def("__rmul__", [](Angle const& self, int other) { return self * other; },
    py::is_operator());
    cls.def("__imul__", [](Angle& self, double other) { return self *= other; });
    cls.def("__imul__", [](Angle& self, int other) { return self *= other; });
    cls.def("__add__", [](Angle const& self, Angle const& other) { return self + other; },
    py::is_operator());
    cls.def("__sub__", [](Angle const& self, Angle const& other) { return self - other; },
    py::is_operator());
    cls.def("__neg__", [](Angle const& self) { return -self; }, py::is_operator());
    cls.def("__iadd__", [](Angle& self, Angle const& other) { return self += other; });
    cls.def("__isub__", [](Angle& self, Angle const& other) { return self -= other; });
    cls.def("__truediv__", [](Angle const& self, double other) { return self / other; },
    py::is_operator());
    // Without an explicit wrapper, Python lets Angle / Angle -> Angle
    cls.def("__truediv__", [](Angle const& self, Angle const& other) {
        throw py::type_error("unsupported operand type(s) for /: 'Angle' and 'Angle'");
    });
    cls.def("__float__", &Angle::operator double);
    cls.def("__abs__", [](Angle const& self) { return std::abs(self.asRadians()) * radians; });

    cls.def("__reduce__", [cls](Angle const& self) {
        return py::make_tuple(cls, py::make_tuple(py::cast(self.asRadians())));
    });
    utils::python::addOutputOp(cls, "__str__");
    utils::python::addOutputOp(cls, "__repr__");
    cls.def("asAngularUnits", &Angle::asAngularUnits);
    cls.def("asRadians", &Angle::asRadians);
    cls.def("asDegrees", &Angle::asDegrees);
    cls.def("asHours", &Angle::asHours);
    cls.def("asArcminutes", &Angle::asArcminutes);
    cls.def("asArcseconds", &Angle::asArcseconds);
    cls.def("asMilliarcseconds", &Angle::asMilliarcseconds);
    cls.def("wrap", &Angle::wrap);
    cls.def("wrapCtr", &Angle::wrapCtr);
    cls.def("wrapNear", &Angle::wrapNear);
    cls.def("separation", &Angle::separation);
    mod.def("isAngle", isAngle<Angle>);
    mod.def("isAngle", isAngle<double>);

    mod.attr("PI") = py::float_(PI);
    mod.attr("TWOPI") = py::float_(TWOPI);
    mod.attr("HALFPI") = py::float_(HALFPI);
    mod.attr("ONE_OVER_PI") = py::float_(ONE_OVER_PI);
    mod.attr("SQRTPI") = py::float_(SQRTPI);
    mod.attr("INVSQRTPI") = py::float_(INVSQRTPI);
    mod.attr("ROOT2") = py::float_(ROOT2);
    mod.def("degToRad", degToRad);
    mod.def("radToDeg", radToDeg);
    mod.def("radToArcsec", radToArcsec);
    mod.def("radToMas", radToMas);
    mod.def("arcsecToRad", arcsecToRad);
    mod.def("masToRad", masToRad);
}

template <typename T, int N>
void declarePointBase(py::module &mod, std::string const &suffix) {
    static std::string const name = "PointBase" + suffix;
    declareCoordinateBase<geom::Point<T, N>, T, N>(mod, "Point" + suffix);
    PyPointBase<T, N>(mod, name.c_str())
            // These are not the usual Python double-underscore operators - they do elementwise comparisons,
            // returning a CoordinateExpr object with boolean x and y values.  NumPy examples to the contrary
            // notwithstanding, true Python comparison operators are expected to return scalar bools.
            .def("eq", [](geom::PointBase<T, N> const &self, geom::Point<T, N> const &rhs) { return self.eq(rhs); })
            .def("ne", [](geom::PointBase<T, N> const &self, geom::Point<T, N> const &rhs) { return self.ne(rhs); })
            .def("lt", [](geom::PointBase<T, N> const &self, geom::Point<T, N> const &rhs) { return self.lt(rhs); })
            .def("le", [](geom::PointBase<T, N> const &self, geom::Point<T, N> const &rhs) { return self.le(rhs); })
            .def("gt", [](geom::PointBase<T, N> const &self, geom::Point<T, N> const &rhs) { return self.gt(rhs); })
            .def("ge", [](geom::PointBase<T, N> const &self, geom::Point<T, N> const &rhs) { return self.ge(rhs); })
            .def("eq", [](geom::PointBase<T, N> const &self, T rhs) { return self.eq(rhs); })
            .def("ne", [](geom::PointBase<T, N> const &self, T rhs) { return self.ne(rhs); })
            .def("lt", [](geom::PointBase<T, N> const &self, T rhs) { return self.lt(rhs); })
            .def("le", [](geom::PointBase<T, N> const &self, T rhs) { return self.le(rhs); })
            .def("gt", [](geom::PointBase<T, N> const &self, T rhs) { return self.gt(rhs); })
            .def("ge", [](geom::PointBase<T, N> const &self, T rhs) { return self.ge(rhs); })
            /* Members */
            .def("asExtent", &geom::PointBase<T, N>::asExtent)
            .def("shift", &geom::PointBase<T, N>::shift)
            .def("scale", &geom::PointBase<T, N>::scale)
            .def("distanceSquared", &geom::PointBase<T, N>::distanceSquared)
            .def("toString", &geom::PointBase<T, N>::toString);
}

// Common functionality
template <typename T, int N>
PyPoint<T, N> declarePoint(py::module &mod, std::string const &suffix) {
    static std::string const name = "Point" + suffix;
    declarePointBase<T, N>(mod, suffix);
    return PyPoint<T, N>(mod, name.c_str())
            /* Constructors */
            .def(py::init<T>(), "value"_a = static_cast<T>(0))
            // Note that we can't use T here because both types are needed
            .def(py::init<geom::Point<double, N> const &>())
            .def(py::init<geom::Point<int, N> const &>())
            .def(py::init<geom::Extent<T, N> const &>())
            .def(py::init<typename geom::Point<T, N>::EigenVector>())
            /* Operators */
            .def("__add__", [](geom::Point<T, N> const &self, geom::Extent<double, N> &other) { return self + other; },
                    py::is_operator())
            .def("__add__", [](geom::Point<T, N> const &self, geom::Extent<int, N> &other) { return self + other; },
                    py::is_operator())
            .def("__sub__", [](geom::Point<T, N> const &self, geom::Point<T, N> &other) { return self - other; },
                    py::is_operator())
            .def("__sub__", [](geom::Point<T, N> const &self, geom::Extent<T, N> &other) { return self - other; },
                    py::is_operator())
            .def("__sub__", [](geom::Point<T, N> const &self, geom::Point<double, N> &other) { return self - other; },
                    py::is_operator())
            .def("__sub__", [](geom::Point<T, N> const &self, geom::Point<int, N> &other) { return self - other; },
                    py::is_operator())
            .def("__sub__", [](geom::Point<T, N> const &self, geom::Extent<double, N> &other) { return self - other; },
                    py::is_operator())
            .def("__sub__", [](geom::Point<T, N> const &self, geom::Extent<int, N> &other) { return self - other; },
                    py::is_operator())
            .def("__eq__", [](geom::Point<T, N> const &self, geom::Point<T, N> const &other) { return self == other; },
                    py::is_operator())
            .def("__ne__", [](geom::Point<T, N> const &self, geom::Point<T, N> const &other) { return self != other; },
                    py::is_operator())
            /* Members */
            .def("clone", [](geom::Point<T, N> const &self) { return geom::Point<T, N>{self}; });

}

// Add functionality only found in N=2 Points
template <typename T>
PyPoint<T, 2> declarePoint2(py::module &mod, std::string const &suffix) {
    auto cls = declarePoint<T, 2>(mod, std::string("2") + suffix);
    /* Member types and enums */
    cls.def_property_readonly_static("dimensions", [](py::object /* cls */) { return 2; });
    /* Constructors */
    cls.def(py::init<int, int>(), "x"_a, "y"_a);
    cls.def(py::init<double, double>(), "x"_a, "y"_a);
    /* Members */
    auto getX = py::overload_cast<>(&geom::Point<T, 2>::getX, py::const_);
    auto getY = py::overload_cast<>(&geom::Point<T, 2>::getY, py::const_);
    cls.def("getX", getX);
    cls.def("getY", getY);
    cls.def("setX", &geom::Point<T, 2>::setX);
    cls.def("setY", &geom::Point<T, 2>::setY);
    cls.def_property("x", getX, &geom::Point<T, 2>::setX);
    cls.def_property("y", getY, &geom::Point<T, 2>::setY);
    return cls;
}

// Add functionality only found in N=3 Points
template <typename T>
PyPoint<T, 3> declarePoint3(py::module &mod, std::string const &suffix) {
    auto cls = declarePoint<T, 3>(mod, std::string("3") + suffix);
    /* Member types and enums */
    cls.def_property_readonly_static("dimensions", [](py::object /* cls */) { return 3; });
    /* Constructors */
    cls.def(py::init<int, int, int>(), "x"_a, "y"_a, "z"_a);
    cls.def(py::init<double, double, double>(), "x"_a, "y"_a, "z"_a);
    /* Members */
    auto getX = py::overload_cast<>(&geom::Point<T, 3>::getX, py::const_);
    auto getY = py::overload_cast<>(&geom::Point<T, 3>::getY, py::const_);
    auto getZ = py::overload_cast<>(&geom::Point<T, 3>::getZ, py::const_);
    cls.def("getX", getX);
    cls.def("getY", getY);
    cls.def("getZ", getZ);
    cls.def("setX", &geom::Point<T, 3>::setX);
    cls.def("setY", &geom::Point<T, 3>::setY);
    cls.def("setZ", &geom::Point<T, 3>::setZ);
    cls.def_property("x", getX, &geom::Point<T, 3>::setX);
    cls.def_property("y", getY, &geom::Point<T, 3>::setY);
    cls.def_property("z", getZ, &geom::Point<T, 3>::setZ);
    return cls;
}

// Declare mixed-type and type-overloaded operators for Point with dimension
// N for both int and double. Because pybind11 tries operators (like any
// overload) `in order', int has to come before double in any overloaded
// operators that dispatch on a scalar, and hence they have to be defined here
// instead of declareExtent.
template <int N>
void declarePointOperators(py::module &mod, PyPoint<int, N> &clsI, PyPoint<double, N> &clsD) {
    clsI.def("__iadd__", [](geom::Point<int, N> &self, geom::Extent<int, N> const &other) {
        self += other;
        return &self;
    });
    clsD.def("__iadd__", [](geom::Point<double, N> &self, geom::Extent<int, N> const &other) {
        self += other;
        return &self;
    });
    clsD.def("__iadd__", [](geom::Point<double, N> &self, geom::Extent<double, N> const &other) {
        self += other;
        return &self;
    });
    clsI.def("__isub__", [](geom::Point<int, N> &self, geom::Extent<int, N> const &other) {
        self -= other;
        return &self;
    });
    clsD.def("__isub__", [](geom::Point<double, N> &self, geom::Extent<int, N> const &other) {
        self -= other;
        return &self;
    });
    clsD.def("__isub__", [](geom::Point<double, N> &self, geom::Extent<double, N> const &other) {
        self -= other;
        return &self;
    });
}



using PySpherePoint = py::class_<geom::SpherePoint, std::shared_ptr<geom::SpherePoint>>;


WRAP(Coordinates) {
    // Only the interface-level classes are defined here, and these functions
    // call others to define their base classes, since those are never shared.

    declareCoordinateExpr<2>(mod, "2");
    declareCoordinateExpr<3>(mod, "3");

    auto clsExtent2I = declareExtent2<int>(mod, "I");
    auto clsExtent2D = declareExtent2<double>(mod, "D");

    auto clsExtent3I = declareExtent3<int>(mod, "I");
    auto clsExtent3D = declareExtent3<double>(mod, "D");

    auto clsPoint2I = declarePoint2<int>(mod, "I");
    auto clsPoint2D = declarePoint2<double>(mod, "D");

    auto clsPoint3I = declarePoint3<int>(mod, "I");
    auto clsPoint3D = declarePoint3<double>(mod, "D");

    declareExtentOperators(mod, clsExtent2I, clsExtent2D);
    declareExtentOperators(mod, clsExtent3I, clsExtent3D);
    declarePointOperators(mod, clsPoint2I, clsPoint2D);
    declarePointOperators(mod, clsPoint3I, clsPoint3D);

}

WRAP(SpherePoint) {
    //wrappers.addSignatureDependency("libcosmicRays.sphgeom");
    PySpherePoint cls(mod, "SpherePoint");
    /* Constructors */
    cls.def(py::init<>());
    cls.def(py::init<Angle const &, Angle const &>(), "longitude"_a, "latitude"_a);
    cls.def(py::init<double, double, AngleUnit>(), "longitude"_a, "latitude"_a, "units"_a);
    cls.def(py::init<sphgeom::Vector3d const &>(), "vector"_a);
    cls.def(py::init<sphgeom::UnitVector3d const &>(), "unitVector"_a);
    cls.def(py::init<sphgeom::LonLat const &>(), "lonLat"_a);
    cls.def(py::init<SpherePoint const &>(), "other"_a);
    py::implicitly_convertible<SpherePoint, sphgeom::LonLat>();
    py::implicitly_convertible<sphgeom::LonLat, SpherePoint>();

    /* Operators */
    cls.def("__getitem__",
            [](SpherePoint const &self, std::ptrdiff_t i) {
        return self[utils::python::cppIndex(2, i)];
    });
    cls.def("__eq__", &SpherePoint::operator==, py::is_operator());
    cls.def("__ne__", &SpherePoint::operator!=, py::is_operator());

    /* Members */
    cls.def("getLongitude", &SpherePoint::getLongitude);
    cls.def("getLatitude", &SpherePoint::getLatitude);
    cls.def("getRa", &SpherePoint::getRa);
    cls.def("getDec", &SpherePoint::getDec);
    cls.def("getVector", &SpherePoint::getVector);
    cls.def("getPosition", &SpherePoint::getPosition, "units"_a);
    cls.def("atPole", &SpherePoint::atPole);
    cls.def("isFinite", &SpherePoint::isFinite);
    cls.def("bearingTo", &SpherePoint::bearingTo, "other"_a);
    cls.def("separation", &SpherePoint::separation, "other"_a);
    cls.def("rotated", &SpherePoint::rotated, "axis"_a, "amount"_a);
    cls.def("offset", &SpherePoint::offset, "bearing"_a, "amount"_a);
    cls.def("getTangentPlaneOffset", &SpherePoint::getTangentPlaneOffset, "other"_a);
    utils::python::addOutputOp(cls, "__str__");
    cls.def("__len__", [](SpherePoint const &) { return 2; });
    cls.def("__reduce__", [cls](SpherePoint const &self) {
        return py::make_tuple(cls, py::make_tuple(py::cast(self.getLongitude()),
                                                  py::cast(self.getLatitude())));
    });

    /* Module level */
    mod.def("averageSpherePoint", averageSpherePoint);

}

WRAP(Interval) {
    py::class_<IntervalI, std::shared_ptr<IntervalI> > Icls(mod, "IntervalI");
    py::enum_<IntervalI::EdgeHandlingEnum>(Icls, "EdgeHandlingEnum")
            .value("EXPAND", IntervalI::EdgeHandlingEnum::EXPAND)
            .value("SHRINK", IntervalI::EdgeHandlingEnum::SHRINK);
    Icls.def(py::init<IntervalD const &, IntervalI::EdgeHandlingEnum>(), "other"_a,
            "edgeHandling"_a = IntervalI::EdgeHandlingEnum::EXPAND);
    Icls.def("getBegin", &IntervalI::getBegin);
    Icls.def_property_readonly("begin", &IntervalI::getBegin);
    Icls.def("getEnd", &IntervalI::getEnd);
    Icls.def_property_readonly("end", &IntervalI::getEnd);
    declareCommonIntervalInterface(Icls);

    py::class_<IntervalD, std::shared_ptr<IntervalD> > Dcls(mod, "IntervalD");
    Dcls.def(py::init<IntervalI const &>());
    Dcls.def("getCenter", &IntervalD::getCenter);
    Dcls.def_property_readonly("center", &IntervalD::getCenter);
    Dcls.def("isFinite", &IntervalD::isFinite);
    declareCommonIntervalInterface(Dcls);
}

WRAP(Box) {
    py::class_<Box2I, std::shared_ptr<Box2I> > Icls(mod, "Box2I");
    Icls.attr("Point") = mod.attr("Point2I");
    Icls.attr("Extent") = mod.attr("Extent2I");

    py::enum_<Box2I::EdgeHandlingEnum>(Icls, "EdgeHandlingEnum")
            .value("EXPAND", Box2I::EdgeHandlingEnum::EXPAND)
            .value("SHRINK", Box2I::EdgeHandlingEnum::SHRINK)
            .export_values();

    Icls.def(py::init<>());
    Icls.def(py::init<Point2I const &, Point2I const &, bool>(), "minimum"_a, "maximum"_a,
            "invert"_a = true);
    Icls.def(py::init<Point2I const &, Extent2I const &, bool>(), "corner"_a, "dimensions"_a,
            "invert"_a = true);
    Icls.def(py::init<IntervalI const &, IntervalI const &>(), "x"_a, "y"_a);
    Icls.def(py::init<Box2D const &, Box2I::EdgeHandlingEnum>(), "other"_a,
            "edgeHandling"_a = Box2I::EXPAND);
    Icls.def(py::init<Box2I const &>(), "other"_a);

    Icls.def("__eq__", [](Box2I const &self, Box2I const &other) { return self == other; },
    py::is_operator());
    Icls.def("__ne__", [](Box2I const &self, Box2I const &other) { return self != other; },
    py::is_operator());

    Icls.def_static("makeCenteredBox", &Box2I::makeCenteredBox, "center"_a, "size"_a);
    Icls.def("swap", &Box2I::swap);
    Icls.def("getMin", &Box2I::getMin);
    Icls.def("getMinX", &Box2I::getMinX);
    Icls.def("getMinY", &Box2I::getMinY);
    Icls.def("getMax", &Box2I::getMax);
    Icls.def("getMaxX", &Box2I::getMaxX);
    Icls.def("getMaxY", &Box2I::getMaxY);
    Icls.def_property_readonly("minX", &Box2I::getMinX);
    Icls.def_property_readonly("minY", &Box2I::getMinY);
    Icls.def_property_readonly("maxX", &Box2I::getMaxX);
    Icls.def_property_readonly("maxY", &Box2I::getMaxY);
    Icls.def("getBegin", &Box2I::getBegin);
    Icls.def("getBeginX", &Box2I::getBeginX);
    Icls.def("getBeginY", &Box2I::getBeginY);
    Icls.def("getEnd", &Box2I::getEnd);
    Icls.def("getEndX", &Box2I::getEndX);
    Icls.def("getEndY", &Box2I::getEndY);
    Icls.def_property_readonly("beginX", &Box2I::getBeginX);
    Icls.def_property_readonly("beginY", &Box2I::getBeginY);
    Icls.def_property_readonly("endX", &Box2I::getEndX);
    Icls.def_property_readonly("endY", &Box2I::getEndY);
    Icls.def("getDimensions", &Box2I::getDimensions);
    Icls.def("getWidth", &Box2I::getWidth);
    Icls.def("getHeight", &Box2I::getHeight);
    Icls.def("getArea", &Box2I::getArea);
    Icls.def_property_readonly("width", &Box2I::getWidth);
    Icls.def_property_readonly("height", &Box2I::getHeight);
    Icls.def_property_readonly("area", &Box2I::getArea);
    Icls.def("getCenter", &Box2I::getCenter);
    Icls.def("getCenterX", &Box2I::getCenterX);
    Icls.def("getCenterY", &Box2I::getCenterY);
    Icls.def_property_readonly("centerX", &Box2I::getCenterX);
    Icls.def_property_readonly("centerY", &Box2I::getCenterY);
    Icls.def("getX", &Box2I::getX);
    Icls.def("getY", &Box2I::getY);
    Icls.def_property_readonly("x", &Box2I::getX);
    Icls.def_property_readonly("y", &Box2I::getY);
    Icls.def("isEmpty", &Box2I::isEmpty);
    Icls.def("contains", py::overload_cast<Point2I const &>(&Box2I::contains, py::const_));
    Icls.def("contains", py::overload_cast<Box2I const &>(&Box2I::contains, py::const_));
    Icls.def("contains",
            py::vectorize(static_cast<bool (Box2I::*)(int x, int y) const>(&Box2I::contains)),
            "x"_a, "y"_a);
    Icls.def("__contains__", py::overload_cast<Point2I const &>(&Box2I::contains, py::const_));
    Icls.def("__contains__", py::overload_cast<Box2I const &>(&Box2I::contains, py::const_));
    Icls.def("overlaps", &Box2I::overlaps);
    Icls.def("intersects", &Box2I::intersects);
    Icls.def("isDisjointFrom", &Box2I::isDisjointFrom);
    Icls.def("grow", py::overload_cast<int>(&Box2I::grow));
    Icls.def("grow", py::overload_cast<Extent2I const &>(&Box2I::grow));
    Icls.def("shift", &Box2I::shift);
    Icls.def("flipLR", &Box2I::flipLR);
    Icls.def("flipTB", &Box2I::flipTB);
    Icls.def("include", py::overload_cast<Point2I const &>(&Box2I::include));
    Icls.def("include", py::overload_cast<Box2I const &>(&Box2I::include));
    Icls.def("clip", &Box2I::clip);
    Icls.def("dilatedBy", py::overload_cast<int>(&Box2I::dilatedBy, py::const_));
    Icls.def("dilatedBy", py::overload_cast<Extent2I const &>(&Box2I::dilatedBy, py::const_));
    Icls.def("erodedBy", py::overload_cast<int>(&Box2I::erodedBy, py::const_));
    Icls.def("erodedBy", py::overload_cast<Extent2I const &>(&Box2I::erodedBy, py::const_));
    Icls.def("shiftedBy", &Box2I::shiftedBy);
    Icls.def("reflectedAboutX", &Box2I::reflectedAboutX);
    Icls.def("reflectedAboutY", &Box2I::reflectedAboutY);
    Icls.def("expandedTo", py::overload_cast<Point2I const &>(&Box2I::expandedTo, py::const_));
    Icls.def("expandedTo", py::overload_cast<Box2I const &>(&Box2I::expandedTo, py::const_));
    Icls.def("clippedTo", &Box2I::clippedTo);
    Icls.def("getCorners", &Box2I::getCorners);
    Icls.def("toString", &Box2I::toString);
    Icls.def("__repr__", [](Box2I const &self) {
        return py::str("Box2I(minimum={}, dimensions={})")
                .format(py::repr(py::cast(self.getMin())), py::repr(py::cast(self.getDimensions())));
    });
    Icls.def("__str__", [](Box2I const &self) {
        return py::str("(minimum={}, maximum={})")
                .format(py::str(py::cast(self.getMin())), py::str(py::cast(self.getMax())));
    });
    Icls.def("__reduce__", [Icls](Box2I const &self) {
        return py::make_tuple(Icls, make_tuple(py::cast(self.getMin()), py::cast(self.getMax())));
    });
    auto getSlices = [](Box2I const &self) {
        return py::make_tuple(py::slice(self.getBeginY(), self.getEndY(), 1),
                              py::slice(self.getBeginX(), self.getEndX(), 1));
    };
    Icls.def("getSlices", getSlices);
    Icls.def_property_readonly("slices", getSlices);

    mod.attr("BoxI") = Icls;


    py::class_<Box2D, std::shared_ptr<Box2D> > Dcls(mod, "Box2D");
    Dcls.attr("Point") = mod.attr("Point2D");
    Dcls.attr("Extent") = mod.attr("Extent2D");

    Dcls.attr("EPSILON") = py::float_(Box2D::EPSILON);
    Dcls.attr("INVALID") = py::float_(Box2D::INVALID);

    Dcls.def(py::init<>());
    Dcls.def(py::init<Point2D const &, Point2D const &, bool>(), "minimum"_a, "maximum"_a,
            "invert"_a = true);
    Dcls.def(py::init<Point2D const &, Extent2D const &, bool>(), "corner"_a, "dimensions"_a,
            "invert"_a = true);
    Dcls.def(py::init<IntervalD const &, IntervalD const &>(), "x"_a, "y"_a);
    Dcls.def(py::init<Box2I const &>());
    Dcls.def(py::init<Box2D const &>());

    Dcls.def("__eq__", [](Box2D const &self, Box2D const &other) { return self == other; },
    py::is_operator());
    Dcls.def("__ne__", [](Box2D const &self, Box2D const &other) { return self != other; },
    py::is_operator());

    Dcls.def_static("makeCenteredBox", &Box2D::makeCenteredBox, "center"_a, "size"_a);
    Dcls.def("swap", &Box2D::swap);
    Dcls.def("getMin", &Box2D::getMin);
    Dcls.def("getMinX", &Box2D::getMinX);
    Dcls.def("getMinY", &Box2D::getMinY);
    Dcls.def("getMax", &Box2D::getMax);
    Dcls.def("getMaxX", &Box2D::getMaxX);
    Dcls.def("getMaxY", &Box2D::getMaxY);
    Dcls.def_property_readonly("minX", &Box2D::getMinX);
    Dcls.def_property_readonly("minY", &Box2D::getMinY);
    Dcls.def_property_readonly("maxX", &Box2D::getMaxX);
    Dcls.def_property_readonly("maxY", &Box2D::getMaxY);
    Dcls.def("getDimensions", &Box2D::getDimensions);
    Dcls.def("getWidth", &Box2D::getWidth);
    Dcls.def("getHeight", &Box2D::getHeight);
    Dcls.def("getArea", &Box2D::getArea);
    Dcls.def_property_readonly("width", &Box2D::getWidth);
    Dcls.def_property_readonly("height", &Box2D::getHeight);
    Dcls.def_property_readonly("area", &Box2D::getArea);
    Dcls.def("getX", &Box2D::getX);
    Dcls.def("getY", &Box2D::getY);
    Dcls.def_property_readonly("x", &Box2D::getX);
    Dcls.def_property_readonly("y", &Box2D::getY);
    Dcls.def("getCenter", &Box2D::getCenter);
    Dcls.def("getCenterX", &Box2D::getCenterX);
    Dcls.def("getCenterY", &Box2D::getCenterY);
    Dcls.def_property_readonly("centerX", &Box2D::getCenterX);
    Dcls.def_property_readonly("centerY", &Box2D::getCenterY);
    Dcls.def("isEmpty", &Box2D::isEmpty);
    Dcls.def("contains", py::overload_cast<Point2D const &>(&Box2D::contains, py::const_));
    Dcls.def("contains", py::overload_cast<Box2D const &>(&Box2D::contains, py::const_));
    Dcls.def("contains",
            py::vectorize(static_cast<bool (Box2D::*)(double x, double y) const>(&Box2D::contains)),
            "x"_a, "y"_a);
    Dcls.def("__contains__", py::overload_cast<Point2D const &>(&Box2D::contains, py::const_));
    Dcls.def("__contains__", py::overload_cast<Box2D const &>(&Box2D::contains, py::const_));
    Dcls.def("intersects", &Box2D::intersects);
    Dcls.def("isDisjointFrom", &Box2D::isDisjointFrom);
    Dcls.def("overlaps", &Box2D::overlaps);
    Dcls.def("grow", py::overload_cast<double>(&Box2D::grow));
    Dcls.def("grow", py::overload_cast<Extent2D const &>(&Box2D::grow));
    Dcls.def("shift", &Box2D::shift);
    Dcls.def("flipLR", &Box2D::flipLR);
    Dcls.def("flipTB", &Box2D::flipTB);
    Dcls.def("include", py::overload_cast<Point2D const &>(&Box2D::include));
    Dcls.def("include", py::overload_cast<Box2D const &>(&Box2D::include));
    Dcls.def("clip", &Box2D::clip);
    Dcls.def("dilatedBy", py::overload_cast<double>(&Box2D::dilatedBy, py::const_));
    Dcls.def("dilatedBy", py::overload_cast<Extent2D const &>(&Box2D::dilatedBy, py::const_));
    Dcls.def("erodedBy", py::overload_cast<double>(&Box2D::erodedBy, py::const_));
    Dcls.def("erodedBy", py::overload_cast<Extent2D const &>(&Box2D::erodedBy, py::const_));
    Dcls.def("shiftedBy", &Box2D::shiftedBy);
    Dcls.def("reflectedAboutX", &Box2D::reflectedAboutX);
    Dcls.def("reflectedAboutY", &Box2D::reflectedAboutY);
    Dcls.def("expandedTo", py::overload_cast<Point2D const &>(&Box2D::expandedTo, py::const_));
    Dcls.def("expandedTo", py::overload_cast<Box2D const &>(&Box2D::expandedTo, py::const_));
    Dcls.def("clippedTo", &Box2D::clippedTo);
    Dcls.def("getCorners", &Box2D::getCorners);
    Dcls.def("toString", &Box2D::toString);
    Dcls.def("__repr__", [](Box2D const &self) {
        return py::str("Box2D(minimum={}, dimensions={})")
                .format(py::repr(py::cast(self.getMin())), py::repr(py::cast(self.getDimensions())));
    });
    Dcls.def("__str__", [](Box2D const &self) {
        return py::str("(minimum={}, maximum={})")
                .format(py::str(py::cast(self.getMin())), py::str(py::cast(self.getMax())));
    });
    Dcls.def("__reduce__", [Dcls](Box2D const &self) {
        return py::make_tuple(Dcls, make_tuple(py::cast(self.getMin()), py::cast(self.getMax())));
    });

    mod.attr("BoxD") = Dcls;
}

WRAP(LinearTransform) {
    py::class_<LinearTransform, std::shared_ptr<LinearTransform> > cls(mod, "LinearTransform");
    // Parameters enum is really only used as integer constants.
    cls.attr("XX") = py::cast(int(LinearTransform::Parameters::XX));
    cls.attr("YX") = py::cast(int(LinearTransform::Parameters::YX));
    cls.attr("XY") = py::cast(int(LinearTransform::Parameters::XY));
    cls.attr("YY") = py::cast(int(LinearTransform::Parameters::YY));

    /* Constructors */
    cls.def(py::init<>());
    cls.def(py::init<LinearTransform::Matrix const &>(), "matrix"_a);

    /* Operators */
    cls.def("__call__",
            py::overload_cast<Point2D const &>(&LinearTransform::operator(), py::const_));
    cls.def("__call__",
            py::overload_cast<Extent2D const &>(&LinearTransform::operator(), py::const_));
    cls.def("__call__",
            // We use pybind11's wrappers for the Python C API to
            // delegate to other wrapped methods because:
            //  - defining this in pure Python is tricky because it's
            //    an overload, not a standalone method;
            //  - we'd rather not add a new pure-Python file just for
            //    this;
            //  - using py::vectorize internal to the method would
            //    involve defining a new internal callable every time
            //    this method is called.
            // The other viable alternative would be to define
            // applyX and applyY as Python callables with py::vectorize
            // outside the lambda as C++ local variables, and then
            // capture them by value in the lambda.  This just seems
            // slightly cleaner, as it's closer to how one would
            // implement this in pure Python, if it wasn't an overload.
            [](py::object self, py::object x, py::object y) {
        return py::make_tuple(self.attr("applyX")(x, y),
                              self.attr("applyY")(x, y));
    },
    "x"_a, "y"_a);
    cls.def("__getitem__",
            [](LinearTransform const &self, int i) { return self[utils::python::cppIndex(4, i)]; });
    cls.def("__getitem__", [](LinearTransform const &self, std::pair<int, int> i) {
        auto row = utils::python::cppIndex(2, i.first);
        auto col = utils::python::cppIndex(2, i.second);
        return self.getMatrix()(row, col);
    });
    cls.def("__mul__", &LinearTransform::operator*, py::is_operator());
    cls.def("__add__", &LinearTransform::operator+, py::is_operator());
    cls.def("__sub__", &LinearTransform::operator-, py::is_operator());
    cls.def("__iadd__", &LinearTransform::operator+=);
    cls.def("__isub__", &LinearTransform::operator-=);

    /* Members */
    cls.def_static("makeScaling",
                   py::overload_cast<double>(&LinearTransform::makeScaling),
                   "scale"_a);
    cls.def_static("makeScaling",
                   py::overload_cast<double, double>(&LinearTransform::makeScaling));
    cls.def_static("makeRotation",
                   py::overload_cast<Angle>(LinearTransform::makeRotation),
                   "angle"_a);
    cls.def("getParameterVector", &LinearTransform::getParameterVector);
    cls.def("getMatrix",
            py::overload_cast<>(& LinearTransform::getMatrix, py::const_));
    cls.def("inverted", &LinearTransform::inverted);
    cls.def("computeDeterminant", &LinearTransform::computeDeterminant);
    cls.def("isIdentity", &LinearTransform::isIdentity);
    cls.def("applyX", py::vectorize(&LinearTransform::applyX), "x"_a, "y"_a);
    cls.def("applyY", py::vectorize(&LinearTransform::applyY), "x"_a, "y"_a);

    cls.def("set",
            [](LinearTransform &self, double xx, double yx, double xy, double yy) {
        self[LinearTransform::XX] = xx;
        self[LinearTransform::XY] = xy;
        self[LinearTransform::YX] = yx;
        self[LinearTransform::YY] = yy;
    },
    "xx"_a, "yx"_a, "xy"_a, "yy"_a);

    cls.def("__str__", [](LinearTransform const &self) {
        return py::str(py::cast(self.getMatrix()));
    });
    cls.def("__repr__", [](LinearTransform const &self) {
        return py::str("LinearTransform(\n{}\n)").format(py::cast(self.getMatrix()));
    });
    cls.def("__reduce__", [cls](LinearTransform const &self) {
                return py::make_tuple(cls, py::make_tuple(py::cast(self.getMatrix())));
    });
}

WRAP(AffineTransform) {
    py::class_<AffineTransform, std::shared_ptr<AffineTransform> > cls(mod, "AffineTransform");
    // Parameters enum is really only used as integer constants.
    cls.attr("XX") = py::cast(int(AffineTransform::Parameters::XX));
    cls.attr("YX") = py::cast(int(AffineTransform::Parameters::YX));
    cls.attr("XY") = py::cast(int(AffineTransform::Parameters::XY));
    cls.attr("YY") = py::cast(int(AffineTransform::Parameters::YY));
    cls.attr("X") = py::cast(int(AffineTransform::Parameters::X));
    cls.attr("Y") = py::cast(int(AffineTransform::Parameters::Y));

    /* Constructors */
    cls.def(py::init<>());
    cls.def(py::init<Eigen::Matrix3d const &>(), "matrix"_a);
    cls.def(py::init<Eigen::Matrix2d const &>(), "linear"_a);
    cls.def(py::init<Eigen::Vector2d const &>(), "translation"_a);
    cls.def(py::init<Eigen::Matrix2d const &, Eigen::Vector2d const &>(),
            "linear"_a, "translation"_a);
    cls.def(py::init<LinearTransform const &>(), "linear"_a);
    cls.def(py::init<Extent2D const &>(), "translation"_a);
    cls.def(py::init<LinearTransform const &, Extent2D const &>(), "linear"_a, "translation"_a);

    /* Operators and special methods */
    cls.def("__mul__", &AffineTransform::operator*, py::is_operator());
    cls.def("__call__",
            py::overload_cast<Point2D const &>(&AffineTransform::operator(), py::const_));
    cls.def("__call__",
            py::overload_cast<Extent2D const &>(&AffineTransform::operator(), py::const_));
    cls.def("__call__",
            // We use pybind11's wrappers for the Python C API to
            // delegate to other wrapped methods because:
            //  - defining this in pure Python is tricky because it's
            //    an overload, not a standalone method;
            //  - we'd rather not add a new pure-Python file just for
            //    this;
            //  - using py::vectorize internal to the method would
            //    involve defining a new internal callable every time
            //    this method is called.
            // The other viable alternative would be to define
            // applyX and applyY as Python callables with py::vectorize
            // outside the lambda as C++ local variables, and then
            // capture them by value in the lambda.  This just seems
            // slightly cleaner, as it's closer to how one would
            // implement this in pure Python, if it wasn't an overload.
            [](py::object self, py::object x, py::object y) {
        return py::make_tuple(self.attr("applyX")(x, y),
                              self.attr("applyY")(x, y));
    },
    "x"_a, "y"_a);
    cls.def("__setitem__", [](AffineTransform &self, int i, double value) {
        if (i < 0 || i > 5) {
            PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d", i);
            throw py::error_already_set();
        }
        self[i] = value;
    });
    cls.def("__getitem__", [](AffineTransform const &self, int row, int col) {
        if (row < 0 || row > 2 || col < 0 || col > 2) {
            PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d, %d", row, col);
            throw py::error_already_set();
        }
        return (self.getMatrix())(row, col);
    });
    cls.def("__getitem__", [](AffineTransform const &self, int i) {
        if (i < 0 || i > 5) {
            PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d", i);
            throw py::error_already_set();
        }
        return self[i];
    });
    cls.def("__str__", [](AffineTransform const &self) {
        return py::str(py::cast(self.getMatrix())); }
    );
    cls.def("__repr__", [](AffineTransform const &self) {
        return py::str("AffineTransform(\n{}\n)").format(py::cast(self.getMatrix()));
    });
    cls.def("__reduce__", [cls](AffineTransform const &self) {
        return py::make_tuple(cls, py::make_tuple(py::cast(self.getMatrix())));
    });

    /* Members */
    cls.def("inverted", &AffineTransform::inverted);
    cls.def("isIdentity", &AffineTransform::isIdentity);
    cls.def("getTranslation", (Extent2D & (AffineTransform::*)()) & AffineTransform::getTranslation);
    cls.def("getLinear", (LinearTransform & (AffineTransform::*)()) & AffineTransform::getLinear);
    cls.def("getMatrix", &AffineTransform::getMatrix);
    cls.def("getParameterVector", &AffineTransform::getParameterVector);
    cls.def("setParameterVector", &AffineTransform::setParameterVector);
    cls.def("applyX", py::vectorize(&AffineTransform::applyX), "x"_a, "y"_a);
    cls.def("applyY", py::vectorize(&AffineTransform::applyY), "x"_a, "y"_a);
    cls.def_static("makeScaling", py::overload_cast<double>(&AffineTransform::makeScaling));
    cls.def_static("makeScaling", py::overload_cast<double, double>(&AffineTransform::makeScaling));
    cls.def_static("makeRotation", &AffineTransform::makeRotation, "angle"_a);
    cls.def_static("makeTranslation", &AffineTransform::makeTranslation, "translation"_a);

    /* Non-members */
    mod.def("makeAffineTransformFromTriple", makeAffineTransformFromTriple);
}

}
WRAP(Geom) {
    wrapAngle(mod);
    wrapCoordinates(mod);
    wrapSpherePoint(mod);
    wrapInterval(mod);
    wrapBox(mod);
    wrapLinearTransform(mod);
    wrapAffineTransform(mod);

}

}
}
