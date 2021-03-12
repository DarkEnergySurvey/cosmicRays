#include "pybind/daf_bind.h"

namespace lsst {
namespace daf {

namespace {
template <typename T, typename C>
void declarePropertySetAccessors(C& cls, std::string const& name) {
    const std::string getName = "get" + name;
    cls.def(getName.c_str(), (T (base::PropertySet::*)(std::string const&) const) & base::PropertySet::get<T>, "name"_a);
    cls.def(getName.c_str(), (T (base::PropertySet::*)(std::string const&, T const&) const) & base::PropertySet::get<T>,
            "name"_a, "defaultValue"_a);

    const std::string getArrayName = "getArray" + name;
    cls.def(getArrayName.c_str(),
            (std::vector<T> (base::PropertySet::*)(std::string const&) const) & base::PropertySet::getArray<T>, "name"_a);

    const std::string setName = "set" + name;
    cls.def(setName.c_str(), (void (base::PropertySet::*)(std::string const&, T const&)) & base::PropertySet::set<T>,
            "name"_a, "value"_a);
    cls.def(setName.c_str(),
            (void (base::PropertySet::*)(std::string const&, std::vector<T> const&)) & base::PropertySet::set<T>,
            "name"_a, "value"_a);

    const std::string addName = "add" + name;
    cls.def(addName.c_str(), (void (base::PropertySet::*)(std::string const&, T const&)) & base::PropertySet::add<T>,
            "name"_a, "value"_a);
    cls.def(addName.c_str(),
            (void (base::PropertySet::*)(std::string const&, std::vector<T> const&)) & base::PropertySet::add<T>,
            "name"_a, "value"_a);

    const std::string typeName = "TYPE_" + name;
    cls.attr(typeName.c_str()) = py::cast(base::PropertySet::typeOfT<T>(), py::return_value_policy::reference);
}

template <typename T, typename C>
void declareAccessors(C& cls, std::string const& name) {
    const std::string getName = "get" + name;
    cls.def(getName.c_str(), (T (base::PropertyList::*)(std::string const&) const) & base::PropertyList::get<T>,
            "name"_a);
    cls.def(getName.c_str(), (T (base::PropertyList::*)(std::string const&, T const&) const) & base::PropertyList::get<T>,
            "name"_a, "defaultValue"_a);

    // Warning: __len__ is ambiguous so do not attempt to define it. It could return
    // the number of unique names or the number of entries (e.g. as returned by toList,
    // a pure Python method). C++ begin and end iterate over unique names, but users often
    // view PropertyList as a representation of a FITS header. When in doubt, refuse to guess.

    const std::string getArrayName = "getArray" + name;
    cls.def(getArrayName.c_str(),
            (std::vector<T> (base::PropertyList::*)(std::string const&) const) & base::PropertyList::getArray<T>,
            "name"_a);

    const std::string setName = "set" + name;
    cls.def(setName.c_str(), (void (base::PropertyList::*)(std::string const&, T const&)) & base::PropertyList::set<T>);
    cls.def(setName.c_str(),
            (void (base::PropertyList::*)(std::string const&, std::vector<T> const&)) & base::PropertyList::set<T>);
    cls.def(setName.c_str(), (void (base::PropertyList::*)(std::string const&, T const&, std::string const&)) &
                                     base::PropertyList::set<T>);
    cls.def(setName.c_str(),
            (void (base::PropertyList::*)(std::string const&, std::vector<T> const&, std::string const&)) &
                    base::PropertyList::set<T>);

    const std::string addName = "add" + name;
    cls.def(addName.c_str(), (void (base::PropertyList::*)(std::string const&, T const&)) & base::PropertyList::add<T>);
    cls.def(addName.c_str(),
            (void (base::PropertyList::*)(std::string const&, std::vector<T> const&)) & base::PropertyList::add<T>);
    cls.def(addName.c_str(), (void (base::PropertyList::*)(std::string const&, T const&, std::string const&)) &
                                     base::PropertyList::add<T>);
    cls.def(addName.c_str(),
            (void (base::PropertyList::*)(std::string const&, std::vector<T> const&, std::string const&)) &
                    base::PropertyList::add<T>);

    const std::string typeName = "TYPE_" + name;
    cls.attr(typeName.c_str()) = py::cast(typeid(T), py::return_value_policy::reference);
}

}
WRAP(Daf) {
    py::class_<base::Persistable, std::shared_ptr<base::Persistable> > (mod, "Persistable");

    py::class_<base::DateTime> dateTimeCls(mod, "DateTime");

    py::enum_<base::DateTime::Timescale>(dateTimeCls, "Timescale")
            .value("TAI", base::DateTime::Timescale::TAI)
            .value("UTC", base::DateTime::Timescale::UTC)
            .value("TT", base::DateTime::Timescale::TT)
            .export_values();

    py::enum_<base::DateTime::DateSystem>(dateTimeCls, "DateSystem")
            .value("JD", base::DateTime::DateSystem::JD)
            .value("MJD", base::DateTime::DateSystem::MJD)
            .value("EPOCH", base::DateTime::DateSystem::EPOCH)
            .export_values();

    dateTimeCls.def(py::init<>())
            .def_readonly_static("invalid_nsecs", &base::DateTime::invalid_nsecs)
            .def(py::init<long long, base::DateTime::Timescale>(), "nsecs"_a, "scale"_a = base::DateTime::Timescale::TAI)
            .def(py::init<double, base::DateTime::DateSystem, base::DateTime::Timescale>(), "date"_a,
                 "system"_a = base::DateTime::DateSystem::MJD, "scale"_a = base::DateTime::Timescale::TAI)
            .def(py::init<int, int, int, int, int, int, base::DateTime::Timescale>())
            .def(py::init<const std::string &, base::DateTime::Timescale>())
            .def("nsecs", &base::DateTime::nsecs, "scale"_a = base::DateTime::Timescale::TAI)
            .def("get", &base::DateTime::get, "system"_a = base::DateTime::DateSystem::MJD,
                 "scale"_a = base::DateTime::Timescale::TAI)
            .def("toString", &base::DateTime::toString)
            .def("gmtime", &base::DateTime::gmtime)
            .def("timespec", &base::DateTime::timespec)
            .def("timeval", &base::DateTime::timeval)
            .def("isValid", &base::DateTime::isValid)
            .def_static("now", &base::DateTime::now)
            .def_static("initializeLeapSeconds", &base::DateTime::initializeLeapSeconds)
            .def("__eq__", [](base::DateTime const &self, base::DateTime const &other) { return self == other; },
    py::is_operator());

    py::class_<std::type_info>(mod, "TypeInfo")
            .def("__eq__",
                 [](std::type_info const& self, std::type_info const& other) { return self == other; })
    .def("__ne__",
         [](std::type_info const& self, std::type_info const& other) { return self != other; })
    .def("name", &std::type_info::name)
            .def("__hash__", &std::type_info::hash_code);

    py::class_<base::PropertySet, std::shared_ptr<base::PropertySet> > propertySetCls(mod, "PropertySet");

    propertySetCls.def(py::init<bool>(), "flat"_a = false);

    propertySetCls.def("deepCopy", &base::PropertySet::deepCopy);
    propertySetCls.def("nameCount", &base::PropertySet::nameCount, "topLevelOnly"_a = true);
    propertySetCls.def("names", &base::PropertySet::names, "topLevelOnly"_a = true);
    propertySetCls.def("paramNames", &base::PropertySet::paramNames, "topLevelOnly"_a = true);
    propertySetCls.def("propertySetNames", &base::PropertySet::propertySetNames, "topLevelOnly"_a = true);
    propertySetCls.def("exists", &base::PropertySet::exists);
    propertySetCls.def("isArray", &base::PropertySet::isArray);
    propertySetCls.def("isUndefined", &base::PropertySet::isUndefined);
    propertySetCls.def("isPropertySetPtr", &base::PropertySet::isPropertySetPtr);
    propertySetCls.def("valueCount",
                       py::overload_cast<>(&base::PropertySet::valueCount, py::const_));
    propertySetCls.def("valueCount",
                       py::overload_cast<std::string const&>(&base::PropertySet::valueCount,
                                                             py::const_));
    propertySetCls.def("typeOf", &base::PropertySet::typeOf, py::return_value_policy::reference);
    propertySetCls.def("toString", &base::PropertySet::toString, "topLevelOnly"_a = false, "indent"_a = "");
    propertySetCls.def("copy", &base::PropertySet::copy, "dest"_a, "source"_a, "name"_a, "asScalar"_a=false);
    propertySetCls.def("combine", &base::PropertySet::combine);
    propertySetCls.def("remove", &base::PropertySet::remove);
    propertySetCls.def("getAsBool", &base::PropertySet::getAsBool);
    propertySetCls.def("getAsInt", &base::PropertySet::getAsInt);
    propertySetCls.def("getAsInt64", &base::PropertySet::getAsInt64);
    propertySetCls.def("getAsUInt64", &base::PropertySet::getAsUInt64);
    propertySetCls.def("getAsDouble", &base::PropertySet::getAsDouble);
    propertySetCls.def("getAsString", &base::PropertySet::getAsString);
    propertySetCls.def("getAsPropertySetPtr", &base::PropertySet::getAsPropertySetPtr);
    propertySetCls.def("getAsPersistablePtr", &base::PropertySet::getAsPersistablePtr);

    declarePropertySetAccessors<bool>(propertySetCls, "Bool");
    declarePropertySetAccessors<short>(propertySetCls, "Short");
    declarePropertySetAccessors<int>(propertySetCls, "Int");
    declarePropertySetAccessors<long>(propertySetCls, "Long");
    declarePropertySetAccessors<long long>(propertySetCls, "LongLong");
    declarePropertySetAccessors<unsigned long long>(propertySetCls, "UnsignedLongLong");
    declarePropertySetAccessors<float>(propertySetCls, "Float");
    declarePropertySetAccessors<double>(propertySetCls, "Double");
    declarePropertySetAccessors<std::nullptr_t>(propertySetCls, "Undef");
    declarePropertySetAccessors<std::string>(propertySetCls, "String");
    declarePropertySetAccessors<base::DateTime>(propertySetCls, "DateTime");
    declarePropertySetAccessors<std::shared_ptr<base::PropertySet> >(propertySetCls, "PropertySet");

    py::class_<base::PropertyList, std::shared_ptr<base::PropertyList>, base::PropertySet> propertyListCls(mod, "PropertyList");

    propertyListCls.def(py::init<>());

    propertyListCls.def("getComment", &base::PropertyList::getComment);
    propertyListCls.def("getOrderedNames", &base::PropertyList::getOrderedNames);
    propertyListCls.def("deepCopy",
                        [](base::PropertyList const& self) { return std::static_pointer_cast<base::PropertySet>(self.deepCopy()); });
    declareAccessors<bool>(propertyListCls, "Bool");
    declareAccessors<short>(propertyListCls, "Short");
    declareAccessors<int>(propertyListCls, "Int");
    declareAccessors<long>(propertyListCls, "Long");
    declareAccessors<long long>(propertyListCls, "LongLong");
    declareAccessors<float>(propertyListCls, "Float");
    declareAccessors<double>(propertyListCls, "Double");
    declareAccessors<std::nullptr_t>(propertyListCls, "Undef");
    declareAccessors<std::string>(propertyListCls, "String");
    declareAccessors<base::DateTime>(propertyListCls, "DateTime");

    propertyListCls.def("setPropertySet",
                        (void (base::PropertyList::*)(std::string const&, base::PropertySet::Ptr const&)) & base::PropertyList::set);

}

}
}
