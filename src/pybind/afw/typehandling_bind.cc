#include <memory>
#include <cstdint>
#include <exception>
#include <iostream>
#include <utility>

#include "pybind/afw_bind.h"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/typehandling/Storable.h"
#include "lsst/afw/typehandling/python.h"
#include "lsst/utils/python/PySharedPtr.h"
#include "lsst/afw/typehandling/SimpleGenericMap.h"
#include "lsst/afw/typehandling/GenericMap.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace python = lsst::utils::python;

namespace pybind11 {
namespace detail {
template <typename... Ts>
struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>> {};
template <>
struct visit_helper<boost::variant> {
    template <typename... Args>
    static auto call(Args&&... args) -> decltype(boost::apply_visitor(args...)) {
        return boost::apply_visitor(args...);
    }
};
}  // namespace detail
}  // namespace pybind11

namespace lsst {
namespace afw {
namespace typehandling {

namespace {
template <typename K>
void declareSimpleGenericMap(py::module_ &mod, std::string const& suffix,
                             std::string const& key) {
    using Class = SimpleGenericMap<K>;

    std::string className = "SimpleGenericMap" + suffix;
    // Give the class a custom docstring to avoid confusing Python users
    std::string docstring =
            "A `dict`-like `~collections.abc.MutableMapping` for use when sharing a map between C++ and "
            "Python.\n" +
            declareGenericMapRestrictions(className, key) +
            R"docstring(
Parameters
----------
mapping : `collections.abc.Mapping`, optional
iterable : iterable, optional
**kwargs
    A ``SimpleGenericMap`` takes the same input arguments as `dict`.
)docstring";
    py::class_<Class, std::shared_ptr<Class>, MutableGenericMap<K>>(mod, className.c_str(), docstring.c_str())
                          // Don't rewrap members of MutableGenericMap

                          /* need __init__(**kw), __init__(mapping, **kw), __init__(iterable, **kw)
                           * can't find a good way to insert a py::handle value from C++
                           * can't call MutableGenericMap.__setitem__ without its class_ object, which
                           *    is in a different pybind11 module but can't be imported (yet)
                           */
            .def(py::init<>())
            .def("copy", [](Class const& self) { return Class(self); });

                          // fromkeys easier to implement in Python

    //wrappers.wrapType(PyClass(wrappers.module, className.c_str(), docstring.c_str()),
    //                  [](auto& mod, auto& cls) {
                          // Don't rewrap members of MutableGenericMap

                          /* need __init__(**kw), __init__(mapping, **kw), __init__(iterable, **kw)
                           * can't find a good way to insert a py::handle value from C++
                           * can't call MutableGenericMap.__setitem__ without its class_ object, which
                           *    is in a different pybind11 module but can't be imported (yet)
                           */
     //                     cls.def(py::init<>());
     //                     cls.def("copy", [](Class const& self) { return Class(self); });

                          // fromkeys easier to implement in Python
     //                 });
}

// Type safety pointless in Python, use unsafe methods to avoid manual type checking
// https://pybind11.readthedocs.io/en/stable/advanced/classes.html#binding-protected-member-functions
template <typename K>
class Publicist : public MutableGenericMap<K> {
public:
    using typename GenericMap<K>::ConstValueReference;
    using GenericMap<K>::unsafeLookup;
    using MutableGenericMap<K>::unsafeErase;
};

template <typename K>
py::object get(GenericMap<K>& self, K const& key) {
    auto callable = static_cast<typename Publicist<K>::ConstValueReference (GenericMap<K>::*)(K) const>(
            &Publicist<K>::unsafeLookup);
    auto variant = (self.*callable)(key);

    // py::cast can't convert PolymorphicValue to Storable; do it by hand
    PolymorphicValue const* storableHolder = boost::get<PolymorphicValue const&>(&variant);
    if (storableHolder) {
        Storable const& value = *storableHolder;
        // Prevent segfaults when assigning a Key<Storable> to Python variable, then deleting from map
        // No existing code depends on being able to modify an item stored by value
        return py::cast(value.cloneStorable());
    } else {
        return py::cast(variant);
    }
};

template <typename K>
void declareGenericMap(py::module_ &mod, std::string const& suffix,
                       std::string const& key) {
    using Class = GenericMap<K>;

    std::string className = "GenericMap" + suffix;
    // Give the class a custom docstring to avoid confusing Python users
    std::string docstring =
            "An abstract `~collections.abc.Mapping` for use when sharing a map between C++ and Python.\n" +
            declareGenericMapRestrictions(className, key);
    py::class_<Class, std::shared_ptr<Class>>(mod, className.c_str(), docstring.c_str())
        // __str__ easier to implement in Python
        // __repr__ easier to implement in Python
        // __eq__ easier to implement in Python
        // __ne__ easier to implement in Python
        // For unknown reasons, overload_cast doesn't work
        // cls.def("__contains__", py::overload_cast<K const&>(&Class::contains, py::const_), "key"_a);
            .def("__contains__", static_cast<bool (Class::*)(K const&) const>(&Class::contains), "key"_a)

            .def("__getitem__",
                 [](Class& self, K const& key) {
                    try {
                        return get(self, key);
                    } catch (pex::exceptions::OutOfRangeError const& e) {
                        // pybind11 doesn't seem to recognize chained exceptions
                        std::stringstream buffer;
                        buffer << "Unknown key: " << key;
                        std::throw_with_nested(py::key_error(buffer.str()));
                    }
                },
                "key"_a)
            .def("get",
                [](Class& self, K const& key, py::object const& def) {
                    try {
                        return get(self, key);
                    } catch (pex::exceptions::OutOfRangeError const& e) {
                        return def;
                    }
                },
                // Prevent segfaults when assigning a key<Storable> to Python variable, then deleting from map
                // No existing code depends on being able to modify an item stored by value
                "key"_a, "default"_a = py::none(), py::return_value_policy::copy)
           .def("__iter__",
                [](Class const& self) { return py::make_iterator(self.keys().begin(), self.keys().end()); },
                py::keep_alive<0, 1>())
           .def("__len__", &Class::size)
           .def("__bool__", [](Class const& self) { return !self.empty(); });
        // Can't wrap keys directly because pybind11 always copies vectors, so it won't be a view
        // items easier to implement in Python
        // values easier to implement in Python

    /*wrappers.wrapType(PyClass(wrappers.module, className.c_str(), docstring.c_str()), [](auto& mod,
                                                                                         auto& cls) {
        // __str__ easier to implement in Python
        // __repr__ easier to implement in Python
        // __eq__ easier to implement in Python
        // __ne__ easier to implement in Python
        // For unknown reasons, overload_cast doesn't work
        // cls.def("__contains__", py::overload_cast<K const&>(&Class::contains, py::const_), "key"_a);
        cls.def("__contains__", static_cast<bool (Class::*)(K const&) const>(&Class::contains), "key"_a)

           .def("__getitem__",
                [](Class& self, K const& key) {
                    try {
                        return get(self, key);
                    } catch (pex::exceptions::OutOfRangeError const& e) {
                        // pybind11 doesn't seem to recognize chained exceptions
                        std::stringstream buffer;
                        buffer << "Unknown key: " << key;
                        std::throw_with_nested(py::key_error(buffer.str()));
                    }
                },
                "key"_a)
           .def("get",
                [](Class& self, K const& key, py::object const& def) {
                    try {
                        return get(self, key);
                    } catch (pex::exceptions::OutOfRangeError const& e) {
                        return def;
                    }
                },
                // Prevent segfaults when assigning a key<Storable> to Python variable, then deleting from map
                // No existing code depends on being able to modify an item stored by value
                "key"_a, "default"_a = py::none(), py::return_value_policy::copy)
           .def("__iter__",
                [](Class const& self) { return py::make_iterator(self.keys().begin(), self.keys().end()); },
                py::keep_alive<0, 1>())
           .def("__len__", &Class::size)
           .def("__bool__", [](Class const& self) { return !self.empty(); });
        // Can't wrap keys directly because pybind11 always copies vectors, so it won't be a view
        // items easier to implement in Python
        // values easier to implement in Python
    });*/

}

template <typename V, class PyClass>
void declareMutableGenericMapTypedMethods(PyClass& cls) {
    using Class = typename PyClass::type;
    cls.def("__setitem__",
            [](Class& self, typename Class::key_type const& key, V const& value) {
                // Need to delete previous key, which may not be of type V
                // TODO: this method provides only basic exception safety
                if (self.contains(key)) {
                    auto callable = &Publicist<typename Class::key_type>::unsafeErase;
                    (self.*callable)(key);
                }
                self.insert(key, value);
            },
            "key"_a, "value"_a);
    // setdefault easier to implement in Python
    // pop easier to implement in Python
}

template <typename K>
void declareMutableGenericMap(py::module_ &mod, std::string const& suffix,
                              std::string const& key) {
    using Class = MutableGenericMap<K>;

    std::string className = "MutableGenericMap" + suffix;
    // Give the class a custom docstring to avoid confusing Python users
    std::string docstring =
            "An abstract `~collections.abc.MutableMapping` for use when sharing a map between C++ and "
            "Python.\n" +
            declareGenericMapRestrictions(className, key);
    py::class_<Class, std::shared_ptr<Class>, GenericMap<K>> cls(mod, className.c_str(), docstring.c_str());
    // Don't rewrap members of GenericMap
    declareMutableGenericMapTypedMethods<std::shared_ptr<Storable const>>(cls);
    declareMutableGenericMapTypedMethods<bool>(cls);
    // TODO: int32 and float are suppressed for now, see DM-21268
    declareMutableGenericMapTypedMethods<std::int64_t>(cls);  // chosen for builtins.int
    declareMutableGenericMapTypedMethods<std::int32_t>(cls);
    declareMutableGenericMapTypedMethods<double>(cls);  // chosen for builtins.float
    declareMutableGenericMapTypedMethods<float>(cls);
    declareMutableGenericMapTypedMethods<std::string>(cls);
    cls.def("__delitem__",
            [](Class& self, K const& key) {
        if (self.contains(key)) {
            auto callable = &Publicist<K>::unsafeErase;
            (self.*callable)(key);
        } else {
            std::stringstream buffer;
            buffer << "Unknown key: " << key;
            throw py::key_error(buffer.str());
        }
    },
    "key"_a);
    cls.def("popitem", [](Class& self) {
        if (!self.empty()) {
            K key = self.keys().back();
            auto result = std::make_pair(key, get(self, key));
            auto callable = &Publicist<K>::unsafeErase;
            (self.*callable)(key);
            return result;
        } else {
            throw py::key_error("Cannot pop from empty GenericMap.");
        }
    });
    cls.def("clear", &Class::clear);
    // cls.def("update",, "other"_a);

    /*wrappers.wrapType(PyClass(wrappers.module, className.c_str(), docstring.c_str()),
                      [](auto& mod, auto& cls) {
                          // Don't rewrap members of GenericMap
                          declareMutableGenericMapTypedMethods<std::shared_ptr<Storable const>>(cls);
                          declareMutableGenericMapTypedMethods<bool>(cls);
                          // TODO: int32 and float are suppressed for now, see DM-21268
                          declareMutableGenericMapTypedMethods<std::int64_t>(cls);  // chosen for builtins.int
                          declareMutableGenericMapTypedMethods<std::int32_t>(cls);
                          declareMutableGenericMapTypedMethods<double>(cls);  // chosen for builtins.float
                          declareMutableGenericMapTypedMethods<float>(cls);
                          declareMutableGenericMapTypedMethods<std::string>(cls);
                          cls.def("__delitem__",
                                  [](Class& self, K const& key) {
                                      if (self.contains(key)) {
                                          auto callable = &Publicist<K>::unsafeErase;
                                          (self.*callable)(key);
                                      } else {
                                          std::stringstream buffer;
                                          buffer << "Unknown key: " << key;
                                          throw py::key_error(buffer.str());
                                      }
                                  },
                                  "key"_a);
                          cls.def("popitem", [](Class& self) {
                              if (!self.empty()) {
                                  K key = self.keys().back();
                                  auto result = std::make_pair(key, get(self, key));
                                  auto callable = &Publicist<K>::unsafeErase;
                                  (self.*callable)(key);
                                  return result;
                              } else {
                                  throw py::key_error("Cannot pop from empty GenericMap.");
                              }
                          });
                          cls.def("clear", &Class::clear);
                          // cls.def("update",, "other"_a);
                      });*/

}

}

WRAP(SimpleGenericMap) {
    declareSimpleGenericMap<std::string>(mod, "S", "strings");
}

WRAP(GenericMap) {
    declareGenericMap<std::string>(mod, "S", "strings");
    declareMutableGenericMap<std::string>(mod, "S", "strings");
}

WRAP(Storable) {
    py::class_<Storable, python::PySharedPtr<Storable>, table::io::Persistable, StorableHelper<>>(mod, "Storable")
            .def(py::init<>())
            .def("__eq__", [](Storable const& self, Storable const& other) { return self.equals(other); }, "other"_a);
}

WRAP(Typehandling) {
    auto typemod = mod.def_submodule("typehandling");
    wrapStorable(typemod);
    wrapGenericMap(typemod);
    wrapSimpleGenericMap(typemod);

}
}
}
}
