#include "pybind/afw_bind.h"

#include <memory>
#include <sstream>

#include "ndarray/pybind11.h"
#include "pybind11/eigen.h"

#include "lsst/afw/cameraGeom/Detector.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/fits.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/geom/polygon/Polygon.h"
#include "lsst/afw/image/ApCorrMap.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/image/TransmissionCurve.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/afw/table/AliasMap.h"
#include "lsst/afw/table/BaseColumnView.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/table/BaseTable.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/Field.h"
#include "lsst/afw/table/Flag.h"
#include "lsst/afw/table/FunctorKey.h"
#include "lsst/afw/table/IdFactory.h"
#include "lsst/afw/table/Key.h"
#include "lsst/afw/table/Match.h"
#include "lsst/afw/table/Schema.h"
#include "lsst/afw/table/SchemaMapper.h"
#include "lsst/afw/table/Simple.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/aggregates.h"
#include "lsst/afw/table/arrays.h"
#include "lsst/afw/table/detail/Access.h"
#include "lsst/afw/table/detail/SchemaImpl.h"
#include "lsst/afw/table/fwd.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/afw/table/python/catalog.h"
#include "lsst/afw/table/python/columnView.h"
#include "lsst/afw/table/python/sortedCatalog.h"
#include "lsst/afw/table/slots.h"
#include "lsst/afw/table/wcsUtils.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/Box.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/SpherePoint.h"
#include "lsst/pex/config/python.h"
#include "lsst/utils/python.h"


namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace afw {
namespace table {

using PySchemaMapper = py::class_<SchemaMapper, std::shared_ptr<SchemaMapper>>;
using utils::python::WrapperCollection;

namespace {

template <typename ReferenceCollection>
void declareUpdateRefCentroids(WrapperCollection &wrappers) {
    wrappers.wrap([](auto &mod) {
        mod.def("updateRefCentroids", updateRefCentroids<ReferenceCollection>, "wcs"_a, "refList"_a);
    });
}

template <typename SourceCollection>
void declareUpdateSourceCoords(WrapperCollection &wrappers) {
    wrappers.wrap([](auto &mod) {
        mod.def("updateSourceCoords", updateSourceCoords<SourceCollection>, "wcs"_a, "sourceList"_a);
    });
}

using PySourceRecord = py::class_<SourceRecord, std::shared_ptr<SourceRecord>, SimpleRecord>;
using PySourceTable = py::class_<SourceTable, std::shared_ptr<SourceTable>, SimpleTable>;
using PySourceColumnView =
        py::class_<SourceColumnViewT<SourceRecord>, std::shared_ptr<SourceColumnViewT<SourceRecord>>,
                   ColumnViewT<SourceRecord>>;

/*
Declare member and static functions for a pybind11 wrapper of SourceRecord
*/
PySourceRecord declareSourceRecord(WrapperCollection &wrappers) {
    return wrappers.wrapType(PySourceRecord(wrappers.module, "SourceRecord"), [](auto &mod, auto &cls) {
        cls.def("getFootprint", &SourceRecord::getFootprint);
        cls.def("setFootprint", &SourceRecord::setFootprint);
        cls.def("getTable", &SourceRecord::getTable);
        cls.def_property_readonly("table", &SourceRecord::getTable);

        cls.def("getParent", &SourceRecord::getParent);
        cls.def("setParent", &SourceRecord::setParent, "id"_a);

        cls.def("getPsfInstFlux", &SourceRecord::getPsfInstFlux);
        cls.def("getPsfInstFluxErr", &SourceRecord::getPsfInstFluxErr);
        cls.def("getPsfFluxFlag", &SourceRecord::getPsfFluxFlag);

        cls.def("getModelInstFlux", &SourceRecord::getModelInstFlux);
        cls.def("getModelInstFluxErr", &SourceRecord::getModelInstFluxErr);
        cls.def("getModelFluxFlag", &SourceRecord::getModelFluxFlag);

        cls.def("getApInstFlux", &SourceRecord::getApInstFlux);
        cls.def("getApInstFluxErr", &SourceRecord::getApInstFluxErr);
        cls.def("getApFluxFlag", &SourceRecord::getApFluxFlag);

        cls.def("getGaussianInstFlux", &SourceRecord::getGaussianInstFlux);
        cls.def("getGaussianInstFluxErr", &SourceRecord::getGaussianInstFluxErr);
        cls.def("getGaussianFluxFlag", &SourceRecord::getGaussianFluxFlag);

        cls.def("getCalibInstFlux", &SourceRecord::getCalibInstFlux);
        cls.def("getCalibInstFluxErr", &SourceRecord::getCalibInstFluxErr);
        cls.def("getCalibFluxFlag", &SourceRecord::getCalibFluxFlag);

        cls.def("getCentroid", &SourceRecord::getCentroid);
        cls.def("getCentroidErr", &SourceRecord::getCentroidErr);
        cls.def("getCentroidFlag", &SourceRecord::getCentroidFlag);

        cls.def("getShape", &SourceRecord::getShape);
        cls.def("getShapeErr", &SourceRecord::getShapeErr);
        cls.def("getShapeFlag", &SourceRecord::getShapeFlag);

        cls.def("getX", &SourceRecord::getX);
        cls.def("getY", &SourceRecord::getY);
        cls.def("getIxx", &SourceRecord::getIxx);
        cls.def("getIyy", &SourceRecord::getIyy);
        cls.def("getIxy", &SourceRecord::getIxy);
        cls.def("updateCoord", (void (SourceRecord::*)(geom::SkyWcs const &)) & SourceRecord::updateCoord,
                "wcs"_a);
        cls.def("updateCoord",
                (void (SourceRecord::*)(geom::SkyWcs const &, PointKey<double> const &)) &
                        SourceRecord::updateCoord,
                "wcs"_a, "key"_a);
    });
}

/*
Declare member and static functions for a pybind11 wrapper of SourceTable
*/
PySourceTable declareSourceTable(WrapperCollection &wrappers) {
    return wrappers.wrapType(PySourceTable(wrappers.module, "SourceTable"), [](auto &mod, auto &cls) {
        cls.def("clone", &SourceTable::clone);
        cls.def_static("make",
                       (std::shared_ptr<SourceTable>(*)(Schema const &, std::shared_ptr<IdFactory> const &)) &
                               SourceTable::make);
        cls.def_static("make", (std::shared_ptr<SourceTable>(*)(Schema const &)) & SourceTable::make);
        cls.def_static("makeMinimalSchema", &SourceTable::makeMinimalSchema);
        cls.def_static("getParentKey", &SourceTable::getParentKey);
        cls.def("copyRecord", (std::shared_ptr<SourceRecord>(SourceTable::*)(BaseRecord const &)) &
                                      SourceTable::copyRecord);
        cls.def("copyRecord",
                (std::shared_ptr<SourceRecord>(SourceTable::*)(BaseRecord const &, SchemaMapper const &)) &
                        SourceTable::copyRecord);
        cls.def("makeRecord", &SourceTable::makeRecord);

        cls.def("getPsfFluxSlot", &SourceTable::getPsfFluxSlot);
        cls.def("definePsfFlux", &SourceTable::definePsfFlux, "name"_a);

        cls.def("getModelFluxSlot", &SourceTable::getModelFluxSlot);
        cls.def("defineModelFlux", &SourceTable::defineModelFlux, "name"_a);

        cls.def("getApFluxSlot", &SourceTable::getApFluxSlot);
        cls.def("defineApFlux", &SourceTable::defineApFlux, "name"_a);

        cls.def("getGaussianFluxSlot", &SourceTable::getGaussianFluxSlot);
        cls.def("defineGaussianFlux", &SourceTable::defineGaussianFlux, "name"_a);

        cls.def("getCalibFluxSlot", &SourceTable::getCalibFluxSlot);
        cls.def("defineCalibFlux", &SourceTable::defineCalibFlux, "name"_a);

        cls.def("getCentroidSlot", &SourceTable::getCentroidSlot);
        cls.def("defineCentroid", &SourceTable::defineCentroid, "name"_a);

        cls.def("getShapeSlot", &SourceTable::getShapeSlot);
        cls.def("defineShape", &SourceTable::defineShape, "name"_a);
    });
}

PySourceColumnView declareSourceColumnView(WrapperCollection &wrappers) {
    table::python::declareColumnView<SourceRecord>(wrappers, "Source", true);
    return wrappers.wrapType(PySourceColumnView(wrappers.module, "SourceColumnView"),
                             [](auto &mod, auto &cls) {
                                 using Class = SourceColumnViewT<SourceRecord>;
                                 cls.def("getPsfInstFlux", &Class::getPsfInstFlux);
                                 cls.def("getPsfInstFluxErr", &Class::getPsfInstFluxErr);
                                 cls.def("getApInstFlux", &Class::getApInstFlux);
                                 cls.def("getApInstFluxErr", &Class::getApInstFluxErr);
                                 cls.def("getModelInstFlux", &Class::getModelInstFlux);
                                 cls.def("getModelInstFluxErr", &Class::getModelInstFluxErr);
                                 cls.def("getGaussianInstFlux", &Class::getGaussianInstFlux);
                                 cls.def("getGaussianInstFluxErr", &Class::getGaussianInstFluxErr);
                                 cls.def("getCalibInstFlux", &Class::getCalibInstFlux);
                                 cls.def("getCalibInstFluxErr", &Class::getCalibInstFluxErr);
                                 cls.def("getX", &Class::getX);
                                 cls.def("getY", &Class::getY);
                                 cls.def("getIxx", &Class::getIxx);
                                 cls.def("getIyy", &Class::getIyy);
                                 cls.def("getIxy", &Class::getIxy);
                             });
}

using PySlotDefinition = py::class_<SlotDefinition>;

void declareSlotDefinition(WrapperCollection &wrappers) {
    wrappers.wrapType(PySlotDefinition(wrappers.module, "SlotDefinition"), [](auto &mod, auto &cls) {
        cls.def("getName", &SlotDefinition::getName);
        cls.def("getAlias", &SlotDefinition::getAlias);
    });
}

/*
Declare standard methods for subclasses of SlotDefinition (but not SlotDefinition itself).
*/
template <typename Class>
void declareSlotDefinitionSubclass(WrapperCollection &wrappers, std::string const &name) {
    wrappers.wrapType(py::class_<Class, SlotDefinition>(wrappers.module, name.c_str()),
                      [](auto &mod, auto &cls) {
                          cls.def(py::init<std::string const &>(), "name"_a);
                          cls.def("isValid", &Class::isValid);
                          cls.def("getMeasKey", &Class::getMeasKey);
                          cls.def("getErrKey", &Class::getErrKey);
                          cls.def("getFlagKey", &Class::getFlagKey);
                          cls.def("setKeys", &Class::setKeys, "alias"_a, "schema"_a);
                      });
}

using PySimpleTable = py::class_<SimpleTable, std::shared_ptr<SimpleTable>, BaseTable>;
using PySimpleRecord = py::class_<SimpleRecord, std::shared_ptr<SimpleRecord>, BaseRecord>;

PySimpleRecord declareSimpleRecord(WrapperCollection &wrappers) {
    return wrappers.wrapType(PySimpleRecord(wrappers.module, "SimpleRecord"), [](auto &mod, auto &cls) {
        cls.def("getId", &SimpleRecord::getId);
        cls.def("setId", &SimpleRecord::setId);
        cls.def("getTable", &SimpleRecord::getTable);
        cls.def_property_readonly("table", &SimpleRecord::getTable);
        cls.def("getCoord", &SimpleRecord::getCoord);
        cls.def("setCoord", &SimpleRecord::setCoord);
        cls.def("getRa", &SimpleRecord::getRa);
        cls.def("setRa", &SimpleRecord::setRa);
        cls.def("getDec", &SimpleRecord::getDec);
        cls.def("setDec", &SimpleRecord::setDec);
    });
}

PySimpleTable declareSimpleTable(WrapperCollection &wrappers) {
    return wrappers.wrapType(PySimpleTable(wrappers.module, "SimpleTable"), [](auto &mod, auto &cls) {
        cls.def_static("make",
                       (std::shared_ptr<SimpleTable>(*)(Schema const &, std::shared_ptr<IdFactory> const &)) &
                               SimpleTable::make);
        cls.def_static("make", (std::shared_ptr<SimpleTable>(*)(Schema const &)) & SimpleTable::make);
        cls.def_static("makeMinimalSchema", &SimpleTable::makeMinimalSchema);
        cls.def_static("checkSchema", &SimpleTable::checkSchema, "schema"_a);
        cls.def_static("getIdKey", &SimpleTable::getIdKey);
        cls.def_static("getCoordKey", &SimpleTable::getCoordKey);

        cls.def("getIdFactory", (std::shared_ptr<IdFactory>(SimpleTable::*)()) & SimpleTable::getIdFactory);
        cls.def("setIdFactory", &SimpleTable::setIdFactory, "idFactory"_a);
        cls.def("clone", &SimpleTable::clone);
        cls.def("makeRecord", &SimpleTable::makeRecord);
        cls.def("copyRecord",
                (std::shared_ptr<SimpleRecord>(SimpleTable::*)(BaseRecord const &)) & SimpleTable::copyRecord,
                "other"_a);
        cls.def("copyRecord",
                (std::shared_ptr<SimpleRecord>(SimpleTable::*)(BaseRecord const &, SchemaMapper const &)) &
                        SimpleTable::copyRecord,
                "other"_a, "mapper"_a);
    });
}

template <typename T>
void declareSchemaMapperOverloads(PySchemaMapper &cls, std::string const &suffix) {
    cls.def("getMapping", (Key<T>(SchemaMapper::*)(Key<T> const &) const) & SchemaMapper::getMapping);
    cls.def("isMapped", (bool (SchemaMapper::*)(Key<T> const &) const) & SchemaMapper::isMapped);
};

using PySchema = py::class_<Schema>;

using PySubSchema = py::class_<SubSchema>;

template <typename T>
using PyFieldBase = py::class_<FieldBase<T>>;

template <typename T>
using PyKeyBase = py::class_<KeyBase<T>>;

template <typename T>
using PyField = py::class_<Field<T>, FieldBase<T>>;

template <typename T>
using PyKey = py::class_<Key<T>, KeyBase<T>, FieldBase<T>>;

template <typename T>
using PySchemaItem = py::class_<SchemaItem<T>>;

// Specializations for FieldBase

template <typename T>
void declareFieldBaseSpecializations(PyFieldBase<T> &cls) {
    cls.def(py::init<>());
}

template <typename T>
void declareFieldBaseSpecializations(PyFieldBase<Array<T>> &cls) {
    cls.def(py::init<int>(), "size"_a = 0);
    cls.def("getSize", &FieldBase<Array<T>>::getSize);
    cls.def("isVariableLength", &FieldBase<Array<T>>::isVariableLength);
}

void declareFieldBaseSpecializations(PyFieldBase<std::string> &cls) {
    cls.def(py::init<int>(), "size"_a = -1);
    cls.def("getSize", &FieldBase<std::string>::getSize);
}

// Specializations for Field

template <typename T>
void declareFieldSpecializations(PyField<T> &cls) {
    cls.def(py::pickle(
            [](Field<T> const &self) {
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(self.getName(), self.getDoc(), self.getUnits());
            },
            [](py::tuple t) {
                int const NPARAMS = 3;
                if (t.size() != NPARAMS) {
                    std::ostringstream os;
                    os << "Invalid number of parameters (" << t.size() << ") when unpickling; expected "
                       << NPARAMS;
                    throw std::runtime_error(os.str());
                }
                return Field<T>(t[0].cast<std::string>(), t[1].cast<std::string>(), t[2].cast<std::string>());
            }));
}

// Field<Array<T>> and Field<std::string> have the same pickle implementation
template <typename T>
void _sequenceFieldSpecializations(PyField<T> &cls) {
    cls.def(py::pickle(
            [](Field<T> const &self) {
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(self.getName(), self.getDoc(), self.getUnits(), self.getSize());
            },
            [](py::tuple t) {
                int const NPARAMS = 4;
                if (t.size() != NPARAMS) {
                    std::ostringstream os;
                    os << "Invalid number of parameters (" << t.size() << ") when unpickling; expected "
                       << NPARAMS;
                    throw std::runtime_error(os.str());
                }
                return Field<T>(t[0].cast<std::string>(), t[1].cast<std::string>(), t[2].cast<std::string>(),
                                t[3].cast<int>());
            }));
}

template <typename T>
void declareFieldSpecializations(PyField<Array<T>> &cls) {
    _sequenceFieldSpecializations(cls);
}

void declareFieldSpecializations(PyField<std::string> &cls) { _sequenceFieldSpecializations(cls); }

// Specializations for KeyBase

template <typename T>
void declareKeyBaseSpecializations(PyKeyBase<T> &) {}

template <typename T>
void declareKeyBaseSpecializations(PyKeyBase<Array<T>> &cls) {
    cls.def("__getitem__", [](Key<Array<T>> const &self, py::object const &index) -> py::object {
        if (py::isinstance<py::slice>(index)) {
            py::slice slice(index);
            py::size_t start = 0, stop = 0, step = 0, length = 0;
            bool valid = slice.compute(self.getSize(), &start, &stop, &step, &length);
            if (!valid) throw py::error_already_set();
            if (step != 1) {
                throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                                  "Step for array Key indexing must be 1.");
            }
            return py::cast(self.slice(start, stop));
        } else {
            return py::cast(self[py::cast<int>(index)]);
        }
    });
    cls.def("slice", &KeyBase<Array<T>>::slice);
}

// Specializations for Key

template <typename T>
void declareKeyAccessors(PyKey<T> &cls) {
    cls.def("get", [](Key<T> const &self, BaseRecord &record) { return record.get(self); });
    cls.def("set", [](Key<T> const &self, BaseRecord &record, typename Key<T>::Value const &value) {
        record.set(self, value);
    });
}

template <typename U>
void declareKeyAccessors(PyKey<Array<U>> &cls) {
    auto getter = [](Key<Array<U>> const &self, BaseRecord &record) -> ndarray::Array<U, 1, 1> {
        return record[self];
    };
    auto setter = [](Key<Array<U>> const &self, BaseRecord &record, py::object const &value) {
        if (self.getSize() == 0) {
            // Variable-length array field: do a shallow copy, which requires a non-const
            // contiguous array.
            record.set(self, py::cast<ndarray::Array<U, 1, 1>>(value));
        } else {
            // Fixed-length array field: do a deep copy, which can work with a const
            // noncontiguous array.  But we need to check the size first, since the
            // penalty for getting that wrong is assert->abort.
            auto v = py::cast<ndarray::Array<U const, 1, 0>>(value);
            ndarray::ArrayRef<U, 1, 1> ref = record[self];
            if (v.size() != ref.size()) {
                throw LSST_EXCEPT(
                        pex::exceptions::LengthError,
                        (boost::format("Array sizes do not agree: %s != %s") % v.size() % ref.size()).str());
            }
            ref = v;
        }
        return;
    };
    cls.def("get", getter);
    cls.def("set", setter);
}

template <typename T>
void declareKeySpecializations(PyKey<T> &cls) {
    declareKeyAccessors(cls);
    cls.def_property_readonly("subfields", [](py::object const &) { return py::none(); });
    cls.def_property_readonly("subkeys", [](py::object const &) { return py::none(); });
    cls.def(py::pickle(
            [](Key<T> const &self) {
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(self.getOffset());
            },
            [](py::tuple t) {
                int const NPARAMS = 1;
                if (t.size() != NPARAMS) {
                    std::ostringstream os;
                    os << "Invalid number of parameters (" << t.size() << ") when unpickling; expected "
                       << NPARAMS;
                    throw std::runtime_error(os.str());
                }
                return detail::Access::makeKey<T>(t[0].cast<int>());
            }));
}

void declareKeySpecializations(PyKey<Flag> &cls) {
    declareKeyAccessors(cls);
    cls.def_property_readonly("subfields", [](py::object const &) { return py::none(); });
    cls.def_property_readonly("subkeys", [](py::object const &) { return py::none(); });
    cls.def("getBit", &Key<Flag>::getBit);
    cls.def(py::pickle(
            [](Key<Flag> const &self) {
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(self.getOffset(), self.getBit());
            },
            [](py::tuple t) {
                int const NPARAMS = 2;
                if (t.size() != NPARAMS) {
                    std::ostringstream os;
                    os << "Invalid number of parameters (" << t.size() << ") when unpickling; expected "
                       << NPARAMS;
                    throw std::runtime_error(os.str());
                }
                return detail::Access::makeKey(t[0].cast<int>(), t[1].cast<int>());
            }));
}

template <typename T>
void declareKeySpecializations(PyKey<Array<T>> &cls) {
    declareKeyAccessors(cls);
    cls.def_property_readonly("subfields", [](Key<Array<T>> const &self) -> py::object {
        py::list result;
        for (int i = 0; i < self.getSize(); ++i) {
            result.append(py::cast(i));
        }
        return py::tuple(result);
    });
    cls.def_property_readonly("subkeys", [](Key<Array<T>> const &self) -> py::object {
        py::list result;
        for (int i = 0; i < self.getSize(); ++i) {
            result.append(py::cast(self[i]));
        }
        return py::tuple(result);
    });
    cls.def(py::pickle(
            [](Key<Array<T>> const &self) {
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(self.getOffset(), self.getElementCount());
            },
            [](py::tuple t) {
                int const NPARAMS = 2;
                if (t.size() != NPARAMS) {
                    std::ostringstream os;
                    os << "Invalid number of parameters (" << t.size() << ") when unpickling; expected "
                       << NPARAMS;
                    throw std::runtime_error(os.str());
                }
                return detail::Access::makeKeyArray<T>(t[0].cast<int>(), t[1].cast<int>());
            }));
}

void declareKeySpecializations(PyKey<std::string> &cls) {
    declareKeyAccessors(cls);
    cls.def_property_readonly("subfields", [](py::object const &) { return py::none(); });
    cls.def_property_readonly("subkeys", [](py::object const &) { return py::none(); });
    cls.def(py::pickle(
            [](Key<std::string> const &self) {
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(self.getOffset(), self.getElementCount());
            },
            [](py::tuple t) {
                int const NPARAMS = 2;
                if (t.size() != NPARAMS) {
                    std::ostringstream os;
                    os << "Invalid number of parameters (" << t.size() << ") when unpickling; expected "
                       << NPARAMS;
                    throw std::runtime_error(os.str());
                }
                return detail::Access::makeKeyString(t[0].cast<int>(), t[1].cast<int>());
            }));
}

// Wrap all helper classes (FieldBase, KeyBase, Key, Field, SchemaItem) declarefor a Schema field type.
template <typename T>
void declareSchemaType(WrapperCollection &wrappers) {
    std::string suffix = FieldBase<T>::getTypeString();
    py::str pySuffix(suffix);

    py::object astropyUnit = py::module::import("astropy.units").attr("Unit");

    // FieldBase
    wrappers.wrapType(PyFieldBase<T>(wrappers.module, ("FieldBase" + suffix).c_str()),
                      [](auto &mod, auto &cls) {
                          cls.def_static("getTypeString", &FieldBase<T>::getTypeString);
                          declareFieldBaseSpecializations(cls);
                      });

    // KeyBase
    wrappers.wrapType(PyKeyBase<T>(wrappers.module, ("KeyBase" + suffix).c_str()), [](auto &mod, auto &cls) {
        cls.def_readonly_static("HAS_NAMED_SUBFIELDS", &KeyBase<T>::HAS_NAMED_SUBFIELDS);
        declareKeyBaseSpecializations(cls);
    });

    // Field
    wrappers.wrapType(PyField<T>(wrappers.module, ("Field" + suffix).c_str()), [pySuffix, astropyUnit](
                                                                                       auto &mod, auto &cls) {
        declareFieldSpecializations(cls);

        mod.attr("_Field")[pySuffix] = cls;

        cls.def(py::init([astropyUnit](  // capture by value to refcount in Python instead of dangle in C++
                                 std::string const &name, std::string const &doc, py::str const &units,
                                 py::object const &size, py::str const &parse_strict) {
                    astropyUnit(units, "parse_strict"_a = parse_strict);
                    std::string u = py::cast<std::string>(units);
                    if (size.is(py::none())) {
                        return new Field<T>(name, doc, u);
                    } else {
                        int s = py::cast<int>(size);
                        return new Field<T>(name, doc, u, s);
                    }
                }),
                "name"_a, "doc"_a = "", "units"_a = "", "size"_a = py::none(), "parse_strict"_a = "raise");
        cls.def("_addTo", [](Field<T> const &self, Schema &schema, bool doReplace) -> Key<T> {
            return schema.addField(self, doReplace);
        });
        cls.def("getName", &Field<T>::getName);
        cls.def("getDoc", &Field<T>::getDoc);
        cls.def("getUnits", &Field<T>::getUnits);
        cls.def("copyRenamed", &Field<T>::copyRenamed);
        utils::python::addOutputOp(cls, "__str__");
        utils::python::addOutputOp(cls, "__repr__");
    });

    // Key
    wrappers.wrapType(PyKey<T>(wrappers.module, ("Key" + suffix).c_str()), [pySuffix](auto &mod, auto &cls) {
        mod.attr("_Key")[pySuffix] = cls;
        cls.def(py::init<>());
        cls.def("__eq__", [](Key<T> const &self, Key<T> const &other) -> bool { return self == other; },
                py::is_operator());
        cls.def("__ne__", [](Key<T> const &self, Key<T> const &other) -> bool { return self != other; },
                py::is_operator());
        cls.def("isValid", &Key<T>::isValid);
        cls.def("getOffset", &Key<T>::getOffset);
        utils::python::addOutputOp(cls, "__str__");
        utils::python::addOutputOp(cls, "__repr__");
        // The Key methods below actually wrap templated methods on Schema and
        // SchemaMapper.  Rather than doing many-type overload resolution by
        // wrapping those methods directly, we use the visitor pattern by having
        // the wrappers for those methods delegate back to these non-templated
        // methods on the templated Key classes.
        cls.def("_findIn", [](Key<T> const &self, Schema const &schema) { return schema.find(self); });
        cls.def("_addMappingTo", [](Key<T> const &self, SchemaMapper &mapper, Field<T> const &field,
                                    bool doReplace) { return mapper.addMapping(self, field, doReplace); });
        cls.def("_addMappingTo", [](Key<T> const &self, SchemaMapper &mapper, std::string const &name,
                                    bool doReplace) { return mapper.addMapping(self, name, doReplace); });
        cls.def("_addMappingTo", [](Key<T> const &self, SchemaMapper &mapper, py::object const &,
                                    bool doReplace) { return mapper.addMapping(self, doReplace); });
        declareKeySpecializations(cls);
    });

    // SchemaItem
    wrappers.wrapType(PySchemaItem<T>(wrappers.module, ("SchemaItem" + suffix).c_str()),
                      [pySuffix](auto &mod, auto &cls) {
                          mod.attr("_SchemaItem")[pySuffix] = cls;
                          cls.def_readonly("key", &SchemaItem<T>::key);
                          cls.def_readonly("field", &SchemaItem<T>::field);
                          cls.def("getKey", [](SchemaItem<T> const &self) { return self.key; });
                          cls.def("getField", [](SchemaItem<T> const &self) { return self.field; });
                          cls.def("__getitem__", [](py::object const &self, int index) -> py::object {
                              if (index == 0) {
                                  return self.attr("key");
                              } else if (index == 1) {
                                  return self.attr("field");
                              }
                              // Have to raise IndexError not some LSST exception to get the
                              // right behavior when unpacking.
                              throw py::index_error("Index to SchemaItem must be 0 or 1.");
                          });
                          cls.def("__len__", [](py::object const &self) -> int { return 2; });
                          cls.def("__str__",
                                  [](py::object const &self) -> py::str { return py::str(py::tuple(self)); });
                          cls.def("__repr__", [](py::object const &self) -> py::str {
                              return py::str("SchemaItem(key={0.key}, field={0.field})").format(self);
                          });
                          cls.def(py::pickle(
                                  [](SchemaItem<T> const &self) {
                                      /* Return a tuple that fully encodes the state of the object */
                                      return py::make_tuple(self.key, self.field);
                                  },
                                  [](py::tuple t) {
                                      int const NPARAMS = 2;
                                      if (t.size() != NPARAMS) {
                                          std::ostringstream os;
                                          os << "Invalid number of parameters (" << t.size()
                                             << ") when unpickling; expected " << NPARAMS;
                                          throw std::runtime_error(os.str());
                                      }
                                      return SchemaItem<T>(t[0].cast<Key<T>>(), t[1].cast<Field<T>>());
                                  }));
                      });
}

// Helper class for Schema::find(name, func) that converts the result to Python.
// In C++14, this should be converted to a universal lambda.
struct MakePythonSchemaItem {
    template <typename T>
    void operator()(SchemaItem<T> const &item) {
        result = py::cast(item);
    }

    py::object result;
};

void declareSchema(WrapperCollection &wrappers) {
    wrappers.wrapType(PySchema(wrappers.module, "Schema"), [](auto &mod, auto &cls) {
        // wrap ComparisonFlags values as ints since we use them as bitflags,
        // not true enums
        cls.attr("EQUAL_KEYS") = py::cast(int(Schema::EQUAL_KEYS));
        cls.attr("EQUAL_NAMES") = py::cast(int(Schema::EQUAL_NAMES));
        cls.attr("EQUAL_DOCS") = py::cast(int(Schema::EQUAL_DOCS));
        cls.attr("EQUAL_UNITS") = py::cast(int(Schema::EQUAL_UNITS));
        cls.attr("EQUAL_FIELDS") = py::cast(int(Schema::EQUAL_FIELDS));
        cls.attr("EQUAL_ALIASES") = py::cast(int(Schema::EQUAL_ALIASES));
        cls.attr("IDENTICAL") = py::cast(int(Schema::IDENTICAL));

        cls.attr("VERSION") = py::cast(int(Schema::VERSION));

        cls.def(py::init<>());
        cls.def(py::init<Schema const &>());
        cls.def("__getitem__", [](Schema &self, std::string const &name) { return self[name]; });
        cls.def("__eq__", [](Schema const &self, Schema const &other) { return self == other; },
                py::is_operator());
        cls.def("__ne__", [](Schema const &self, Schema const &other) { return self != other; },
                py::is_operator());
        cls.def("getRecordSize", &Schema::getRecordSize);
        cls.def("getFieldCount", &Schema::getFieldCount);
        cls.def("getFlagFieldCount", &Schema::getFlagFieldCount);
        cls.def("getNonFlagFieldCount", &Schema::getNonFlagFieldCount);
        cls.def("find", [](py::object const &self, py::object const &key) -> py::object {
            try {
                if (py::isinstance<py::str>(key) || py::isinstance<py::bytes>(key)) {
                    Schema const &s = py::cast<Schema const &>(self);
                    std::string name = py::cast<std::string>(key);
                    MakePythonSchemaItem func;
                    s.findAndApply(name, func);
                    return func.result;
                }
                return key.attr("_findIn")(self);
            } catch (pex::exceptions::NotFoundError &err) {
                // Avoid API change by re-throwing as KeyError.
                PyErr_SetString(PyExc_KeyError, err.what());
                throw py::error_already_set();
            }
        });
        cls.def("getNames", &Schema::getNames, "topOnly"_a = false);
        cls.def("getAliasMap", &Schema::getAliasMap);
        cls.def("setAliasMap", &Schema::setAliasMap, "aliases"_a);
        cls.def("disconnectAliases", &Schema::disconnectAliases);
        cls.def("forEach", [](Schema &self, py::object &obj) { self.forEach(obj); });
        cls.def("compare", &Schema::compare, "other"_a, "flags"_a = int(Schema::EQUAL_KEYS));
        cls.def("contains", (int (Schema::*)(Schema const &, int) const) & Schema::contains, "other"_a,
                "flags"_a = int(Schema::EQUAL_KEYS));
        cls.def("__contains__", [](py::object const &self, py::object const &key) {
            try {
                self.attr("find")(key);
            } catch (py::error_already_set &err) {
                err.restore();
                PyErr_Clear();
                return false;
            }
            return true;
        });
        cls.def_static("readFits", (Schema(*)(std::string const &, int)) & Schema::readFits, "filename"_a,
                       "hdu"_a = fits::DEFAULT_HDU);
        cls.def_static("readFits", (Schema(*)(fits::MemFileManager &, int)) & Schema::readFits, "manager"_a,
                       "hdu"_a = fits::DEFAULT_HDU);

        cls.def("join",
                (std::string(Schema::*)(std::string const &, std::string const &) const) & Schema::join,
                "a"_a, "b"_a);
        cls.def("join",
                (std::string(Schema::*)(std::string const &, std::string const &, std::string const &)
                         const) &
                        Schema::join,
                "a"_a, "b"_a, "c"_a);
        cls.def("join",
                (std::string(Schema::*)(std::string const &, std::string const &, std::string const &,
                                        std::string const &) const) &
                        Schema::join,
                "a"_a, "b"_a, "c"_a, "d"_a);
        utils::python::addOutputOp(cls, "__str__");
        utils::python::addOutputOp(cls, "__repr__");
    });
}

void declareSubSchema(WrapperCollection &wrappers) {
    wrappers.wrapType(PySubSchema(wrappers.module, "SubSchema"), [](auto &mod, auto &cls) {
        cls.def("getNames", &SubSchema::getNames, "topOnly"_a = false);
        cls.def("getPrefix", &SubSchema::getPrefix);
        cls.def("asKey", [](SubSchema const &self) -> py::object {
            MakePythonSchemaItem func;
            self.apply(func);
            return func.result.attr("key");
        });
        cls.def("asField", [](SubSchema const &self) -> py::object {
            MakePythonSchemaItem func;
            self.apply(func);
            return func.result.attr("field");
        });
        cls.def("find", [](SubSchema const &self, std::string const &name) -> py::object {
            MakePythonSchemaItem func;
            self.findAndApply(name, func);
            return func.result;
        });
        cls.def("__getitem__", [](SubSchema &self, std::string const &name) { return self[name]; });
    });
}

/// @internal Declare match code templated on two types of catalog
template <typename Catalog1, typename Catalog2>
void declareMatch2(WrapperCollection &wrappers, std::string const &prefix) {
    typedef typename Catalog1::Record Record1;
    typedef typename Catalog2::Record Record2;
    typedef std::vector<Match<typename Catalog1::Record, typename Catalog2::Record>> MatchList;

    using Class = Match<Record1, Record2>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;
    wrappers.wrapType(PyClass(wrappers.module, (prefix + "Match").c_str()), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<std::shared_ptr<Record1> const &, std::shared_ptr<Record2> const &, double>(),
                "first"_a, "second"_a, "distance"_a);

        // struct fields
        cls.def_readwrite("first", &Match<Record1, Record2>::first);
        cls.def_readwrite("second", &Match<Record1, Record2>::second);
        cls.def_readwrite("distance", &Match<Record1, Record2>::distance);
    });

    // Free Functions
    wrappers.wrap([](auto &mod) {
        mod.def("unpackMatches", &unpackMatches<Catalog1, Catalog2>, "matches"_a, "cat1"_a, "cat2"_a);

        mod.def("matchRaDec",
                (MatchList(*)(Catalog1 const &, Catalog2 const &, lsst::geom::Angle,
                              MatchControl const &))matchRaDec<Catalog1, Catalog2>,
                "cat1"_a, "cat2"_a, "radius"_a, "mc"_a = MatchControl());
    });
};

/// @internal Declare match code templated on one type of catalog
template <typename Catalog>
void declareMatch1(WrapperCollection &wrappers) {
    typedef std::vector<Match<typename Catalog::Record, typename Catalog::Record>> MatchList;
    wrappers.wrap([](auto &mod) {
        mod.def("matchRaDec",
                (MatchList(*)(Catalog const &, lsst::geom::Angle, MatchControl const &))matchRaDec<Catalog>,
                "cat"_a, "radius"_a, "mc"_a = MatchControl());
    });
}

using PyExposureRecord = py::class_<ExposureRecord, std::shared_ptr<ExposureRecord>, BaseRecord>;
using PyExposureTable = py::class_<ExposureTable, std::shared_ptr<ExposureTable>, BaseTable>;
using PyExposureCatalog =
        py::class_<ExposureCatalogT<ExposureRecord>, std::shared_ptr<ExposureCatalogT<ExposureRecord>>,
                   SortedCatalogT<ExposureRecord>>;

PyExposureRecord declareExposureRecord(WrapperCollection &wrappers) {
    return wrappers.wrapType(PyExposureRecord(wrappers.module, "ExposureRecord"), [](auto &mod, auto &cls) {
        cls.def("getId", &ExposureRecord::getId);
        cls.def("setId", &ExposureRecord::setId, "id"_a);
        cls.def("getBBox", &ExposureRecord::getBBox);
        cls.def("setBBox", &ExposureRecord::setBBox, "bbox"_a);
        cls.def("getTable", &ExposureRecord::getTable);
        cls.def_property_readonly("table", &ExposureRecord::getTable);
        cls.def("contains",
                (bool (ExposureRecord::*)(lsst::geom::SpherePoint const &, bool) const) &
                        ExposureRecord::contains,
                "coord"_a, "includeValidPolygon"_a = false);
        cls.def("contains",
                (bool (ExposureRecord::*)(lsst::geom::Point2D const &, geom::SkyWcs const &, bool) const) &
                        ExposureRecord::contains,
                "point"_a, "wcs"_a, "includeValidPolygon"_a = false);
        cls.def("getWcs", &ExposureRecord::getWcs);
        cls.def("setWcs", &ExposureRecord::setWcs, "wcs"_a);
        cls.def("getPsf", &ExposureRecord::getPsf);
        cls.def("setPsf", &ExposureRecord::setPsf, "psf"_a);

        cls.def("getPhotoCalib", &ExposureRecord::getPhotoCalib);
        cls.def("setPhotoCalib", &ExposureRecord::setPhotoCalib, "photoCalib"_a);
        cls.def("getApCorrMap", &ExposureRecord::getApCorrMap);
        cls.def("setApCorrMap", &ExposureRecord::setApCorrMap, "apCorrMap"_a);
        cls.def("getValidPolygon", &ExposureRecord::getValidPolygon);

        // Workaround for DM-10289.
        cls.def("setValidPolygon",
                [](ExposureRecord &self, py::object polygon) {
                    if (polygon.is(py::none())) {
                        self.setValidPolygon(nullptr);
                    } else {
                        self.setValidPolygon(py::cast<std::shared_ptr<afw::geom::polygon::Polygon>>(polygon));
                    }
                },
                "polygon"_a);

        cls.def("getVisitInfo", &ExposureRecord::getVisitInfo);
        cls.def("setVisitInfo", &ExposureRecord::setVisitInfo, "visitInfo"_a);
        cls.def("getTransmissionCurve", &ExposureRecord::getTransmissionCurve);
        cls.def("setTransmissionCurve", &ExposureRecord::setTransmissionCurve, "transmissionCurve"_a);
        cls.def("getDetector", &ExposureRecord::getDetector);
        cls.def("setDetector", &ExposureRecord::setDetector, "detector"_a);
    });
}

PyExposureTable declareExposureTable(WrapperCollection &wrappers) {
    return wrappers.wrapType(PyExposureTable(wrappers.module, "ExposureTable"), [](auto &mod, auto &cls) {
        cls.def_static("make", &ExposureTable::make);
        cls.def_static("makeMinimalSchema", &ExposureTable::makeMinimalSchema);
        cls.def_static("checkSchema", &ExposureTable::checkSchema, "schema"_a);

        cls.def_static("getIdKey", &ExposureTable::getIdKey);
        cls.def_static("getBBoxMinKey", &ExposureTable::getBBoxMinKey);
        cls.def_static("getBBoxMaxKey", &ExposureTable::getBBoxMaxKey);

        cls.def("clone", &ExposureTable::clone);
        cls.def("makeRecord", &ExposureTable::makeRecord);
        cls.def("copyRecord", (std::shared_ptr<ExposureRecord>(ExposureTable::*)(BaseRecord const &)) &
                                      ExposureTable::copyRecord);
        cls.def("copyRecord", (std::shared_ptr<ExposureRecord>(ExposureTable::*)(BaseRecord const &,
                                                                                 SchemaMapper const &)) &
                                      ExposureTable::copyRecord);
    });
}

PyExposureCatalog declareExposureCatalog(WrapperCollection &wrappers) {
    using Catalog = ExposureCatalogT<ExposureRecord>;
    table::python::declareSortedCatalog<ExposureRecord>(wrappers, "Exposure", true);

    // We need py::dynamic_attr() in class definition to support our Python-side caching
    // of the associated ColumnView.
    return wrappers.wrapType(
            PyExposureCatalog(wrappers.module, "ExposureCatalog", py::dynamic_attr()),
            [](auto &mod, auto &cls) {
                cls.def(py::init<Schema const &>(), "schema"_a);
                cls.def(py::init<std::shared_ptr<ExposureTable> const &>(),
                        "table"_a = std::shared_ptr<ExposureTable>());
                cls.def(py::init<Catalog const &>(), "other"_a);
                // Constructor taking C++ iterators not wrapped; we recommend .extend() (defined in pure
                // Python) instead.
                cls.def_static("readFits", (Catalog(*)(std::string const &, int, int)) & Catalog::readFits,
                               "filename"_a, "hdu"_a = fits::DEFAULT_HDU, "flags"_a = 0);
                cls.def_static("readFits", (Catalog(*)(fits::MemFileManager &, int, int)) & Catalog::readFits,
                               "manager"_a, "hdu"_a = fits::DEFAULT_HDU, "flags"_a = 0);
                // readFits taking Fits objects not wrapped, because Fits objects are not wrapped.

                cls.def("subset",
                        (Catalog(Catalog::*)(ndarray::Array<bool const, 1> const &) const) & Catalog::subset,
                        "mask"_a);
                cls.def("subset",
                        (Catalog(Catalog::*)(std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t) const) &
                                Catalog::subset,
                        "startd"_a, "stopd"_a, "step"_a);
                cls.def("subsetContaining",
                        (Catalog(Catalog::*)(lsst::geom::SpherePoint const &, bool) const) &
                                Catalog::subsetContaining,
                        "coord"_a, "includeValidPolygon"_a = false);
                cls.def("subsetContaining",
                        (Catalog(Catalog::*)(lsst::geom::Point2D const &, geom::SkyWcs const &, bool) const) &
                                Catalog::subsetContaining,
                        "point"_a, "wcs"_a, "includeValidPolygon"_a = false);
            });
};

// We don't expose base classes (e.g. FunctorKey) to Python, since they're just used to
// define a CRTP interface in C++ and in Python that's just duck-typing.

template <typename T>
using PyPointKey = py::class_<PointKey<T>, std::shared_ptr<PointKey<T>>>;

template <typename Box>
using PyBoxKey = py::class_<BoxKey<Box>, std::shared_ptr<BoxKey<Box>>>;

using PyCoordKey = py::class_<CoordKey, std::shared_ptr<CoordKey>>;

using PyQuadrupoleKey = py::class_<QuadrupoleKey, std::shared_ptr<QuadrupoleKey>>;

using PyEllipseKey = py::class_<EllipseKey, std::shared_ptr<EllipseKey>>;

template <typename T, int N>
using PyCovarianceMatrixKey =
        py::class_<CovarianceMatrixKey<T, N>, std::shared_ptr<CovarianceMatrixKey<T, N>>>;

template <typename T>
static void declarePointKey(WrapperCollection &wrappers, std::string const &suffix) {
    wrappers.wrapType(
            PyPointKey<T>(wrappers.module, ("Point" + suffix + "Key").c_str()), [](auto &mod, auto &cls) {
                cls.def(py::init<>());
                cls.def(py::init<Key<T> const &, Key<T> const &>(), "x"_a, "y"_a);
                cls.def(py::init<SubSchema const &>());
                cls.def("__eq__", &PointKey<T>::operator==, py::is_operator());
                cls.def("__ne__", &PointKey<T>::operator!=, py::is_operator());
                cls.def("getX", &PointKey<T>::getX);
                cls.def("getY", &PointKey<T>::getY);
                cls.def("isValid", &PointKey<T>::isValid);
                cls.def_static("addFields", &PointKey<T>::addFields, "schema"_a, "name"_a, "doc"_a, "unit"_a);
                cls.def("set", [](PointKey<T> &self, BaseRecord &record,
                                  lsst::geom::Point<T, 2> const &value) { return self.set(record, value); });
                cls.def("get", &PointKey<T>::get);
            });
};

template <typename Box>
static void declareBoxKey(WrapperCollection &wrappers, std::string const &suffix) {
    wrappers.wrapType(
            PyBoxKey<Box>(wrappers.module, ("Box" + suffix + "Key").c_str()), [](auto &mod, auto &cls) {
                using Element = typename Box::Element;
                cls.def(py::init<>());
                cls.def(py::init<PointKey<Element> const &, PointKey<Element> const &>(), "min"_a, "max"_a);
                cls.def(py::init<SubSchema const &>());
                cls.def("__eq__", &BoxKey<Box>::operator==, py::is_operator());
                cls.def("__ne__", &BoxKey<Box>::operator!=, py::is_operator());
                cls.def("getMin", &BoxKey<Box>::getMin);
                cls.def("getMax", &BoxKey<Box>::getMax);
                cls.def("isValid", &BoxKey<Box>::isValid);
                cls.def_static("addFields", &BoxKey<Box>::addFields, "schema"_a, "name"_a, "doc"_a, "unit"_a);
                cls.def("set", &BoxKey<Box>::set);
                cls.def("get", &BoxKey<Box>::get);

            });
};

static void declareCoordKey(WrapperCollection &wrappers) {
    wrappers.wrapType(PyCoordKey(wrappers.module, "CoordKey"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<Key<lsst::geom::Angle>, Key<lsst::geom::Angle>>(), "ra"_a, "dec"_a);
        cls.def(py::init<SubSchema const &>());
        cls.def("__eq__", &CoordKey::operator==, py::is_operator());
        cls.def("__ne__", &CoordKey::operator!=, py::is_operator());
        cls.def_static("addFields", &CoordKey::addFields, "schema"_a, "name"_a, "doc"_a);
        cls.def("getRa", &CoordKey::getRa);
        cls.def("getDec", &CoordKey::getDec);
        cls.def("isValid", &CoordKey::isValid);
        cls.def("get", [](CoordKey &self, BaseRecord const &record) { return self.get(record); });
        cls.def("set", &CoordKey::set);
    });
}

static void declareQuadrupoleKey(WrapperCollection &wrappers) {
    wrappers.wrapType(PyQuadrupoleKey(wrappers.module, "QuadrupoleKey"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<Key<double> const &, Key<double> const &, Key<double> const &>(), "ixx"_a, "iyy"_a,
                "ixy"_a);
        cls.def(py::init<SubSchema const &>());
        cls.def("__eq__", &QuadrupoleKey::operator==, py::is_operator());
        cls.def("__nq__", &QuadrupoleKey::operator!=, py::is_operator());
        cls.def_static("addFields", &QuadrupoleKey::addFields, "schema"_a, "name"_a, "doc"_a,
                       "coordType"_a = CoordinateType::PIXEL);
        cls.def("getIxx", &QuadrupoleKey::getIxx);
        cls.def("getIyy", &QuadrupoleKey::getIyy);
        cls.def("getIxy", &QuadrupoleKey::getIxy);
        cls.def("isValid", &QuadrupoleKey::isValid);
        cls.def("set", &QuadrupoleKey::set);
        cls.def("get", &QuadrupoleKey::get);
    });
}

static void declareEllipseKey(WrapperCollection &wrappers) {
    wrappers.wrapType(PyEllipseKey(wrappers.module, "EllipseKey"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<QuadrupoleKey const &, PointKey<double> const &>(), "qKey"_a, "pKey"_a);
        cls.def(py::init<SubSchema const &>());
        cls.def("__eq__", &EllipseKey::operator==, py::is_operator());
        cls.def("__nq__", &EllipseKey::operator!=, py::is_operator());
        cls.def_static("addFields", &EllipseKey::addFields, "schema"_a, "name"_a, "doc"_a, "unit"_a);
        cls.def("get", &EllipseKey::get);
        cls.def("set", &EllipseKey::set);
        cls.def("isValid", &EllipseKey::isValid);
        cls.def("getCore", &EllipseKey::getCore);
        cls.def("getCenter", &EllipseKey::getCenter);
    });
}

template <typename T, int N>
static void declareCovarianceMatrixKey(WrapperCollection &wrappers, const ::std::string &suffix) {
    wrappers.wrapType(
            PyCovarianceMatrixKey<T, N>(wrappers.module, ("CovarianceMatrix" + suffix + "Key").c_str()),
            [](auto &mod, auto &cls) {
                typedef std::vector<Key<T>> ErrKeyArray;
                typedef std::vector<Key<T>> CovarianceKeyArray;
                typedef std::vector<std::string> NameArray;

                cls.def(py::init<>());
                // Ordering of the next two ctor declaration matters, as a workaround for DM-8580.
                cls.def(py::init<SubSchema const &, NameArray const &>());
                cls.def(py::init<ErrKeyArray const &, CovarianceKeyArray const &>(), "err"_a,
                        "cov"_a = CovarianceKeyArray());
                cls.def("__eq__", &CovarianceMatrixKey<T, N>::operator==, py::is_operator());
                cls.def("__ne__", &CovarianceMatrixKey<T, N>::operator!=, py::is_operator());
                cls.def_static("addFields",
                               (CovarianceMatrixKey<T, N>(*)(Schema &, std::string const &, NameArray const &,
                                                             std::string const &, bool)) &
                                       CovarianceMatrixKey<T, N>::addFields,
                               "schema"_a, "prefix"_a, "names"_a, "unit"_a, "diagonalOnly"_a = false);
                cls.def_static("addFields",
                               (CovarianceMatrixKey<T, N>(*)(Schema &, std::string const &, NameArray const &,
                                                             NameArray const &, bool)) &
                                       CovarianceMatrixKey<T, N>::addFields,
                               "schema"_a, "prefix"_a, "names"_a, "units"_a, "diagonalOnly"_a = false);
                cls.def("set", [](CovarianceMatrixKey<T, N> &cov, BaseRecord &record,
                                  Eigen::Matrix<T, N, N> const &value) { return cov.set(record, value); });
                cls.def("get", [](CovarianceMatrixKey<T, N> &cov, BaseRecord const &record) {
                    return cov.get(record);
                });
                cls.def("isValid", &CovarianceMatrixKey<T, N>::isValid);
                cls.def("setElement", &CovarianceMatrixKey<T, N>::setElement);
                cls.def("getElement", &CovarianceMatrixKey<T, N>::getElement);
            });
}


// We don't expose base classes (e.g. FunctorKey) to Python, since they're just used to
// define a CRTP interface in C++ and in Python that's just duck-typing.

template <typename T>
using PyArrayKey = py::class_<ArrayKey<T>, std::shared_ptr<ArrayKey<T>>>;

template <typename T>
void declareArrayKey(WrapperCollection &wrappers, std::string const &suffix) {
    wrappers.wrapType(
            PyArrayKey<T>(wrappers.module, ("Array" + suffix + "Key").c_str()), [](auto &mod, auto &cls) {
                cls.def(py::init<>());
                cls.def(py::init<Key<Array<T>> const &>());
                cls.def(py::init<SubSchema const &>());
                cls.def(py::init<std::vector<Key<T>> const &>());

                cls.def_static("addFields",
                               (ArrayKey<T>(*)(Schema &, std::string const &, std::string const &,
                                               std::string const &, std::vector<T> const &)) &
                                       ArrayKey<T>::addFields,
                               "schema"_a, "name"_a, "doc"_a, "unit"_a, "docData"_a);
                cls.def_static("addFields",
                               (ArrayKey<T>(*)(Schema &, std::string const &, std::string const &,
                                               std::string const &, int size)) &
                                       ArrayKey<T>::addFields,
                               "schema"_a, "name"_a, "doc"_a, "unit"_a, "size"_a);
                cls.def("get", &ArrayKey<T>::get);
                cls.def("set", &ArrayKey<T>::set);
                cls.def("isValid", &ArrayKey<T>::isValid);
                cls.def("__eq__", &ArrayKey<T>::operator==, py::is_operator());
                cls.def("__ne__", &ArrayKey<T>::operator!=, py::is_operator());
                cls.def("__getitem__",
                        // Special implementation of __getitem__ to support ints and slices
                        [](ArrayKey<T> const &self, py::object const &index) -> py::object {
                            if (py::isinstance<py::slice>(index)) {
                                py::slice slice(index);
                                py::size_t start = 0, stop = 0, step = 0, length = 0;
                                bool valid = slice.compute(self.getSize(), &start, &stop, &step, &length);
                                if (!valid) throw py::error_already_set();
                                if (step != 1) {
                                    throw py::index_error("Step for ArrayKey must be 1.");
                                }
                                return py::cast(self.slice(start, stop));
                            } else {
                                std::size_t n = utils::python::cppIndex(self.getSize(),
                                                                        py::cast<std::ptrdiff_t>(index));
                                return py::cast(self[n]);
                            }
                        });
                cls.def("getSize", &ArrayKey<T>::getSize);
                cls.def("slice", &ArrayKey<T>::slice);
            });
};

using PyBaseRecord = py::class_<BaseRecord, std::shared_ptr<BaseRecord>>;
using PyBaseTable = py::class_<BaseTable, std::shared_ptr<BaseTable>>;

template <typename T>
void declareBaseRecordOverloads(PyBaseRecord &cls, std::string const &suffix) {
    typedef typename Field<T>::Value (BaseRecord::*Getter)(Key<T> const &) const;
    typedef void (BaseRecord::*Setter)(Key<T> const &, typename Field<T>::Value const &);
    cls.def(("get" + suffix).c_str(), (Getter)&BaseRecord::get);
    cls.def(("set" + suffix).c_str(), (Setter)&BaseRecord::set);
}

template <typename T>
void declareBaseRecordArrayOverloads(PyBaseRecord &cls, std::string const &suffix) {
    auto getter = [](BaseRecord &self, Key<Array<T>> const &key) -> ndarray::Array<T, 1, 1> {
        return self[key];
    };
    auto setter = [](BaseRecord &self, Key<Array<T>> const &key, py::object const &value) {
        if (key.getSize() == 0) {
            // Variable-length array field: do a shallow copy, which requires a non-const
            // contiguous array.
            self.set(key, py::cast<ndarray::Array<T, 1, 1>>(value));
        } else {
            // Fixed-length array field: do a deep copy, which can work with a const
            // noncontiguous array.  But we need to check the size first, since the
            // penalty for getting that wrong is assert->abort.
            auto v = py::cast<ndarray::Array<T const, 1, 0>>(value);
            ndarray::ArrayRef<T, 1, 1> ref = self[key];
            if (v.size() != ref.size()) {
                throw LSST_EXCEPT(
                        pex::exceptions::LengthError,
                        (boost::format("Array sizes do not agree: %s != %s") % v.size() % ref.size()).str());
            }
            ref = v;
        }
        return;
    };
    cls.def(("get" + suffix).c_str(), getter);
    cls.def(("set" + suffix).c_str(), setter);
}

PyBaseRecord declareBaseRecord(WrapperCollection &wrappers) {
    return wrappers.wrapType(PyBaseRecord(wrappers.module, "BaseRecord"), [](auto &mod, auto &cls) {
        utils::python::addSharedPtrEquality<BaseRecord>(cls);
        cls.def("assign", (void (BaseRecord::*)(BaseRecord const &)) & BaseRecord::assign);
        cls.def("assign",
                (void (BaseRecord::*)(BaseRecord const &, SchemaMapper const &)) & BaseRecord::assign);
        cls.def("getSchema", &BaseRecord::getSchema);
        cls.def("getTable", &BaseRecord::getTable);
        cls.def_property_readonly("schema", &BaseRecord::getSchema);
        cls.def_property_readonly("table", &BaseRecord::getTable);

        declareBaseRecordOverloads<double>(cls, "D");
        declareBaseRecordOverloads<float>(cls, "F");
        declareBaseRecordOverloads<lsst::afw::table::Flag>(cls, "Flag");
        declareBaseRecordOverloads<std::uint8_t>(cls, "B");
        declareBaseRecordOverloads<std::uint16_t>(cls, "U");
        declareBaseRecordOverloads<std::int32_t>(cls, "I");
        declareBaseRecordOverloads<std::int64_t>(cls, "L");
        declareBaseRecordOverloads<std::string>(cls, "String");
        declareBaseRecordOverloads<lsst::geom::Angle>(cls, "Angle");
        declareBaseRecordArrayOverloads<std::uint8_t>(cls, "ArrayB");
        declareBaseRecordArrayOverloads<std::uint16_t>(cls, "ArrayU");
        declareBaseRecordArrayOverloads<int>(cls, "ArrayI");
        declareBaseRecordArrayOverloads<float>(cls, "ArrayF");
        declareBaseRecordArrayOverloads<double>(cls, "ArrayD");
        utils::python::addOutputOp(cls, "__str__");  // __repr__ is defined in baseContinued.py

        // These are master getters and setters that can take either strings, Keys, or
        // FunctorKeys, and dispatch to key.get.
        auto getter = [](py::object const &self, py::object key) -> py::object {
            py::object schema = self.attr("schema");
            if (py::isinstance<py::str>(key) || py::isinstance<py::bytes>(key)) {
                key = schema.attr("find")(key).attr("key");
            }
            return key.attr("get")(self);
        };
        auto setter = [](py::object const &self, py::object key, py::object const &value) -> void {
            py::object schema = self.attr("schema");
            if (py::isinstance<py::str>(key) || py::isinstance<py::bytes>(key)) {
                key = schema.attr("find")(key).attr("key");
            }
            key.attr("set")(self, value);
        };

        // The distinction between get/set and operator[] is meaningful in C++, because "record[k] = v"
        // operates by returning an object that can be assigned to.
        // But there's no meaningful difference between get/set and __getitem__/__setitem__.
        cls.def("get", getter);
        cls.def("__getitem__", getter);
        cls.def("set", setter);
        cls.def("__setitem__", setter);
    });
}

PyBaseTable declareBaseTable(WrapperCollection &wrappers) {
    return wrappers.wrapType(PyBaseTable(wrappers.module, "BaseTable"), [](auto &mod, auto &cls) {
        utils::python::addSharedPtrEquality<BaseTable>(cls);
        cls.def_static("make", &BaseTable::make);
        cls.def("getMetadata", &BaseTable::getMetadata);
        cls.def("setMetadata", &BaseTable::setMetadata, "metadata"_a);
        cls.def("popMetadata", &BaseTable::popMetadata);
        cls.def("makeRecord", &BaseTable::makeRecord);
        cls.def("copyRecord",
                (std::shared_ptr<BaseRecord>(BaseTable::*)(BaseRecord const &)) & BaseTable::copyRecord);
        cls.def("copyRecord",
                (std::shared_ptr<BaseRecord>(BaseTable::*)(BaseRecord const &, SchemaMapper const &)) &
                        BaseTable::copyRecord);
        cls.def("getSchema", &BaseTable::getSchema);
        cls.def_property_readonly("schema", &BaseTable::getSchema);
        cls.def("getBufferSize", &BaseTable::getBufferSize);
        cls.def("clone", &BaseTable::clone);
        cls.def("preallocate", &BaseTable::preallocate);
    });
}

using PyBaseColumnView = py::class_<BaseColumnView, std::shared_ptr<BaseColumnView>>;

using PyBitsColumn = py::class_<BitsColumn, std::shared_ptr<BitsColumn>>;

template <typename T, typename PyClass>
static void declareBaseColumnViewOverloads(PyClass &cls) {
    cls.def("_basicget", [](BaseColumnView & self, Key<T> const &key) -> typename ndarray::Array<T, 1> const {
        return self[key];
    });
};

template <typename U, typename PyClass>
static void declareBaseColumnViewArrayOverloads(PyClass &cls) {
    cls.def("_basicget",
            [](BaseColumnView & self, Key<lsst::afw::table::Array<U>> const &key) ->
            typename ndarray::Array<U, 2, 1> const { return self[key]; });
};

template <typename PyClass>
static void declareBaseColumnViewFlagOverloads(PyClass &cls) {
    cls.def("_basicget",
            [](BaseColumnView &self, Key<Flag> const &key) -> ndarray::Array<bool const, 1, 1> const {
                return ndarray::copy(self[key]);
            });
};

static void declareBaseColumnView(WrapperCollection &wrappers) {
    // We can't call this "BaseColumnView" because that's the typedef for "ColumnViewT<BaseRecord>".
    // This is just a mostly-invisible implementation base class, so we use the same naming convention
    // we use for those.
    wrappers.wrapType(PyBaseColumnView(wrappers.module, "_BaseColumnViewBase"), [](auto &mod, auto &cls) {
        cls.def("getTable", &BaseColumnView::getTable);
        cls.def_property_readonly("table", &BaseColumnView::getTable);
        cls.def("getSchema", &BaseColumnView::getSchema);
        cls.def_property_readonly("schema", &BaseColumnView::getSchema);
        // _getBits supports a Python version of getBits that accepts None and field names as keys
        cls.def("_getBits", &BaseColumnView::getBits);
        cls.def("getAllBits", &BaseColumnView::getAllBits);
        declareBaseColumnViewOverloads<std::uint8_t>(cls);
        declareBaseColumnViewOverloads<std::uint16_t>(cls);
        declareBaseColumnViewOverloads<std::int32_t>(cls);
        declareBaseColumnViewOverloads<std::int64_t>(cls);
        declareBaseColumnViewOverloads<float>(cls);
        declareBaseColumnViewOverloads<double>(cls);
        declareBaseColumnViewFlagOverloads(cls);
        // std::string columns are not supported, because numpy string arrays
        // do not have the same memory model as ours.
        declareBaseColumnViewArrayOverloads<std::uint8_t>(cls);
        declareBaseColumnViewArrayOverloads<std::uint16_t>(cls);
        declareBaseColumnViewArrayOverloads<int>(cls);
        declareBaseColumnViewArrayOverloads<float>(cls);
        declareBaseColumnViewArrayOverloads<double>(cls);
        // lsst::geom::Angle requires custom wrappers, because ndarray doesn't
        // recognize it natively; we just return a double view
        // (e.g. radians).
        using AngleArray = ndarray::Array<lsst::geom::Angle, 1>;
        using DoubleArray = ndarray::Array<double, 1>;
        cls.def("_basicget", [](BaseColumnView &self, Key<lsst::geom::Angle> const &key) -> DoubleArray {
            ndarray::Array<lsst::geom::Angle, 1, 0> a = self[key];
            return ndarray::detail::ArrayAccess<DoubleArray>::construct(
                    reinterpret_cast<double *>(a.getData()),
                    ndarray::detail::ArrayAccess<AngleArray>::getCore(a));
        });
    });
}

static void declareBitsColumn(WrapperCollection &wrappers) {
    wrappers.wrapType(PyBitsColumn(wrappers.module, "BitsColumn"), [](auto &mod, auto &cls) {
        cls.def("getArray", &BitsColumn::getArray);
        cls.def_property_readonly("array", &BitsColumn::getArray);
        cls.def("getBit", (BitsColumn::IntT(BitsColumn::*)(Key<Flag> const &) const) & BitsColumn::getBit,
                "key"_a);
        cls.def("getBit", (BitsColumn::IntT(BitsColumn::*)(std::string const &) const) & BitsColumn::getBit,
                "name"_a);
        cls.def("getMask", (BitsColumn::IntT(BitsColumn::*)(Key<Flag> const &) const) & BitsColumn::getMask,
                "key"_a);
        cls.def("getMask", (BitsColumn::IntT(BitsColumn::*)(std::string const &) const) & BitsColumn::getMask,
                "name"_a);
    });
}

void wrapWcsUtils(WrapperCollection &wrappers) {
    declareUpdateRefCentroids<std::vector<std::shared_ptr<lsst::afw::table::SimpleRecord>>>(wrappers);
    declareUpdateRefCentroids<lsst::afw::table::SimpleCatalog>(wrappers);

    declareUpdateSourceCoords<std::vector<std::shared_ptr<lsst::afw::table::SourceRecord>>>(wrappers);
    declareUpdateSourceCoords<lsst::afw::table::SourceCatalog>(wrappers);
}

void wrapSlots(WrapperCollection &wrappers) {

    declareSlotDefinition(wrappers);
    declareSlotDefinitionSubclass<FluxSlotDefinition>(wrappers, "Flux");
    declareSlotDefinitionSubclass<CentroidSlotDefinition>(wrappers, "Centroid");
    declareSlotDefinitionSubclass<ShapeSlotDefinition>(wrappers, "Shape");
}

void wrapSource(WrapperCollection &wrappers) {

    // SourceFitsFlags enum values are used as integer masks, so wrap as attributes instead of an enum
    // static_cast is required to avoid an import error (py::cast and py::int_ do not work by themselves
    // and are not required with the static_cast)
    auto &mod = wrappers.module;
    mod.attr("SOURCE_IO_NO_FOOTPRINTS") = static_cast<int>(SourceFitsFlags::SOURCE_IO_NO_FOOTPRINTS);
    mod.attr("SOURCE_IO_NO_HEAVY_FOOTPRINTS") =
            static_cast<int>(SourceFitsFlags::SOURCE_IO_NO_HEAVY_FOOTPRINTS);

    auto clsSourceRecord = declareSourceRecord(wrappers);
    auto clsSourceTable = declareSourceTable(wrappers);
    auto clsSourceColumnView = declareSourceColumnView(wrappers);
    auto clsSourceCatalog = table::python::declareSortedCatalog<SourceRecord>(wrappers, "Source");

    clsSourceRecord.attr("Table") = clsSourceTable;
    clsSourceRecord.attr("ColumnView") = clsSourceColumnView;
    clsSourceRecord.attr("Catalog") = clsSourceCatalog;
    clsSourceTable.attr("Record") = clsSourceRecord;
    clsSourceTable.attr("ColumnView") = clsSourceColumnView;
    clsSourceTable.attr("Catalog") = clsSourceCatalog;
    clsSourceCatalog.attr("Record") = clsSourceRecord;
    clsSourceCatalog.attr("Table") = clsSourceTable;
    clsSourceCatalog.attr("ColumnView") = clsSourceColumnView;
}

void wrapSimple(WrapperCollection &wrappers) {
    auto clsSimpleRecord = declareSimpleRecord(wrappers);
    auto clsSimpleTable = declareSimpleTable(wrappers);
    auto clsSimpleColumnView = table::python::declareColumnView<SimpleRecord>(wrappers, "Simple");
    auto clsSimpleCatalog = table::python::declareSortedCatalog<SimpleRecord>(wrappers, "Simple");

    clsSimpleRecord.attr("Table") = clsSimpleTable;
    clsSimpleRecord.attr("ColumnView") = clsSimpleColumnView;
    clsSimpleRecord.attr("Catalog") = clsSimpleCatalog;
    clsSimpleTable.attr("Record") = clsSimpleRecord;
    clsSimpleTable.attr("ColumnView") = clsSimpleColumnView;
    clsSimpleTable.attr("Catalog") = clsSimpleCatalog;
    clsSimpleCatalog.attr("Record") = clsSimpleRecord;
    clsSimpleCatalog.attr("Table") = clsSimpleTable;
    clsSimpleCatalog.attr("ColumnView") = clsSimpleColumnView;
}

void wrapSchemaMapper(WrapperCollection &wrappers) {
    wrappers.wrapType(PySchemaMapper(wrappers.module, "SchemaMapper"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<Schema const &, Schema const &>());
        cls.def(py::init<Schema const &, bool>(), "input"_a, "shareAliasMap"_a = false);

        cls.def("getInputSchema", &SchemaMapper::getInputSchema);
        cls.def("getOutputSchema", &SchemaMapper::getOutputSchema);
        cls.def("editOutputSchema", &SchemaMapper::editOutputSchema,
                py::return_value_policy::reference_internal);
        cls.def("addMinimalSchema", &SchemaMapper::addMinimalSchema, "minimal"_a, "doMap"_a = true);
        cls.def_static("removeMinimalSchema", &SchemaMapper::removeMinimalSchema);
        cls.def_static("join", &SchemaMapper::join, "inputs"_a, "prefixes"_a = std::vector<std::string>());

        declareSchemaMapperOverloads<std::uint8_t>(cls, "B");
        declareSchemaMapperOverloads<std::uint16_t>(cls, "U");
        declareSchemaMapperOverloads<std::int32_t>(cls, "I");
        declareSchemaMapperOverloads<std::int64_t>(cls, "L");
        declareSchemaMapperOverloads<float>(cls, "F");
        declareSchemaMapperOverloads<double>(cls, "D");
        declareSchemaMapperOverloads<std::string>(cls, "String");
        declareSchemaMapperOverloads<lsst::geom::Angle>(cls, "Angle");
        declareSchemaMapperOverloads<lsst::afw::table::Flag>(cls, "Flag");
        declareSchemaMapperOverloads<lsst::afw::table::Array<std::uint8_t>>(cls, "ArrayB");
        declareSchemaMapperOverloads<lsst::afw::table::Array<std::uint16_t>>(cls, "ArrayU");
        declareSchemaMapperOverloads<lsst::afw::table::Array<int>>(cls, "ArrayI");
        declareSchemaMapperOverloads<lsst::afw::table::Array<float>>(cls, "ArrayF");
        declareSchemaMapperOverloads<lsst::afw::table::Array<double>>(cls, "ArrayD");
    });
}

void wrapSchema(WrapperCollection &wrappers) {
    // We'll add instantiations of Field, Key, and SchemaItem to these private
    // dicts, and then in schemaContinued.py we'll add them to a TemplateMeta
    // ABC.
    auto &mod = wrappers.module;
    mod.attr("_Field") = py::dict();
    mod.attr("_Key") = py::dict();
    mod.attr("_SchemaItem") = py::dict();

    declareSchemaType<std::uint8_t>(wrappers);
    declareSchemaType<std::uint16_t>(wrappers);
    declareSchemaType<std::int32_t>(wrappers);
    declareSchemaType<std::int64_t>(wrappers);
    declareSchemaType<float>(wrappers);
    declareSchemaType<double>(wrappers);
    declareSchemaType<std::string>(wrappers);
    declareSchemaType<lsst::geom::Angle>(wrappers);
    declareSchemaType<Array<std::uint8_t>>(wrappers);
    declareSchemaType<Array<std::uint16_t>>(wrappers);
    declareSchemaType<Array<int>>(wrappers);
    declareSchemaType<Array<float>>(wrappers);
    declareSchemaType<Array<double>>(wrappers);
    declareSchemaType<Flag>(wrappers);

    declareSchema(wrappers);
    declareSubSchema(wrappers);
}

void wrapMatch(WrapperCollection &wrappers) {
    wrappers.wrapType(py::class_<MatchControl>(wrappers.module, "MatchControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        LSST_DECLARE_CONTROL_FIELD(cls, MatchControl, findOnlyClosest);
        LSST_DECLARE_CONTROL_FIELD(cls, MatchControl, symmetricMatch);
        LSST_DECLARE_CONTROL_FIELD(cls, MatchControl, includeMismatches);
    });

    declareMatch2<SimpleCatalog, SimpleCatalog>(wrappers, "Simple");
    declareMatch2<SimpleCatalog, SourceCatalog>(wrappers, "Reference");
    declareMatch2<SourceCatalog, SourceCatalog>(wrappers, "Source");
    declareMatch1<SimpleCatalog>(wrappers);
    declareMatch1<SourceCatalog>(wrappers);

    wrappers.wrap([](auto &mod) {
        mod.def("matchXy",
                (SourceMatchVector(*)(SourceCatalog const &, SourceCatalog const &, double,
                                      MatchControl const &))matchXy,
                "cat1"_a, "cat2"_a, "radius"_a, "mc"_a = MatchControl());
        mod.def("matchXy", (SourceMatchVector(*)(SourceCatalog const &, double, MatchControl const &))matchXy,
                "cat"_a, "radius"_a, "mc"_a = MatchControl());
    });
}

using PyIdFactory = py::class_<IdFactory, std::shared_ptr<IdFactory>>;

void wrapIdFactory(utils::python::WrapperCollection& wrappers) {
    wrappers.wrapType(PyIdFactory(wrappers.module, "IdFactory"), [](auto& mod, auto& cls) {
        cls.def("__call__", &IdFactory::operator());
        cls.def("notify", &IdFactory::notify, "id"_a);
        cls.def("clone", &IdFactory::clone);
        cls.def_static("makeSimple", IdFactory::makeSimple);
        cls.def_static("makeSource", IdFactory::makeSource, "expId"_a, "reserved"_a);
        cls.def_static("computeReservedFromMaxBits", IdFactory::computeReservedFromMaxBits, "maxBits"_a);
    });
}

void wrapExposure(WrapperCollection &wrappers) {

    auto clsExposureRecord = declareExposureRecord(wrappers);
    auto clsExposureTable = declareExposureTable(wrappers);
    auto clsExposureColumnView = table::python::declareColumnView<ExposureRecord>(wrappers, "Exposure");
    auto clsExposureCatalog = declareExposureCatalog(wrappers);

    clsExposureRecord.attr("Table") = clsExposureTable;
    clsExposureRecord.attr("ColumnView") = clsExposureColumnView;
    clsExposureRecord.attr("Catalog") = clsExposureCatalog;
    clsExposureTable.attr("Record") = clsExposureRecord;
    clsExposureTable.attr("ColumnView") = clsExposureColumnView;
    clsExposureTable.attr("Catalog") = clsExposureCatalog;
    clsExposureCatalog.attr("Record") = clsExposureRecord;
    clsExposureCatalog.attr("Table") = clsExposureTable;
    clsExposureCatalog.attr("ColumnView") = clsExposureColumnView;
}

void wrapAggregates(WrapperCollection &wrappers) {
    wrappers.wrapType(py::enum_<CoordinateType>(wrappers.module, "CoordinateType"), [](auto &mod, auto &enm) {
        enm.value("PIXEL", CoordinateType::PIXEL);
        enm.value("CELESTIAL", CoordinateType::CELESTIAL);
        enm.export_values();
    });

    declarePointKey<double>(wrappers, "2D");
    declarePointKey<int>(wrappers, "2I");

    declareBoxKey<lsst::geom::Box2D>(wrappers, "2D");
    declareBoxKey<lsst::geom::Box2I>(wrappers, "2I");

    declareCoordKey(wrappers);
    declareQuadrupoleKey(wrappers);
    declareEllipseKey(wrappers);

    declareCovarianceMatrixKey<float, 2>(wrappers, "2f");
    declareCovarianceMatrixKey<float, 3>(wrappers, "3f");
    declareCovarianceMatrixKey<float, 4>(wrappers, "4f");
    declareCovarianceMatrixKey<float, Eigen::Dynamic>(wrappers, "Xf");
    declareCovarianceMatrixKey<double, 2>(wrappers, "2d");
    declareCovarianceMatrixKey<double, 3>(wrappers, "3d");
    declareCovarianceMatrixKey<double, 4>(wrappers, "4d");
    declareCovarianceMatrixKey<double, Eigen::Dynamic>(wrappers, "Xd");
}

using PyAliasMap = py::class_<AliasMap, std::shared_ptr<AliasMap>>;

void wrapAliasMap(utils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyAliasMap(wrappers.module, "AliasMap"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<AliasMap const &>());

        cls.def("__len__", &AliasMap::size);
        cls.def("empty", &AliasMap::empty);
        cls.def("apply", &AliasMap::apply, "name"_a);
        cls.def("get", &AliasMap::get, "alias"_a);
        cls.def("__getitem__", &AliasMap::get, "alias"_a);
        cls.def("set", &AliasMap::set, "alias"_a, "target"_a);
        cls.def("__setitem__", &AliasMap::set);
        cls.def("erase", &AliasMap::erase, "alias"_a);
        cls.def("__delitem__", &AliasMap::erase, "alias"_a);
        cls.def("__eq__", [](AliasMap &self, AliasMap &other) { return self == other; });
        cls.def("__ne__", [](AliasMap &self, AliasMap &other) { return self != other; });
        cls.def("contains", &AliasMap::contains, "other"_a);
        cls.def("__contains__", &AliasMap::contains);
        cls.def("items", [](AliasMap &self) { return py::make_iterator(self.begin(), self.end()); },
                py::keep_alive<0, 1>());
    });
}

void wrapArrays(WrapperCollection &wrappers) {
    declareArrayKey<float>(wrappers, "F");
    declareArrayKey<double>(wrappers, "D");
}

void wrapBase(WrapperCollection &wrappers) {

    auto clsBaseTable = declareBaseTable(wrappers);
    auto clsBaseRecord = declareBaseRecord(wrappers);
    auto clsBaseCatalog = table::python::declareCatalog<BaseRecord>(wrappers, "Base");
    auto clsBaseColumnView = table::python::declareColumnView<BaseRecord>(wrappers, "Base");

    clsBaseRecord.attr("Table") = clsBaseTable;
    clsBaseRecord.attr("ColumnView") = clsBaseColumnView;
    clsBaseRecord.attr("Catalog") = clsBaseCatalog;
    clsBaseTable.attr("Record") = clsBaseRecord;
    clsBaseTable.attr("ColumnView") = clsBaseColumnView;
    clsBaseTable.attr("Catalog") = clsBaseCatalog;
    clsBaseCatalog.attr("Record") = clsBaseRecord;
    clsBaseCatalog.attr("Table") = clsBaseTable;
    clsBaseCatalog.attr("ColumnView") = clsBaseColumnView;
}

void wrapBaseColumnView(WrapperCollection &wrappers) {
    declareBaseColumnView(wrappers);
    declareBitsColumn(wrappers);
}
}  // namespace

namespace io {
namespace {
void wrapPersistable(utils::python::WrapperCollection &wrappers) {
    // TODO: uncomment once afw.fits uses WrapperCollection
    // wrappers.addSignatureDependency("lsst.afw.fits");

    wrappers.wrapType(py::class_<Persistable, std::shared_ptr<Persistable>>(wrappers.module, "Persistable"), [](auto &mod, auto &cls) {
        cls.def("writeFits",
                (void (Persistable::*)(std::string const &, std::string const &) const) &
                        Persistable::writeFits,
                "fileName"_a, "mode"_a = "w");
        cls.def("writeFits",
                (void (Persistable::*)(fits::MemFileManager &, std::string const &) const) &
                        Persistable::writeFits,
                "manager"_a, "mode"_a = "w");
        cls.def("isPersistable", &Persistable::isPersistable);
    });
}

void wrapFits(utils::python::WrapperCollection& wrappers) {
    wrappers.wrap([](auto& mod) {
        mod.def("setPreppedRowsFactor",
                [](std::size_t n) { FitsSchemaInputMapper::PREPPED_ROWS_FACTOR = n; });
        mod.def("getPreppedRowsFactor", []() { return FitsSchemaInputMapper::PREPPED_ROWS_FACTOR; });
    });
}

}
}
void wrapTable(py::module_ &mod) {
    WrapperCollection wrappers(mod, "lsst.afw.table");
    wrapAliasMap(wrappers);
    wrapSchema(wrappers);
    wrapSchemaMapper(wrappers);
    wrapBaseColumnView(wrappers);
    wrapBase(wrappers);
    wrapIdFactory(wrappers);
    wrapArrays(wrappers);
    wrapAggregates(wrappers);
    wrapSlots(wrappers);
    wrapSimple(wrappers);
    wrapSource(wrappers);
    wrapExposure(wrappers);
    wrapMatch(wrappers);
    wrapWcsUtils(wrappers);
    wrappers.finish();
    WrapperCollection iowrappers(mod, "lsst.afw.table.io");
    io::wrapPersistable(iowrappers);
    io::wrapFits(iowrappers);
    iowrappers.finish();

}

}
}
}
