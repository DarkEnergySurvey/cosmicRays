#include "pybind/pex_bind.h"
#include "lsst/pex/policy.h"
#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace lsst {
namespace pex {
namespace {
//void tryLsstExceptionWarn(const char *message) {
//    // Try to warn that exception translation failed, if we fail, clear the exception raised by the
//    // warning attempt so we can raise a less-informative exception based on the original.
//    int s = PyErr_WarnEx(PyExc_Warning, message, 1);
//    if (s) {
//        PyErr_Clear();
//    }
//}

/*
 * Raise a Python exception that wraps the given C++ exception instance.
 *
 * Most of the work is delegated to the pure-Python function pex.exceptions.wrappers.translate(),
 * which looks up the appropriate Python exception class from a dict that maps C++ exception
 * types to their custom Python wrappers.  Everything else here is basically just importing that
 * module, preparing the arguments, and calling that function, along with the very verbose error
 * handling required by the Python C API.
 *
 * If any point we fail to translate the exception, we print a Python warning and raise the built-in
 * Python RuntimeError exception with the same message as the C++ exception.
 *
 * @param pyex a wrapped instance of pex::exceptions::Exception
 *
void raiseLsstException(py::object &pyex) {
    static auto module =
            py::reinterpret_borrow<py::object>(PyImport_ImportModule("exceptionWrappers"));
    if (!module.ptr()) {
        tryLsstExceptionWarn("Failed to import C++ Exception wrapper module.");
    } else {
        static auto translate =
                py::reinterpret_borrow<py::object>(PyObject_GetAttrString(module.ptr(), "translate"));
        if (!translate.ptr()) {
            tryLsstExceptionWarn("Failed to find translation function for C++ Exceptions.");
        } else {
            // Calling the Python translate() returns an instance of the appropriate Python
            // exception that wraps the C++ exception instance that we give it.
            auto instance = py::reinterpret_steal<py::object>(
                    PyObject_CallFunctionObjArgs(translate.ptr(), pyex.ptr(), NULL));
            if (!instance.ptr()) {
                // We actually expect a null return here, as translate() should raise an exception
                tryLsstExceptionWarn("Failed to translate C++ Exception to Python.");
            } else {
                auto type = py::reinterpret_borrow<py::object>(PyObject_Type(instance.ptr()));
                PyErr_SetObject(type.ptr(), instance.ptr());
            }
        }
    }
}*/

WRAP(Exceptions) {

    py::class_<exceptions::Exception>(mod, "Exception")
            .def(py::init<std::string const &>())
            .def("addMessage", &exceptions::Exception::addMessage)
            .def("getTraceback", &exceptions::Exception::getTraceback)
            .def("addToStream", &exceptions::Exception::addToStream)
            .def("what", &exceptions::Exception::what)
            .def("getType", &exceptions::Exception::getType)
            .def("clone", &exceptions::Exception::clone)
            .def("asString",
                 [](exceptions::Exception &self) -> std::string {
        std::ostringstream stream;
        self.addToStream(stream);
        return stream.str();
    })
    .def("__repr__", [](exceptions::Exception &self) -> std::string {
        std::stringstream s;
        s << "Exception('" << self.what() << "')";
        return s.str();
    });

    py::class_<exceptions::LogicError, exceptions::Exception>(mod, "LogicError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::NotFoundError, exceptions::Exception>(mod, "NotFoundError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::RuntimeError, exceptions::Exception>(mod, "RuntimeError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::IoError, exceptions::RuntimeError>(mod, "IoError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::OverflowError, exceptions::RuntimeError>(mod, "OverflowError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::RangeError, exceptions::RuntimeError>(mod, "RangeError")
            .def(py::init<std::string const &>());

    //py::class_<exceptions::TypeError, exceptions::LogicError>(mod, "TypeError")
    //        .def(py::init<std::string const &>());

    py::class_<exceptions::UnderflowError, exceptions::RuntimeError>(mod, "UnderflowError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::DomainError, exceptions::LogicError>(mod, "DomainError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::InvalidParameterError, exceptions::LogicError>(mod, "InvalidParameterError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::LengthError, exceptions::LogicError>(mod, "LengthError")
            .def(py::init<std::string const &>());

    py::class_<exceptions::OutOfRangeError, exceptions::LogicError>(mod, "OutOfRangeError")
            .def(py::init<std::string const &>());

    //py::register_exception_translator([](std::exception_ptr p) {
    //    try {
    //        if (p) std::rethrow_exception(p);
    //    } catch (const exceptions::Exception &e) {
     //       py::object current_exception;
    //        current_exception = py::cast(e.clone(), py::return_value_policy::take_ownership);
    //        raiseLsstException(current_exception);
    //    }
    //});

    py::register_exception<policy::BadNameError>(mod, "BadNameError", PyExc_RuntimeError);
    //exceptions::python::declareException<policy::BadNameError, lsst::pex::exceptions::RuntimeError>(mod, "BadNameError",
    //                                                                                                           "RuntimeError");
    py::register_exception<policy::DictionaryError>(mod, "DictionaryError", PyExc_RuntimeError);
    //lsst::pex::exceptions::python::declareException<policy::DictionaryError, lsst::pex::exceptions::DomainError>(mod, "DictionaryError",
    //                                                                                                             "DomainError");
    py::register_exception<policy::NameNotFound>(mod, "NameNotFoundError", PyExc_NameError);
    //lsst::pex::exceptions::python::declareException<policy::NameNotFound, lsst::pex::exceptions::NotFoundError>(mod, "NameNotFound",
    //                                                                                                            "NotFoundError");
    py::register_exception<policy::TypeError>(mod, "TypeError", PyExc_TypeError);
    //lsst::pex::exceptions::python::declareException<policy::TypeError, lsst::pex::exceptions::DomainError>(mod, "TypeError", "DomainError");
    py::class_<policy::ValidationError, exceptions::LogicError>clsValidationError(mod, "ValidationError");

    py::enum_<policy::ValidationError::ErrorType>(clsValidationError, "ErrorType")
            .value("OK", policy::ValidationError::ErrorType::OK)
            .value("WRONG_TYPE", policy::ValidationError::ErrorType::WRONG_TYPE)
            .value("MISSING_REQUIRED", policy::ValidationError::ErrorType::MISSING_REQUIRED)
            .value("NOT_AN_ARRAY", policy::ValidationError::ErrorType::NOT_AN_ARRAY)
            .value("ARRAY_TOO_SHORT", policy::ValidationError::ErrorType::ARRAY_TOO_SHORT)
            .value("TOO_FEW_VALUES", policy::ValidationError::ErrorType::TOO_FEW_VALUES)
            .value("TOO_MANY_VALUES", policy::ValidationError::ErrorType::TOO_MANY_VALUES)
            .value("WRONG_OCCURRENCE_COUNT", policy::ValidationError::ErrorType::WRONG_OCCURRENCE_COUNT)
            .value("VALUE_DISALLOWED", policy::ValidationError::ErrorType::VALUE_DISALLOWED)
            .value("VALUE_OUT_OF_RANGE", policy::ValidationError::ErrorType::VALUE_OUT_OF_RANGE)
            .value("BAD_VALUE", policy::ValidationError::ErrorType::BAD_VALUE)
            .value("UNKNOWN_NAME", policy::ValidationError::ErrorType::UNKNOWN_NAME)
            .value("BAD_DEFINITION", policy::ValidationError::ErrorType::BAD_DEFINITION)
            .value("NOT_LOADED", policy::ValidationError::ErrorType::NOT_LOADED)
            .value("UNKNOWN_ERROR", policy::ValidationError::ErrorType::UNKNOWN_ERROR)
            .export_values();

    clsValidationError.def(py::init<const std::string&>());
    clsValidationError.def(py::init<char const*, int, char const*>());

    clsValidationError.def_readonly_static("EMPTY", &policy::ValidationError::EMPTY);
    clsValidationError.def("getErrorMessageFor", &policy::ValidationError::getErrorMessageFor);
    clsValidationError.def("getParamCount", &policy::ValidationError::getParamCount);
    clsValidationError.def("paramNames", &policy::ValidationError::paramNames);
    clsValidationError.def("getParamNames", &policy::ValidationError::getParamNames);
    clsValidationError.def("getErrors",
                           (int (policy::ValidationError::*)(const std::string&) const) & policy::ValidationError::getErrors);
    clsValidationError.def("getErrors", (int (policy::ValidationError::*)() const) & policy::ValidationError::getErrors);
    clsValidationError.def("describe", &policy::ValidationError::describe);
    clsValidationError.def("what", &policy::ValidationError::what);

}
}
WRAP(Pex){
    wrapExceptions(mod);

    py::class_<policy::SupportedFormats, std::shared_ptr<policy::SupportedFormats> >(mod, "SupportedFormats")
            .def(py::init<>())
            .def("registerFormat", &policy::SupportedFormats::registerFormat)
            .def("recognizeType", &policy::SupportedFormats::recognizeType)
            .def("supports", &policy::SupportedFormats::supports)
            .def("size", &policy::SupportedFormats::size);

    py::class_<policy::PolicyStringDestination>(mod, "PolicyStringDestination")
            .def(py::init<>())
            .def(py::init<const std::string&>())
            .def("getData", &policy::PolicyStringDestination::getData);

    py::class_<policy::PolicySource, std::shared_ptr<policy::PolicySource> > cls(mod, "PolicySource");

    py::class_<policy::PolicyString, std::shared_ptr<policy::PolicyString>, policy::PolicySource>(mod, "PolicyString")
            .def(py::init<const std::string&, const policy::SupportedFormats::Ptr&>(), "data"_a,
                 "fmts"_a = policy::PolicyString::defaultFormats)
            .def(py::init<const policy::SupportedFormats::Ptr&>(), "fmts"_a = policy::PolicyString::defaultFormats);

    py::class_<policy::PolicyFile, std::shared_ptr<policy::PolicyFile>, policy::PolicySource> (mod, "PolicyFile")
            // SupportedFormats is not exposed to Python so don't export the default argument
            .def(py::init([](std::string const& filepath) { return new policy::PolicyFile(filepath); }))

            .def("getPath", &policy::PolicyFile::getPath)
            .def("exists", &policy::PolicyFile::exists)
            .def("getFormatName", (const std::string& (policy::PolicyFile::*)()) & policy::PolicyFile::getFormatName)
            .def("load", (void (policy::PolicyFile::*)(policy::Policy&)) & policy::PolicyFile::load)

            .def_readonly_static("EXT_PAF", &policy::PolicyFile::EXT_PAF)
            .def_readonly_static("EXT_XML", &policy::PolicyFile::EXT_XML);

    py::class_<policy::DefaultPolicyFile, std::shared_ptr<policy::DefaultPolicyFile>, policy::PolicyFile> (mod,
                                                                                      "DefaultPolicyFile")
            .def(py::init<const char* const, const std::string&, const std::string&, bool>(), "productName"_a,
                 "filepath"_a, "repos"_a = "", "strict"_a = true)

            .def("load", &policy::DefaultPolicyFile::load)
            .def("getRepositoryPath",
                 [](policy::DefaultPolicyFile const& self) -> std::string { return self.getRepositoryPath().native(); });


    py::class_<policy::UrnPolicyFile, std::shared_ptr<policy::UrnPolicyFile>, policy::DefaultPolicyFile> (mod, "UrnPolicyFile")
            .def(py::init<const std::string&, bool, bool>(), "urn"_a, "strictUrn"_a = false,
                 "strictLoads"_a = true)

            .def_static("productNameFromUrn", &policy::UrnPolicyFile::productNameFromUrn)
            .def_static("filePathFromUrn", &policy::UrnPolicyFile::filePathFromUrn)
            .def_static("reposFromUrn", &policy::UrnPolicyFile::reposFromUrn)
            .def_readonly_static("URN_PREFIX", &policy::UrnPolicyFile::URN_PREFIX)
            .def_readonly_static("URN_PREFIX_ABBREV", &policy::UrnPolicyFile::URN_PREFIX_ABBREV)
            .def_static("looksLikeUrn", &policy::UrnPolicyFile::looksLikeUrn);

    py::class_<policy::paf::PAFWriter> (mod, "PAFWriter")
            .def(py::init<>())
            .def(py::init<const std::string&>())

            .def("writeBools", &policy::paf::PAFWriter::writeBools)
            .def("writeInts", &policy::paf::PAFWriter::writeInts)
            .def("writeDoubles", &policy::paf::PAFWriter::writeDoubles)
            .def("writeStrings", &policy::paf::PAFWriter::writeStrings)
            .def("writePolicies", &policy::paf::PAFWriter::writePolicies)
            .def("writeFiles", &policy::paf::PAFWriter::writeFiles)

            /* Inherited from PolicyWriter */
            .def("write", (void (policy::paf::PAFWriter::*)(const lsst::pex::policy::Policy&, bool)) & policy::paf::PAFWriter::write,
                 "policy"_a, "doDecl"_a = false)
            .def("close", (void (policy::paf::PAFWriter::*)()) & policy::paf::PAFWriter::close)
            .def("toString", (std::string (policy::paf::PAFWriter::*)()) & policy::paf::PAFWriter::toString);


    py::class_<policy::Policy, std::shared_ptr<policy::Policy> > clsPolicy(mod, "Policy");

    py::enum_<policy::Policy::ValueType>(clsPolicy, "ValueType")
            .value("UNDETERMINED", policy::Policy::ValueType::UNDETERMINED)
            .value("UNDEF", policy::Policy::ValueType::UNDEF)
            .value("BOOL", policy::Policy::ValueType::BOOL)
            .value("INT", policy::Policy::ValueType::INT)
            .value("DOUBLE", policy::Policy::ValueType::DOUBLE)
            .value("STRING", policy::Policy::ValueType::STRING)
            .value("POLICY", policy::Policy::ValueType::POLICY)
            .value("FILE", policy::Policy::ValueType::FILE)
            .export_values();

    clsPolicy.def(py::init<>())
            .def(py::init<bool, const policy::Dictionary&>())
            .def(py::init<const std::string&>())
            .def(py::init<policy::Policy&, bool>(), "pol"_a, "deep"_a = false)
            .def(py::init<const policy::PolicySource&>())

            .def_static("createPolicyFromUrn", &policy::Policy::createPolicyFromUrn, "urn"_a, "validate"_a = true)
            .def_static("createPolicy", (policy::Policy * (*)(policy::PolicySource&, bool, bool)) & policy::Policy::createPolicy,
                        "input"_a, "doIncludes"_a = true, "validate"_a = true)
            .def_static("createPolicy",
                        (policy::Policy * (*)(const std::string&, bool, bool)) & policy::Policy::createPolicy, "input"_a,
                        "doIncludes"_a = true, "validate"_a = true)
            .def_static("createPolicy",
                        (policy::Policy * (*)(policy::PolicySource&, const std::string&, bool)) & policy::Policy::createPolicy,
                        "input"_a, "repos"_a, "validate"_a = true)
            .def_static("createPolicy",
                        (policy::Policy * (*)(const std::string&, const std::string&, bool)) & policy::Policy::createPolicy,
                        "input"_a, "repos"_a, "validate"_a = true)
            .def("getValueType",
                 (policy::Policy::ValueType (policy::Policy::*)(const std::string&) const) & policy::Policy::getValueType)
            .def("nameCount", &policy::Policy::nameCount)
            .def("names", (int (policy::Policy::*)(std::list<std::string>&, bool, bool) const) & policy::Policy::names,
                 "names"_a, "topLevelOnly"_a = false, "append"_a = false)
            .def("paramNames",
                 (int (policy::Policy::*)(std::list<std::string>&, bool, bool) const) & policy::Policy::paramNames,
                 "names"_a, "topLevelOnly"_a = false, "append"_a = false)
            .def("policyNames",
                 (int (policy::Policy::*)(std::list<std::string>&, bool, bool) const) & policy::Policy::policyNames,
                 "names"_a, "topLevelOnly"_a = false, "append"_a = false)
            .def("fileNames",
                 (int (policy::Policy::*)(std::list<std::string>&, bool, bool) const) & policy::Policy::fileNames, "names"_a,
                 "topLevelOnly"_a = false, "append"_a = false)
            .def("names", (policy::Policy::StringArray (policy::Policy::*)(bool) const) & policy::Policy::names,
                 "topLevelOnly"_a = false)
            .def("paramNames", (policy::Policy::StringArray (policy::Policy::*)(bool) const) & policy::Policy::paramNames,
                 "topLevelOnly"_a = false)
            .def("policyNames", (policy::Policy::StringArray (policy::Policy::*)(bool) const) & policy::Policy::policyNames,
                 "topLevelOnly"_a = false)
            .def("fileNames", (policy::Policy::StringArray (policy::Policy::*)(bool) const) & policy::Policy::fileNames,
                 "topLevelOnly"_a = false)
            .def("isDictionary", &policy::Policy::isDictionary)
            .def("canValidate", &policy::Policy::canValidate)
            .def("getDictionary", &policy::Policy::getDictionary)
            .def("setDictionary", &policy::Policy::setDictionary)
            // Somehow default arguments don't work here
            .def("validate", [](policy::Policy const& self) { return self.validate(); })
            .def("validate", [](policy::Policy const& self, policy::ValidationError* errs) { return self.validate(errs); })
            .def("valueCount", &policy::Policy::valueCount)
            .def("isArray", &policy::Policy::isArray)
            .def("exists", &policy::Policy::exists)
            .def("isBool", &policy::Policy::isBool)
            .def("isInt", &policy::Policy::isInt)
            .def("isDouble", &policy::Policy::isDouble)
            .def("isString", &policy::Policy::isString)
            .def("isPolicy", &policy::Policy::isPolicy)
            .def("isFile", &policy::Policy::isFile)
            .def("getTypeInfo", &policy::Policy::getTypeInfo)
            .def("getPolicy", (policy::Policy::Ptr (policy::Policy::*)(const std::string&)) & policy::Policy::getPolicy)
            .def("getFile", &policy::Policy::getFile)
            .def("getBool", &policy::Policy::getBool)
            .def("getInt", &policy::Policy::getInt)
            .def("getDouble", &policy::Policy::getDouble)
            .def("getString", &policy::Policy::getString)
            .def("getPolicyArray", &policy::Policy::getPolicyArray)
            .def("getFileArray", &policy::Policy::getFileArray)
            .def("getBoolArray", &policy::Policy::getBoolArray)
            .def("getIntArray", &policy::Policy::getIntArray)
            .def("getDoubleArray", &policy::Policy::getDoubleArray)
            .def("getStringArray", &policy::Policy::getStringArray)
            // _set is called from Python set
            .def("_set", (void (policy::Policy::*)(const std::string&, const policy::Policy::Ptr&)) & policy::Policy::set)
            .def("_set", (void (policy::Policy::*)(const std::string&, bool)) & policy::Policy::set)
            .def("_set", (void (policy::Policy::*)(const std::string&, int)) & policy::Policy::set)
            .def("_set", (void (policy::Policy::*)(const std::string&, double)) & policy::Policy::set)
            .def("_set", (void (policy::Policy::*)(const std::string&, const std::string&)) & policy::Policy::set)
            .def("add", (void (policy::Policy::*)(const std::string&, const policy::Policy::Ptr&)) & policy::Policy::add)
            .def("add", (void (policy::Policy::*)(const std::string&, bool)) & policy::Policy::add)
            .def("add", (void (policy::Policy::*)(const std::string&, int)) & policy::Policy::add)
            .def("add", (void (policy::Policy::*)(const std::string&, double)) & policy::Policy::add)
            .def("add", (void (policy::Policy::*)(const std::string&, const std::string&)) & policy::Policy::add)
            .def("remove", &policy::Policy::remove)
            .def("loadPolicyFiles", (int (policy::Policy::*)(bool)) & policy::Policy::loadPolicyFiles, "strict"_a = true)

            // boost::filesystem::path is equivalent to str in Python
            .def("loadPolicyFiles",
                 [](policy::Policy& self, const std::string& path, bool strict = true) -> int {
        return self.loadPolicyFiles(boost::filesystem::path(path), strict);
    },
    "repository"_a, "strict"_a = true)

            .def("mergeDefaults", &policy::Policy::mergeDefaults, "defaultPol"_a, "keepForValidation"_a = true,
                 "errs"_a = nullptr)
            .def("str", &policy::Policy::str, "name"_a, "indent"_a = "")
            .def("toString", &policy::Policy::toString)
            .def("__str__", &policy::Policy::toString)  // Cleanup stringification later
            .def("asPropertySet", &policy::Policy::asPropertySet);

    py::class_<policy::Dictionary, std::shared_ptr<policy::Dictionary>, policy::Policy>(mod, "Dictionary")
            .def(py::init<>())
            .def(py::init<std::string &>())
            .def(py::init<const policy::Dictionary &>())
            .def(py::init<policy::PolicyFile &>())
            .def(py::init<const policy::Policy &>())
            .def("getDefinitions", (policy::Policy::Ptr (policy::Dictionary::*)()) & policy::Dictionary::getDefinitions)
            .def("check", &policy::Dictionary::check)
            .def("definedNames",
                 (int (policy::Dictionary::*)(std::list<std::string> &, bool) const) & policy::Dictionary::definedNames,
                 "names"_a, "append"_a = false)
            .def("definedNames",
                 (policy::Dictionary::StringArray (policy::Dictionary::*)() const) & policy::Dictionary::definedNames)
            .def("getDef", &policy::Dictionary::getDef)
            .def("makeDef", &policy::Dictionary::makeDef)
            .def("hasSubDictionary", &policy::Dictionary::hasSubDictionary)
            .def("getSubDictionary", &policy::Dictionary::getSubDictionary)
            // For some strange reason pybind11 doesn't like it if we use the default argument here
            .def("validate",
                 [](policy::Dictionary const &self, policy::Policy const &pol) { return self.validate(pol); }, "pol"_a)
            .def("validate", [](policy::Dictionary const &self, policy::Policy const &pol,
                 policy::ValidationError *errs) { self.validate(pol, errs); }, "pol"_a, "errs"_a)
            .def("loadPolicyFiles", (int (policy::Dictionary::*)(bool)) & policy::Dictionary::loadPolicyFiles,
                 "strict"_a = true)
            .def("loadPolicyFiles", [](policy::Dictionary &self, std::string const &repository) {
                 return self.loadPolicyFiles(repository);
                })
            .def("loadPolicyFiles", [](policy::Dictionary &self, std::string const &repository, bool strict) {
                 return self.loadPolicyFiles(repository, strict);
                })
            .def("getPrefix", &policy::Dictionary::getPrefix);

    py::class_<policy::Definition>(mod, "Definition")
            .def(py::init<const std::string &>(), "paramName"_a = "")
            .def("getName", &policy::Definition::getName)
            .def("getPrefix", &policy::Definition::getPrefix)
            .def("setPrefix", &policy::Definition::setPrefix)
            .def("isChildDefinition", &policy::Definition::isChildDefinition)
            .def("setChildDefinition", &policy::Definition::setChildDefinition)
            .def("isWildcard", &policy::Definition::isWildcard)
            .def("setWildcard", &policy::Definition::setWildcard)
            .def("setName", &policy::Definition::setName)
            .def("getData", &policy::Definition::getData)
            .def("setData", &policy::Definition::setData)
            .def("getType", &policy::Definition::getType)
            .def("getTypeName", &policy::Definition::getTypeName)
            .def("getDefault", &policy::Definition::getDefault)
            .def("getDescription", &policy::Definition::getDescription)
            .def("getMaxOccurs", &policy::Definition::getMaxOccurs)
            .def("getMinOccurs", &policy::Definition::getMinOccurs)
            .def("check", &policy::Definition::check);
    py::print("X");
}

}
}
