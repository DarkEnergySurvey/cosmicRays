#include "pybind/base_bind.h"
#include "lsst/base/versions.h"
//#include "lsst/base/threads.h"
#include "lsst/base/ModuleImporter.h"

namespace lsst {
namespace base {
namespace {
WRAP(Versions) {
    mod.def("getRuntimeVersions", &lsst::base::getRuntimeVersions);
    mod.def("getCfitsioVersion", &lsst::base::getCfitsioVersion);
    mod.def("getFftwVersion", &lsst::base::getFftwVersion);
    mod.def("getWcslibVersion", &lsst::base::getWcslibVersion);
    mod.def("getGslVersion", &lsst::base::getGslVersion);
}

WRAP(Threads) {
//    mod.def("haveThreads", &lsst::base::haveThreads);
//    mod.def("setNumThreads", &lsst::base::setNumThreads);
//    mod.def("getNumThreads", &lsst::base::getNumThreads);
//    mod.def("disableImplicitThreading", &lsst::base::disableImplicitThreading);
}

class PythonModuleImporter : public ModuleImporter {
public:
    static ModuleImporter const* get() {
        static PythonModuleImporter const instance;
        return &instance;
    }

private:
    PythonModuleImporter() {}

protected:
    virtual bool _import(std::string const& name) const;
};

bool PythonModuleImporter::_import(std::string const& name) const {
    PyObject* mod = PyImport_ImportModule(name.c_str());
    if (mod) {
        Py_DECREF(mod);
        return true;
    } else {
        // If the Python C API call returned a null pointer, it will
        // also have set an exception.  We don't want that, because
        // this isn't necessarily an error (that's up to the caller).
        PyErr_Clear();
    }
    return false;
}

void installPythonModuleImporter() { ModuleImporter::install(PythonModuleImporter::get()); }

WRAP(Cppimport) {
    lsst::base::installPythonModuleImporter();
}

}

WRAP(Base) {
    wrapVersions(mod);
    wrapThreads(mod);
    wrapCppimport(mod);
}
}
}
