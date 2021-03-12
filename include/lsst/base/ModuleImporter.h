// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef LSST_BASE_ModuleImporter_h_INCLUDED
#define LSST_BASE_ModuleImporter_h_INCLUDED

/**
 *  @file
 *
 *  Mechanism for safely importing Python modules from C++; should not be included
 *  except by its own implementation file, the ioLib.i file, and Persistable.cc.
 */

#include <string>

namespace lsst { namespace base {

/**
 *  @brief Base class that defines an interface for importing Python modules.
 *
 *  The default implementation (defined in the source file) simply returns
 *  false, indicating that it can't import the given module.  The functional
 *  implementation is in the ioLib Swig module, which is installed when that
 *  module is imported.  That machinery keeps us from calling Python C-API
 *  functions from standalone C++ binaries that aren't linked with Python.
 */
class ModuleImporter {
public:

    /// Import the given Python module, and return true if successful.
    static bool import(std::string const & name);

    // No copying
    ModuleImporter (const ModuleImporter&) = delete;
    ModuleImporter& operator=(const ModuleImporter&) = delete;

    // No moving
    ModuleImporter (ModuleImporter&&) = delete;
    ModuleImporter& operator=(ModuleImporter&&) = delete;

protected:

    ModuleImporter() {}

    virtual bool _import(std::string const & name) const = 0;

    virtual ~ModuleImporter() {}

private:

    friend void installPythonModuleImporter();

    static void install(ModuleImporter const * importer);

};

}} // namespace lsst::base

#endif // !LSST_BASE_ModuleImporter_h_INCLUDED
