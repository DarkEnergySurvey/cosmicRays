// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

/**
 * @file lsst/pex/policy/exceptions.h
 * @ingroup pex
 * @brief definition of Policy-specific exceptions classes
 * @author Ray Plante
 */

#ifndef LSST_PEX_POLICY_EXCEPTIONS_H
#define LSST_PEX_POLICY_EXCEPTIONS_H

#include "lsst/base.h"
#include "lsst/pex/exceptions.h"

#define POL_EARGS_TYPED char const *ex_file, int ex_line, char const *ex_func
#define POL_EARGS_UNTYPED ex_file, ex_line, ex_func
#define POL_EXCEPT_VIRTFUNCS(etn)                                      \
    char const *etn::getType(void) const throw() { return #etn " *"; } \
    lsst::pex::exceptions::Exception *etn::clone(void) const { return new etn(*this); }

namespace lsst {
namespace pex {
namespace policy {

/**
 * an exception thrown when Policy parameter name has an illegal form.  This
 * usually means that it contains zero-length fields.  In other words, it
 * either starts with a period, ends with a period, or contains two or more
 * consecutive periods.
 */
class LSST_EXPORT BadNameError : public lsst::pex::exceptions::RuntimeError {
public:
    BadNameError(POL_EARGS_TYPED)
            : lsst::pex::exceptions::RuntimeError(POL_EARGS_UNTYPED, "Illegal Policy parameter name") {}
    BadNameError(POL_EARGS_TYPED, const std::string &badname)
            : lsst::pex::exceptions::RuntimeError(POL_EARGS_UNTYPED,
                                                  std::string("Illegal Policy parameter name: ") + badname) {}
    virtual char const *getType(void) const throw();
    virtual lsst::pex::exceptions::Exception *clone() const;
};

/**
 * There is a problem with a dictionary.
 */
class LSST_EXPORT DictionaryError : public lsst::pex::exceptions::DomainError {
public:
    DictionaryError(POL_EARGS_TYPED)
            : lsst::pex::exceptions::DomainError(POL_EARGS_UNTYPED, "Malformed dictionary") {}
    DictionaryError(POL_EARGS_TYPED, const std::string &msg)
            : lsst::pex::exceptions::DomainError(POL_EARGS_UNTYPED,
                                                 std::string("Malformed dictionary: ") + msg) {}
    virtual char const *getType(void) const throw();
    virtual lsst::pex::exceptions::Exception *clone() const;
};

/**
 * an exception indicating that a policy parameter of a given name can
 * not be found in a Policy object.
 */
class LSST_EXPORT NameNotFound : public lsst::pex::exceptions::NotFoundError {
public:
    NameNotFound(POL_EARGS_TYPED)
            : lsst::pex::exceptions::NotFoundError(POL_EARGS_UNTYPED, "Policy parameter name not found") {}
    NameNotFound(POL_EARGS_TYPED, const std::string &parameter)
            : lsst::pex::exceptions::NotFoundError(
                      POL_EARGS_UNTYPED, std::string("Policy parameter name not found: ") + parameter) {}
    virtual char const *getType(void) const throw();
    virtual lsst::pex::exceptions::Exception *clone() const;
};

/**
 * an exception indicating that a policy parameter with a given name has a
 * type different from the one that was requested.
 */
class LSST_EXPORT TypeError : public lsst::pex::exceptions::DomainError {
public:
    TypeError(POL_EARGS_TYPED)
            : lsst::pex::exceptions::DomainError(POL_EARGS_UNTYPED, "Parameter has wrong type") {}
    TypeError(POL_EARGS_TYPED, const std::string &parameter, const std::string &expected)
            : lsst::pex::exceptions::DomainError(POL_EARGS_UNTYPED, std::string("Parameter \"") + parameter +
                                                                            "\" has wrong type; expecting " +
                                                                            expected + ".") {}
    virtual char const *getType(void) const throw();
    virtual lsst::pex::exceptions::Exception *clone() const;
};

}  // namespace policy
}  // namespace pex
}  // namespace lsst

#endif  // LSST_PEX_POLICY_EXCEPTIONS_H
