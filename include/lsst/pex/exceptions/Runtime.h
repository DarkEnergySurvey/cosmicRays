/*
 * This file is part of pex_exceptions.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#if !defined(LSST_RUNTIME_EXCEPTION)
#define LSST_RUNTIME_EXCEPTION 1

#include "lsst/pex/exceptions/Exception.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace pex {
namespace exceptions {

/**
 * Reports errors in the logical structure of the program.
 *
 * LogicError and its subclasses should be thrown to represent problems, such
 * as violation of logical preconditions or class invariants, that are in
 * principle preventable using defensive programming or other good practices.
 * In many cases, it may not be appropriate to catch them.
 *
 * @see RuntimeError
 * @see std::logic_error
 */
LSST_EXCEPTION_TYPE(LogicError, Exception, lsst::pex::exceptions::LogicError)

/**
 * Reports arguments outside the domain of an operation.
 *
 * This exception should be reserved for mathematical operations that are
 * defined on a limited range of inputs. InvalidParameterError is more
 * appropriate for non-mathematical operations.
 *
 * @see std::domain_error
 */
LSST_EXCEPTION_TYPE(DomainError, LogicError, lsst::pex::exceptions::DomainError)

/**
 * Reports invalid arguments.
 *
 * This exception reports errors that arise because an argument value has not been accepted.
 *
 * @see std::invalid_argument
 */
LSST_EXCEPTION_TYPE(InvalidParameterError, LogicError, lsst::pex::exceptions::InvalidParameterError)

/**
 * Reports attempts to exceed implementation-defined length limits for some classes.
 *
 * For example, some collection classes might not be able to handle more than
 * some number of elements, or bit fields might support a limited number of flags.
 *
 * @see std::length_error
 */
LSST_EXCEPTION_TYPE(LengthError, LogicError, lsst::pex::exceptions::LengthError)

/**
 * Reports attempts to access elements outside a valid range of indices.
 *
 * @see NotFoundError
 * @see std::out_of_range
 *
 * @note pybind11 wrappers should manually translate this exception to
 * `py::index_error` when appropriate. Some Python language constructs check
 * for exceptions that are exactly `IndexError` rather than a
 * sub- or superclass.
 */
LSST_EXCEPTION_TYPE(OutOfRangeError, LogicError, lsst::pex::exceptions::OutOfRangeError)

/**
 * Reports errors that are due to events beyond the control of the program.
 *
 * RuntimeError and its subclasses represent problems that cannot be easily
 * predicted or prevented. In other words, a RuntimeError is a possible outcome
 * of calling a function or method even in well-written programs, and should be
 * handled at the appropriate level.
 *
 * In Python, this exception inherits from `builtins.RuntimeError`.
 *
 * @see LogicError
 * @see std::runtime_error
 */
LSST_EXCEPTION_TYPE(RuntimeError, Exception, lsst::pex::exceptions::RuntimeError)

/**
 * Reports when the result of an operation cannot be represented by the destination type.
 *
 * Situations covered by this exception include lossy type conversions.
 *
 * @see OverflowError
 * @see UnderflowError
 * @see std::range_error
 */
LSST_EXCEPTION_TYPE(RangeError, RuntimeError, lsst::pex::exceptions::RangeError)

/**
 * Reports when the result of an arithmetic operation is too large for the destination type.
 *
 * In Python, this exception inherits from `builtins.OverflowError`.
 *
 * @see std::overflow_error
 */
LSST_EXCEPTION_TYPE(OverflowError, RuntimeError, lsst::pex::exceptions::OverflowError)

/**
 * Reports when the result of an arithmetic operation is too small for the destination type.
 *
 * In Python, this exception inherits from `builtins.ArithmeticError`.
 *
 * @see std::underflow_error
 */
LSST_EXCEPTION_TYPE(UnderflowError, RuntimeError, lsst::pex::exceptions::UnderflowError)

/**
 * Reports attempts to access elements using an invalid key.
 *
 * This exception may represent lookup failures in classes that resemble C++
 * maps or Python dictionaries, but it may also be used when the relationship
 * between an identifier and a resource is more abstract.
 *
 * In Python, this exception inherits from `builtins.LookupError`.
 *
 * @see OutOfRangeError
 *
 * @note pybind11 wrappers should manually translate this exception to
 * `py::key_error` when appropriate. Some Python language constructs check
 * for exceptions that are exactly `KeyError` rather than a
 * sub- or superclass.
 */
LSST_EXCEPTION_TYPE(NotFoundError, Exception, lsst::pex::exceptions::NotFoundError)

/**
 * Reports errors in external input/output operations.
 *
 * In Python, this exception inherits from `builtins.IOError`.
 *
 * @see std::ios_base::failure
 */
LSST_EXCEPTION_TYPE(IoError, RuntimeError, lsst::pex::exceptions::IoError)

/**
 * Reports errors from accepting an object of an unexpected or inappropriate type.
 *
 * In Python, this exception inherits from `builtins.TypeError`.
 */
LSST_EXCEPTION_TYPE(TypeError, LogicError, lsst::pex::exceptions::TypeError)

}  // namespace exceptions
}  // namespace pex
}  // namespace lsst

#endif
