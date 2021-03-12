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

#ifndef LSST_PEX_EXCEPTIONS_EXCEPTION_H
#define LSST_PEX_EXCEPTIONS_EXCEPTION_H

#include <exception>
#include <ostream>
#include <string>
#include <vector>

#include "lsst/base.h"
#include "boost/current_function.hpp"

namespace lsst {
namespace pex {
namespace exceptions {

/// For internal use; gathers the file, line, and function for a tracepoint.
#define LSST_EXCEPT_HERE __FILE__, __LINE__, BOOST_CURRENT_FUNCTION

/**
 * Create an exception with a given type.
 *
 * @param[in] type C++ type of the exception to be thrown.
 * @param[in] ... The message, and optionally other arguments (dependent on the type).
 */
#define LSST_EXCEPT(type, ...) type(LSST_EXCEPT_HERE, __VA_ARGS__)

/**
 * @brief Add the current location and a message to an existing exception before
 * rethrowing it.
 */
#define LSST_EXCEPT_ADD(e, m) e.addMessage(LSST_EXCEPT_HERE, m)

/// The initial arguments required for new exception subclasses.
#define LSST_EARGS_TYPED char const *ex_file, int ex_line, char const *ex_func, std::string const &ex_message

/// The initial arguments to the base class constructor for new subclasses.
#define LSST_EARGS_UNTYPED ex_file, ex_line, ex_func, ex_message

/**
 * Macro used to define new types of exceptions without additional data.
 *
 * @param[in] t Type of the exception.
 * @param[in] b Base class of the exception.
 * @param[in] c C++ class of the exception (fully specified).
 */
#define LSST_EXCEPTION_TYPE(t, b, c)                                                          \
    class LSST_EXPORT t : public b {                                                                      \
    public:                                                                                   \
        t(LSST_EARGS_TYPED) : b(LSST_EARGS_UNTYPED){};                                        \
        t(std::string const& message) : b(message){};                                         \
        virtual char const* getType(void) const noexcept { return #c " *"; };                 \
        virtual lsst::pex::exceptions::Exception* clone(void) const { return new t(*this); }; \
    };

/// One point in the Traceback vector held by Exception
struct Tracepoint {
    /**
     * Standard constructor, intended for C++ use.
     *
     * @param[in] file Filename.
     * @param[in] line Line number.
     * @param[in] func Function name.
     * @param[in] message Informational string attached to exception.
     */
    Tracepoint(char const* file, int line, char const* func, std::string const& message);

    char const* _file;  // Compiled strings only; does not need deletion
    int _line;
    char const* _func;  // Compiled strings only; does not need deletion
    std::string _message;
};
typedef std::vector<Tracepoint> Traceback;

/**
 * Provides consistent interface for LSST exceptions.
 *
 * All exceptions defined by the LSST Stack are derived from this class. Code
 * should not throw or catch Exception directly, but should instead be written
 * in terms of the appropriate subclasses (e.g., catch RuntimeError to handle
 * all unknown errors).
 *
 * In Python, this exception inherits from `builtins.Exception`.
 */
class LSST_EXPORT Exception : public std::exception {
public:
    /**
     * Standard constructor, intended for C++ use via the LSST_EXCEPT() macro.
     *
     * @param[in] file Filename (automatically passed in by macro).
     * @param[in] line Line number (automatically passed in by macro).
     * @param[in] func Function name (automatically passed in by macro).
     * @param[in] message Informational string attached to exception.
     */
    Exception(char const* file, int line, char const* func,
              std::string const& message);  // Should use LSST_EARGS_TYPED, but that confuses doxygen.

    /**
     * Message-only constructor, intended for use from Python only.
     *
     * While this constructor can be called from C++, it's better to use the @ref LSST_EXCEPT
     * macro there, which stores file/line/function information as well.  In Python, however,
     * that information is stored outside the exception, so we don't want to duplicate it,
     * and hence this constructor is invoked instead.
     *
     * @param[in] message Informational string attached to exception.
     */
    explicit Exception(std::string const& message);

    virtual ~Exception(void) noexcept;

    /**
     * Add a tracepoint and a message to an exception before rethrowing it (access via @ref LSST_EXCEPT_ADD).
     *
     * @param[in] file Filename (automatically passed in by macro).
     * @param[in] line Line number (automatically passed in by macro).
     * @param[in] func Function name (automatically passed in by macro).
     * @param[in] message Additional message to associate with this rethrow.
     */
    void addMessage(char const* file, int line, char const* func, std::string const& message);

    /// Retrieve the list of tracepoints associated with an exception.
    Traceback const& getTraceback(void) const noexcept;

    /**
     * @brief Add a text representation of this exception, including its traceback with
     * messages, to a stream.
     *
     * @param[in] stream Reference to an output stream.
     * @returns Reference to the output stream after adding the text.
     */
    virtual std::ostream& addToStream(std::ostream& stream) const;

    /**
     * Return a character string summarizing this exception.
     *
     * This combines all the messages added to the exception, but not the type or
     * traceback (use the stream operator to get this more detailed information).
     *
     * Not allowed to throw any exceptions.
     *
     * @returns String representation; does not need to be freed/deleted.
     */
    virtual char const* what(void) const noexcept;

    /**
     * Return the fully-specified C++ type of a pointer to the exception.  This
     * is overridden by derived classes (automatically if the @ref LSST_EXCEPTION_TYPE
     * macro is used).  It is used by the Python interface.
     *
     * @returns String with the C++ type; does not need to be freed/deleted.
     */
    virtual char const* getType(void) const noexcept;

    /**
     * Return a copy of the exception as an Exception pointer. Can be overridden by
     * derived classes that add data or methods.
     *
     * @returns Pointer to a copy of the exception.
     */
    virtual Exception* clone(void) const;

private:
    std::string _message;
    Traceback _traceback;
};

/**
 * Push the text representation of an exception onto a stream.
 *
 * @param[in] stream Reference to an output stream.
 * @param[in] e Exception to output.
 * @returns Reference to the output stream after adding the text.
 */
std::ostream& operator<<(std::ostream& stream, Exception const& e);
}
}  // namespace pex
}  // namespace lsst

#endif
