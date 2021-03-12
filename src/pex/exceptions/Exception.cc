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

#ifndef __GNUC__
#define __attribute__(x) /*NOTHING*/
#endif

#include <cstring>
#include <ostream>
#include <sstream>
#include <string>

#include "lsst/pex/exceptions/Exception.h"

namespace lsst {
namespace pex {
namespace exceptions {

Tracepoint::Tracepoint(char const* file, int line, char const* func, std::string const& message)
        : _file(file), _line(line), _func(func), _message(message) {}

Exception::Exception(char const* file, int line, char const* func, std::string const& message)
        : _message(message), _traceback(1, Tracepoint(file, line, func, message)) {}

Exception::Exception(std::string const& message) : _message(message), _traceback() {}

Exception::~Exception(void) noexcept {}

void Exception::addMessage(char const* file, int line, char const* func, std::string const& message) {
    std::ostringstream stream;
    stream << _message;
    if (_traceback.empty()) {
        // This means the message-only constructor was used, which should only happen
        // from Python...but this method isn't accessible from Python, so maybe
        // this Exception was thrown in Python, then passed back to C++.  Or, more
        // likely, someone didn't read the warning and threw this exception in C++
        // without using LSST_EXCEPT.
        // But, because we don't have a traceback for that first message, we can't
        // add one here without messing up the stringification logic later.  Since
        // this is a rare case (and should be considered a bug, but we don't want
        // exception code throwing its own exceptions unless it absolutely has to),
        // we'll proceed by just appending the message and ignoring the traceback.
        stream << "; " << message;
    } else {
        if (_traceback.size() == static_cast<std::size_t>(1)) {
            // The original message doesn't have an index (because it's faster,
            // and there's no need if there's only one), so when we add the second,
            // we have to give it an index.
            stream << " {0}; ";
        } else {
            stream << "; ";
        }
        stream << message << " {" << _traceback.size() << "}";
        _traceback.push_back(Tracepoint(file, line, func, message));
    }
    _message = stream.str();
}

Traceback const& Exception::getTraceback(void) const noexcept { return _traceback; }

std::ostream& Exception::addToStream(std::ostream& stream) const {
    if (_traceback.empty()) {
        // The exception was raised in Python, so we don't include the traceback, the type, or any
        // newlines, because Python will print those itself.
        stream << _message;
    } else {
        stream << std::endl;  // Start with a newline to separate our stuff from Pythons "<type>: " prefix.
        for (std::size_t i = 0; i != _traceback.size(); ++i) {
            stream << "  File \"" << _traceback[i]._file << "\", line " << _traceback[i]._line << ", in "
                   << _traceback[i]._func << std::endl;
            stream << "    " << _traceback[i]._message << " {" << i << "}" << std::endl;
        }
        std::string type(getType(), 0, std::strlen(getType()) - 2);
        stream << type << ": '" << _message << "'" << std::endl;
    }
    return stream;
}

char const* Exception::what(void) const noexcept { return _message.c_str(); }

char const* Exception::getType(void) const noexcept { return "lsst::pex::exceptions::Exception *"; }

Exception* Exception::clone(void) const { return new Exception(*this); }

std::ostream& operator<<(std::ostream& stream, Exception const& e) { return e.addToStream(stream); }

}  // namespace exceptions
}  // namespace pex
}  // namespace lsst
