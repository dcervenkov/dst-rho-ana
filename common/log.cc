/**
 *  @file    log.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2018-09-12
 *
 *  A logging class that uses printf-like formatting while using C++ iostreams
 *  in the background. It also allows adding color tags, time-stamps, etc. to
 *  the messages.
 *
 */

#include "log.h"

// Standard includes
#include <stdarg.h>
#include <iostream>

static Log::LogLevel log_level;

namespace Log {

// ANSI terminal formatting sequences
const char* reset = "\033[0m";
const char* red = "\033[31m";
const char* green = "\033[32m";
const char* yellow = "\033[33m";
const char* blue = "\033[34m";

/**
 * Set the lowest message level that should be printed. Messages with lower
 * level will be silently ignored.
 *
 * @param level Lowest message level to be outputted
 */
void setLogLevel(LogLevel level) { log_level = level; }

/**
 * Print function that behaves very much like printf, with the added log level
 * argument.
 *
 * @param level Log level of the message
 * @param fmt printf-like formatting string
 * @param ... Optional parameters, much like printf has
 */
void print(LogLevel level, const char* fmt, ...) {
    if (level >= log_level) {
        va_list ap;
        va_start(ap, fmt);
        char buf[1000];
        vsnprintf(buf, 1000, fmt, ap);
        switch (level) {
            case debug:
                std::cout << blue << "[DEBUG] " << reset;
                break;
            case info:
                std::cout << green << "[INFO] " << reset;
                break;
            case warning:
                std::cout << yellow << "[WARNING] " << reset;
                break;
            case error:
                std::cout << red << "[ERROR] " << reset;
                break;
            default:
                break;
        }
        std::cout << buf;
        va_end(ap);
    }
}

/**
 * std::cout-like constructor. To be used the same way as std::cout.
 *
 * @param level Log level of the message
 * @param out Stream to which the message should be sent
 */
LogLine::LogLine(LogLevel level, std::ostream& out) : out_(out), level_(level) {
    if (level_ >= log_level) {
        switch (level) {
            case debug:
                out << blue << "[DEBUG] " << reset;
                break;
            case info:
                out << green << "[INFO] " << reset;
                break;
            case warning:
                out << yellow << "[WARNING] " << reset;
                break;
            case error:
                out << red << "[ERROR] " << reset;
                break;
            default:
                break;
        }
    }
}

LogLine::~LogLine() {
    if (level_ >= log_level) {
        stream_ << "\n";
        out_ << stream_.rdbuf();
        out_.flush();
    }
}
}  // namespace Log
