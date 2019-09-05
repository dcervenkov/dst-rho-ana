/**
 *  @file    log.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2018-09-12
 *
 *  A logging framework that supports two different approaches.
 *  1) printf-like formatting while using C++ iostreams in the background.
 *     Usage: Log::print(Log::debug, "Message: %s", c_str);
 *  2) std::cout like approach which can print things like std::string
 *     Usage: Log::LogLine(Log::debug) << string;
 *
 * It also allows adding color tags, time-stamps,
 *  etc. to the messages.
 *
 */

#pragma once

#include <iostream>
#include <sstream>

namespace Log {

enum LogLevel { debug, info, warning, error };
void setLogLevel(LogLevel level);
void print(LogLevel level, const char* format, ...);

/**
 * Class used for std::cout like printing of stamped log messages.
 *
 * Usage: Log::LogLine(Log::debug) << string;
 *
 */
class LogLine {
   public:
    LogLine(LogLevel level, std::ostream& out = std::cout);
    ~LogLine();

    template <class T>
    LogLine& operator<<(const T& thing) {
        stream_ << thing;
        return *this;
    }

   private:
    std::stringstream stream_;
    std::ostream& out_;
    // static LogFilter...
    LogLevel level_;
};

}  // namespace Log
