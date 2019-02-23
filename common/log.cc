#include "log.h"

// Standard includes
#include <stdarg.h>
#include <iostream>

static Log::LogLevel log_level;

namespace Log {

const char* reset = "\033[0m";
const char* red = "\033[31m";
const char* green = "\033[32m";
const char* yellow = "\033[33m";
const char* blue = "\033[34m";

void setLogLevel(LogLevel level) { log_level = level; }

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

}  // namespace Log