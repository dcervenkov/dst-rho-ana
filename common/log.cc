#include "log.h"

// Standard includes
#include <stdarg.h>
#include <iostream>

static Log::LogLevel log_level;

namespace Log {

void setLogLevel(LogLevel level) { log_level = level; }

void print(LogLevel level, const char* fmt, ...) {
    if (level >= log_level) {
        va_list ap;
        va_start(ap, fmt);
        char buf[1000];
        vsnprintf(buf, 1000, fmt, ap);
        std::cout << buf;
        va_end(ap);
    }
}

}  // namespace Log