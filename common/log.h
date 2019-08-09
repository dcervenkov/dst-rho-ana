/**
 *  @file    log.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2018-09-12
 *
 *  A logging class that uses printf-like formatting while using C++ iostreams
 *  in the background. It also allows adding color tags, time-stamps, etc. to
 *  the messages.
 *
 */

#ifndef LOG_H_
#define LOG_H_

namespace Log {

enum LogLevel { debug, info, warning, error };
void setLogLevel(LogLevel level);
void print(LogLevel level, const char* format, ...);

}  // namespace Log

#endif  // LOG_H_
