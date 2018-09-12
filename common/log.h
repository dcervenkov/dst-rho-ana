#ifndef LOG_H_
#define LOG_H_

namespace Log {

enum LogLevel { debug, info, warning, error };
void setLogLevel(LogLevel level);
void print(LogLevel level, const char* format, ...);

}  // namespace Log

#endif  // LOG_H_
