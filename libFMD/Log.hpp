#ifndef LOG_HPP
#define LOG_HPP

/**
 * Static logging system.
 */

#include <string>
#include <iostream>
#include <ostream>
#include <ctime>
#include <functional>
#include <sstream>


// Configure log levels to compile in and out. Levels with OutStreams go to
// stdout, levels with NullStreams go nowhere.

// Critical: Log messages that stop the program.
#define CRITICAL_STREAM OutStream
// Error: Log messages that indicate something is broken
#define ERROR_STREAM OutStream
// Warning: Log messages that indicate something is unlikely.
#define WARNING_STREAM OutStream
// Output: Messages the user is supposed to see
#define OUTPUT_STREAM OutStream
// Info: Messages about the operation of the program
#define INFO_STREAM OutStream
// Debug: Messages about the internals of the program.
#define DEBUG_STREAM OutStream
// Trace: Messages that pedantically describe what the program is doing.
#define TRACE_STREAM NullStream

// Also we have some macros to facilitate lazy logging. If you use these, your
// things you log are not evaluated unless that log level is enabled.

// Capture everything by reference and when evaluated send the expression to a
// stringstream. Make sure not to put the expression in parens so it can daisy-
// chain << operators. We also need to make sure to convert to std::function so
// later we can use the correct template (not the generic one for all types).
// TODO: Maybe use <http://stackoverflow.com/a/6928769/402891> to detect
// callableness instead and save a std::function construction?
// Use like:
// Log::debug() << "I know " << LOG_LAZY("Some" << "Stuff") << std::endl;
#define LOG_LAZY(exp) (std::function<void(std::stringstream&)>)\
    ([&](std::stringstream& stream) -> void { stream << exp; })

// We also have some log level macros to make that easy. Use like:
// LOG_INFO("i is " << i << std::endl);
// TODO: Make these no-ops early if using NullStream?
#define LOG_CRITICAL(exp) Log::critical() << LOG_LAZY(exp);
#define LOG_ERROR(exp) Log::error() << LOG_LAZY(exp);
#define LOG_WARNING(exp) Log::warning() << LOG_LAZY(exp);
#define LOG_OUTPUT(exp) Log::output() << LOG_LAZY(exp);
#define LOG_INFO(exp) Log::info() << LOG_LAZY(exp);
#define LOG_DEBUG(exp) Log::debug() << LOG_LAZY(exp);
#define LOG_TRACE(exp) Log::trace() << LOG_LAZY(exp);

/**
 * Define a typedef for the type of ostream manipulators like std::endl. They
 * are really templates, but they can template down to this type, and this one
 * works on cout, so it's OK to just use this.
 */
typedef std::ostream& (*OstreamManipulator)(std::ostream&);

/**
 * A pretend stream discards all output.
 */
class NullStream {
public:
    /**
     * Define a template operator that eats anything you try to throw at it.
     */
    template<typename T>
    NullStream& operator<<(T thing) {
        return *this;
    }
    
    /**
     * Define another template that passes stream manipulators (which are really
     * functions) on to the underlying stream.
     */
    inline NullStream& operator<<(OstreamManipulator manipulator) {
        // Ignore the manipulator
        return *this;
    }
    
    /**
     * Define another function that takes lambdas of stringstreams and ignores
     * them.
     */
    inline NullStream& operator<<(
        std::function<void(std::stringstream&)> message) {
        
        // Ignore the lambda
        return *this;
    }
};

/**
 * A stream that sends output to stdout.
 */
class OutStream {
public:
    /**
     * Define a template operator that sends everything to standard output.
     */
    template<typename T>
    OutStream& operator<<(T thing) {
        std::cout << thing;
        return *this;
    }
    
    /**
     * Define another function that passes stream manipulators (which are really
     * functions) on to the underlying stream.
     */
    inline OutStream& operator<<(OstreamManipulator manipulator) {
        // Send the manipulator to standard output.
        std::cout << manipulator;
        return *this;
    }
    
    /**
     * Define another function that takes lambdas of stringstreams and evaluates
     * them.
     */
    inline OutStream& operator<<(
        std::function<void(std::stringstream&)> message) {
        
        // Send the message to a stringstream and then to the output.
        
        std::stringstream stream;
        message(stream);
        std::cout << stream.str();
        
        return *this;
    }
     
    
};

/**
 * A class that holds static methods for logging in a stream-like way.
 */
class Log {
public:
    /**
     * Function to get a stream to log CRITICAL-level massages to. Call once per
     * line.
     */
    static inline CRITICAL_STREAM& critical() {
        return criticalStream << timestamp << "CRITICAL: ";
    }
    
    /**
     * Function to get a stream to log ERROR-level massages to. Call once per
     * line.
     */
    static inline ERROR_STREAM& error() {
        return errorStream << timestamp << "ERROR: ";
    }
    
    /**
     * Function to get a stream to log WARNING-level massages to. Call once per
     * line.
     */
    static inline WARNING_STREAM& warning() {
        return errorStream << timestamp << "WARNING: ";
    }
    
    /**
     * Function to get a stream to log OUTPUT-level massages to. Call once per
     * line.
     */
    static inline OUTPUT_STREAM& output() {
        return outputStream << timestamp << "OUTPUT: ";
    }
    
    /**
     * Function to get a stream to log INFO-level massages to. Call once per
     * line.
     */
    static inline INFO_STREAM& info() {
        return infoStream << timestamp << "INFO: ";
    }
    
    /**
     * Function to get a stream to log INFO-level massages to. Call once per
     * line.
     */
    static inline DEBUG_STREAM& debug() {
        return debugStream << timestamp << "DEBUG: ";
    }
    
    /**
     * Function to get a stream to log TRACE-level massages to. Call once per
     * line.
     */
    static inline TRACE_STREAM& trace() {
        return traceStream << timestamp << "TRACE: ";
    }
    
private:
    
    // Static streams for all the log levels.
    static CRITICAL_STREAM criticalStream;
    static ERROR_STREAM errorStream;
    static WARNING_STREAM warningStream;
    static OUTPUT_STREAM outputStream;
    static INFO_STREAM infoStream;
    static DEBUG_STREAM debugStream;
    static TRACE_STREAM traceStream;
    
    // How long can a time be?
    static const int TIME_CHARS = 80;
    static const std::string TIME_FORMAT;
    
    /**
     * We define a manipulator which manipulates a stream into printing the
     * current time. This way we only do the time formatting when an actual
     * stream is getting printed to.
     */
    static inline std::ostream& timestamp(std::ostream& stream) {
        // Get the current time. See <http://stackoverflow.com/a/16358264>
        
        // This holds the global time
        time_t globalTime;
        time(&globalTime);
        
        // This holds the local time
        struct tm * localTime = localtime(&globalTime);
        
        // This holds the formatted time
        char buffer[TIME_CHARS];
        if(strftime(buffer, TIME_CHARS, TIME_FORMAT.c_str(), localTime)) {
            // Send the time to the stream. It is null-terminated.
            stream << buffer;
        } else {
            // Complain (in the log) that we can't format the time. Maybe it's
            // now the year 70 billion.
            stream << "(Time too long) ";
        }
        
        // Pass on the stream.
        return stream;

    }
    
    // Don't ever let anyone make a Log object.
    Log();
};

#endif
