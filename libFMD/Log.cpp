#include "Log.hpp"

// Set the time format.
const std::string Log::TIME_FORMAT = "[%m-%d-%Y %H:%M:%S] ";

// Have static initializations of all the streams.
CRITICAL_STREAM Log::criticalStream;
ERROR_STREAM Log::errorStream;
WARNING_STREAM Log::warningStream;
OUTPUT_STREAM Log::outputStream;
INFO_STREAM Log::infoStream;
DEBUG_STREAM Log::debugStream;
TRACE_STREAM Log::traceStream;
