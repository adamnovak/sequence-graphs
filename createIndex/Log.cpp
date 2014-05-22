#include "Log.hpp"

// Set the time format.
const std::string Log::TIME_FORMAT = "[%m-%d-%Y %I:%M:%S] ";

// Have static initializations of all the streams.
CRITICAL_STREAM Log::criticalStream;
ERROR_STREAM Log::errorStream;
OUTPUT_STREAM Log::outputStream;
INFO_STREAM Log::infoStream;
DEBUG_STREAM Log::debugStream;
TRACE_STREAM Log::traceStream;
