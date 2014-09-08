#ifndef UNIXUTIL_HPP
#define UNIXUTIL_HPP

/**
 * unixUtil.hpp: utility functions for doing useful things with Unix.
 */

/**
 * Log current memory usage at INFO level.
 */
void logMemory();

/**
 * Signal handler to exit with exit(), preserving gprof profiling data if
 * applicable.
 */
void exitOnSignal(int signalNumber);

/**
 * Signal handler for printing a stack trace and exiting on a signal. See
 * <http://stackoverflow.com/a/77336/402891>
 */
void stacktraceOnSignal(int signalNumber);
#endif
