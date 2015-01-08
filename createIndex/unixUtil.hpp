#ifndef UNIXUTIL_HPP
#define UNIXUTIL_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

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

/**
 * Save a vector of numbers as a single-column TSV.
 * TODO: Is this UNIX-y enough? Or do we need another util file just for this?
 */
template<typename T>
void
writeColumn(
    std::vector<T> numbers,
    std::string filename
) {

    // Open up the file to write.
    std::ofstream file(filename.c_str());
    
    for(auto number : numbers) {
        // Write each number on its own line
        file << number << std::endl;
    }    
    
    file.close();
}

#endif
