#ifndef INDEXUTIL_HPP
#define INDEXUTIL_HPP

#include <FMDIndex.hpp>
#include <Log.hpp>

#include <string>
#include <vector>

#include <boost/filesystem.hpp>

// indexUtil.hpp: Utility functions for working with FMDIndexes.

/**
 * Start a new index in the given directory (by replacing it), and index the
 * given FASTAs for the bottom level FMD index. Optionally takes a suffix array
 * sample rate to use. Returns the FMD index that gets created.
 */
FMDIndex*
buildIndex(
    std::string indexDirectory,
    std::vector<std::string> fastas,
    int sampleRate = 128
);

#endif
