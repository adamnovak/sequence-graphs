#ifndef CREATEINDEX_FMDINDEXBUILDER_HPP
#define CREATEINDEX_FMDINDEXBUILDER_HPP

#include <string>
#include <vector>

#include <rlcsa/rlcsa_builder.h>

/**
 * A class for building an FMD Index with RLCSA. Every index has a "basename",
 * which is a filename prefix to which extensions are appended for actual files.
 *
 * Depends on the RLCSA binaries being available in the 
 */
class FMDIndexBuilder {

    public:
        /**
         * Create a new FMDIndexBuilder using the specified basename for its
         * index. Optionally takes a suffix array sample rate. If an index with
         * that basename already exists, the builder will merge new things into
         * it.
         */
        FMDIndexBuilder(const std::string& basename, int sampleRate = 128);
        
        /**
         * Add the contents of the given FASTA file to the index, both forwards
         * and in reverse complement.
         */
        void add(const std::string& filename);
        
        /**
         * Close all files, sync to disk, and shut down the Builder. Must be
         * called before the index can be read. After this is called, no other
         * method on the same object may be called.
         */
        void close();
    protected:
        /**
         * Keep track of our index basename.
         */
        std::string basename;
        
        /**
         * How often do we make suffix array samples? Really a period.
         */
        int sampleRate;
        
        /**
         * Keep around a vector into which we will attempt to stuff our entire
         * input.
         */
        std::vector<char> input;
        
        /**
         * How many threads should we use when building RLCSAs?
         */
        static const size_t THREADS = 64;
};

#endif
