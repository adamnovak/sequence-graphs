#ifndef CREATEINDEX_FMDINDEXBUILDER_HPP
#define CREATEINDEX_FMDINDEXBUILDER_HPP

#include <string>

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
         * Keep around an RLCSA builder to build the RLCSA index.
         */
        CSA::RLCSABuilder builder;
        
        /** 
         * How big of a buffer do we reserve in our RLCSA for indexing new
         * sequences? Probably ought to be bigger than the largest single
         * sequence. chr1 is 247,249,719 bp in hg18, so let's put 300,000,000 =
         * 300 mb.
         *
         * Should be about the amount we want to index in a single step, but we
         * can't break sequences, so any sequences that don't fit get their own
         * steps anyway.
         */
        static const size_t BUFFER_SIZE = 300000000;
        
        /**
         * How many threads should we use when building RLCSAs?
         */
        static const size_t THREADS = 10;
};

#endif
