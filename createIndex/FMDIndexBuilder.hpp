#ifndef CREATEINDEX_FMDINDEXBUILDER_HPP
#define CREATEINDEX_FMDINDEXBUILDER_HPP

#include <string>
#include <iostream>

/**
 * A class for building an FMD Index with libsuffixtools. Every index has a
 * "basename", which is a filename prefix to which extensions are appended for
 * actual files.
 */
class FMDIndexBuilder {

    public:
        /**
         * Create a new FMDIndexBuilder using the specified basename for its
         * index. If an index with that basename already exists, it will be
         * replaced. Optionally, you can specify a suffix array sample rate.
         */
        FMDIndexBuilder(const std::string& basename, int sampleRate = 64);
        
        /**
         * Add the contents of the given FASTA file to the index, both forwards
         * and in reverse complement.
         */
        void add(const std::string& filename);
        
        /**
         * Build the final index, close all files, sync to disk, and shut down
         * the FMDIndexBuilder. Must be called before the index can be read.
         * After this is called, no other method on the same object may be
         * called.
         */
        void close();
    protected:
        /**
         * Keep track of our index basename.
         */
        std::string basename;
        
        /**
         * Keep a temporary directory for intermediate files.
         */
        std::string tempDir;
        
        /**
         * Keep the name of the file we're saving the contigs in.
         */
        std::string tempFastaName;
        
        /**
         * Keep around a file to save the contigs in.
         */
        std::ofstream tempFasta;
        
        /**
         * Keep around a file to save the contig names and lengths in.
         */
        std::ofstream contigFile;
        
        /**
         * Keep track of the sample rate to use when we produce the sampled
         * suffix array (when close() is called).
         */
        int sampleRate;
        
        /**
         * How many threads should we use when building the index?
         */
        static const size_t NUM_THREADS = 10;
};

#endif
