#ifndef CREATEINDEX_FMDINDEXBUILDER_HPP
#define CREATEINDEX_FMDINDEXBUILDER_HPP

#include <string>

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
         * index. If an index with that basename already exists, the builder
         * will merge new things into it.
         */
        FMDIndexBuilder(const std::string& basename);
        
        /**
         * Add the contents of the given FASTA file to the index, both forwards
         * and in reverse complement.
         */
        void add(const std::string& filename);
    protected:
        /**
         * Merge the index pointed to by otherBasename into our index.
         */
        void merge(const std::string& otherBasename);
        
        /**
         * Keep track of our index basename.
         */
        std::string basename;

};

#endif
