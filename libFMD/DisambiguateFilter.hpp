#ifndef DISAMBIGUATEFILTER_HPP
#define DISAMBIGUATEFILTER_HPP

#include <vector>
#include "Mapping.hpp"
#include "FMDIndex.hpp"

/**
 * Represents a filter which disambiguates left and right Mappings into a single
 * Mapping.
 */
class DisambiguateFilter {

public:

    /**
     * Make a new DisambiguateAFilter using the given index for contig lengths.
     */
    DisambiguateFilter(const FMDIndex& index);

    /**
     * Apply the filter. If only mapped on the left, take that. If only mapped
     * on the right, flip it and take that. If mapped on the left and the right
     * to the same place, take that. Otherwise, return unmapped.
     */
    std::vector<Mapping> apply(std::vector<Mapping> leftMappings, 
        std::vector<Mapping> rightMappings);
        
protected:

    // What index knows how long all our contigs are?
    const FMDIndex& index;
        
};

#endif
