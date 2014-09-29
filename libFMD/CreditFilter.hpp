#ifndef CREDITFILTER_HPP
#define CREDITFILTER_HPP

#include <vector>
#include "Mapping.hpp"
#include "FMDIndex.hpp"
#include "DisambiguateFilter.hpp"

/**
 * Represents a filter for mappings which applies credit. Lives in its own class
 * so it is easy to unit test.
 */
class CreditFilter {

public:
    
    /**
     * Make a new credit filter that uses the given index to get context
     * lengths.
     */
    CreditFilter(const FMDIndex& index);
    
    /**
     * Given a vector of left mappings (in right semantics with left context
     * set, in left to right order) and an equal-length vector of right mappings
     * (also in right semantics, with right context set, in left-to-right
     * order), apply credit. Produce a vector of disambiguated mappings with
     * credit. TODO: needs to check base identity.
     */
    std::vector<Mapping> apply(std::vector<Mapping> leftMappings, 
        std::vector<Mapping> rightMappings);
        
protected:

    // What index knows how long all our contigs are?
    const FMDIndex& index;
    
    // We need to be able to disambiguate things.
    DisambiguateFilter disambiguate;

};

#endif
