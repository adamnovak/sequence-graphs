#ifndef CREDITFILTER2_HPP
#define CREDITFILTER2_HPP

#include <vector>
#include "Mapping.hpp"
#include "FMDIndex.hpp"
#include "DisambiguateFilter.hpp"

/**
 * Represents a filter for mappings which applies credit. Lives in its own class
 * so it is easy to unit test.
 */
class CreditFilter2 {

public:
    
    /**
     * Make a new credit filter that uses the given index to get context
     * lengths, with the given bitvector of ranges for counting words, and using
     * the given max number of mismatches for credit assignment.
     */
    CreditFilter2(const FMDIndex& index, const GenericBitvector& ranges,
        size_t z_max);
    
    /**
     * Given a vector of left mappings (in right semantics with left context
     * set, in left to right order) and an equal-length vector of right mappings
     * (also in right semantics, with right context set, in left-to-right
     * order), and also the qquery string that was mapped, apply credit. Produce
     * a vector of disambiguated mappings with credit. TODO: needs to check base
     * identity.
     */
    std::vector<Mapping> apply(const std::vector<Mapping>& leftMappings, 
        const std::vector<Mapping>& rightMappings, const std::string& query);
        
protected:

    // What index knows how long all our contigs are?
    const FMDIndex& index;
    
    // Where are the range boundaries for counting word occurrences?
    const GenericBitvector& ranges;
    
    // How many mismatches are we allowing in sentinel words or maximal left or
    // right contexts?
    size_t z_max;
    
    // We need to be able to disambiguate things.
    DisambiguateFilter disambiguate;

};

#endif
