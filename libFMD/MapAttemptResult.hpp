#ifndef MAPATTEMPTRESULT_HPP
#define MAPATTEMPTRESULT_HPP

#include "GenericBitVector.hpp"

/**
 * A triple to hold the return values from FMD::mapPosition() or
 * FMD::partialMap(). Holds a flag for whether the mapping succeeded or not, an
 * FMDPosition corresponding either to where the character mapped or the longest
 * search starting at the character that did actually return results, and the
 * number of characters in the FMDPosition's search pattern.
 */
struct MapAttemptResult
{
    bool is_mapped;
    FMDPosition position;
    size_t characters;
};

struct creditMapAttemptResult
{
    bool is_mapped;
    FMDPosition position;
    size_t characters;
    size_t maxCharacters;
};

struct MisMatchAttemptResults
{

    MisMatchAttemptResults(): is_mapped(0), positions(), characters(0) {
        // Nothing to do. No positions at all is fine.
    }

    bool is_mapped;
    
    // Holds pairs of result set and mismatch count.
    std::vector<std::pair<FMDPosition,size_t>> positions;
    
    // How many characters have been searched?
    size_t characters;
    
    
    /**
     * Return the actual number of matches represented by a
     * MisMatchAttemptResults. If a mask is specified, only counts matches with
     * 1s in the mask.
     */
    inline size_t getLength(const GenericBitVector* mask = NULL) {
        
        // How many total results are there?
        size_t total = 0;
        
         for(auto position : positions) {
            // Add in the length from every position.
            total += position.first.getLength(mask);
         }
         
         // Return the sum.
         return total;
    }
    
    /**
     * Returns true if there are no for the search (when restricting ourselves
     * to the BWT positions covered by the mask), false otherwise.
     */
    inline bool isEmpty(const GenericBitVector* mask = NULL) const {
        for(auto position : positions) {
            // Go through all the ranges
            if(!position.first.isEmpty(mask)) {
                // We found a range that includes something
                return false;
            }
        }
        // No range included anything
        return true;
    }
    
    /**
     * Return the range number that all search results belong to, or -1 if there
     * is no such range number.
     *
     * May only be called if the result set is nonempty.
     */
    inline int64_t range(const GenericBitVector& ranges, 
        const GenericBitVector* mask = NULL) const {
        
        // Check the first search result (which must exist because we said so.)
        int64_t candidate = positions[0].first.range(ranges, mask);
        
        if(candidate == -1) {
            // Already failed.
            return candidate;
        }
        
        for(size_t i = 1; i < positions.size(); i++) {
            // Check to make sure everything else has the same range as the
            // first result.
            if(positions[i].first.range(ranges, mask) != candidate) {
                // They don't all match
                return -1;
            }
        }
        
        // Everybody has the same range.
        return candidate;
    }
    
    
};

#endif
