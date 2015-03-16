#ifndef ZIPMAPPINGSCHEME_HPP
#define ZIPMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"
#include "Mapping.hpp"
#include "IntervalIndex.hpp"
#include "Log.hpp"
#include "Matching.hpp"

#include <iomanip>

#include <unordered_map>

/**
 * Mapping scheme supporting mapping to graphs, where you zip together maximal
 * matches looking out in both directions from each position. Finds a subset of
 * the unique substrings for each base that the natural mapping scheme finds,
 * but will never map a base that the actual natural to-graph mapping scheme
 * would call as conflicted.
 */
class ZipMappingScheme: public MappingScheme {

public:
    /**
     * Inherit the constructor.
     */
    using MappingScheme::MappingScheme;
    
    /**
     * Map the given query string according to the natural mapping scheme with
     * optional inexact credit. When a mapping is found, the callback function
     * will be called with the query base index, and the TextPosition to which
     * it maps in the forward direction.
     *
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) const override;
        
    
protected:

    /**
     * Use the inchworm algorithm to find the longest right context present in
     * the reference for each base in the query. Results are in the same order
     * as the characters in the string.
     */
    std::vector<FMDPosition> findRightContexts(const std::string& query) const;
    
    /**
     * Explore all retractions of the two FMDPositions. Return either an empty
     * set (if no explored retraction finds overlapping positions between the
     * two sides), a set with one element (if we find exactly one such
     * overlap), or a set with two (or more) elements (if we find multiple
     * overlaps).
     */
    std::set<TextPosition> exploreRetractions(const FMDPosition& left,
        size_t patternLengthLeft, const FMDPosition& right,
        size_t patternLengthRight) const;    
    
};

#endif
