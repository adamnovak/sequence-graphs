#ifndef RESULTRANGE_HPP
#define RESULTRANGE_HPP

#include <queue>
#include <array>
#include "GenericBitVector.hpp"
#include "FMDIndex.hpp"
#include "util.hpp"


/**
 * Represents a range of search results, and some metadata about the search
 * string used to find them. Allows us to track and retract mismatches.
 *
 * Search proceeds from right to left through the query string, extending on the
 * left, retracting on the right, and mapping positions to strings based on
 * their left context. This means that positions with similar contexts (and
 * that are thus likely to be merged together) will be near each other in BWT
 * space, so BWT ranges make sense for merged positions.
 */
class ResultRange {
public:
    /**
     * Make a new, empty ResultRange.
     */
    ResultRange();

    /**
     * Make a ResultRange covering the whole BWT from the given FMDIndex,
     * starting at the end of the given query.
     */
    ResultRange(const FMDIndex& field, size_t queryLength);
    
    /**
     * Try extending on the left with a match. Return the extended ResultRange
     * (which may be empty). The range being extended may not be empty.
     */
    ResultRange extendLeftMatch(const FMDIndex& index,
        const std::string& query);
    
    /**
     * Try extending on the left with all 3 possible mismatches. Return all 3
     * ResultRanges produced, some of which may be empty. The range being
     * extended must not be empty.
     */
    std::array<ResultRange, 3> extendLeftMismatch(const FMDIndex& index, 
        const std::string& query);
        
    /**
     * Retract one character on the right, forgetting its mismatch if necessary.
     * The range being retracted must not be empty, and it must have a nonempty
     * search string.
     */
    ResultRange retractRight(const FMDIndex& index);
    
    /**
     * Return whether the range is empty or not under the given mask if any.
     */
    bool isEmpty(GenericBitVector* mask = NULL);
    
    /**
     * Return the number of items in the range, under the given bitvector mask
     * if any.
     */
    size_t getLength(GenericBitVector* mask = NULL);
    
    /**
     * Return the number of characters in the search string for this result
     * range.
     */
    size_t getSearchStringLength();

    /**
     * Copy is OK.
     */
    ResultRange(const ResultRange& other) = default;

    /**
     * Move is OK.
     */
    ResultRange(ResultRange&& other) = default;
    
    /**
     * Assignment is OK.
     */
    ResultRange& operator=(const ResultRange& other) = default;
    
    /**
     * Move assignment is OK.
     */
    ResultRange& operator=(ResultRange&& other) = default;
    
private:
    // For now we just wrap an FMDPosition
    FMDPosition position;
    
    // Where in the overall query string does the part we have searched start?
    // Inclusive.
    size_t searchStringStart;
    // Where in the overall query string does the part we have searched end?
    // Exclusive.
    size_t searchStringEnd;
    
    // Which positions in the search string are included as mismatches? These
    // are sorted in ascending order; new ones are added at the back and old
    // ones are popped off at the front.
    std::queue<size_t> mismatches;
    

};

#endif
