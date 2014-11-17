#ifndef MISMATCHRESULT_HPP
#define MISMATCHRESULT_HPP

#include "FMDPosition.hpp"
#include <deque>

// We use a reference to the index to do extensions.
class FMDIndex;

/**
 * Represents a single BWT range in a search that allows mismatches. Keeps track
 * of the positions of its mismatches, and knows how to extend on the left and
 * retract on the right.
 */ 
class MismatchResult {
public:

    // Default copy constructor is OK, but assignment won't work because we
    // contain a reference.
    
    /**
     * Make a new MismatchResult covering an entire index.
     */
    MismatchResult(const FMDIndex& index);
    
    /**
     * Return true if this MismatchResult is the same as the other one. Needed
     * so we can make an std::set of these.
     */
    inline bool operator==(const MismatchResult& other) const {
        
        return result.getForwardStart() == other.result.getForwardStart() &&
            // Note that we ignore the reverse start on purpose.
            result.getEndOffset() == other.result.getEndOffset() &&
            startIndex == other.startIndex && endIndex == other.endIndex &&
            // We can just use the deque's built-in equality test.
            mismatches == other.mismatches;
    }
    
    /**
     * Return true if this MismatchResult is the same as the other one. Needed
     * for completeness.
     */
    inline bool operator!=(const MismatchResult& other) const {
        
        return !(*this == other);
    }
    
    /**
     * Return true if this MismatchResult belongs before the other one when
     * sorting. Needed so we can make an std::set of these.
     */
    inline bool operator<(const MismatchResult& other) const {
        return result.getForwardStart() < other.result.getForwardStart() ||
            // Note that we ignore the reverse start on purpose.
            result.getEndOffset() < other.result.getEndOffset() ||
            startIndex < other.startIndex || endIndex < other.endIndex ||
            // We can also use the deque's built-in less than operator.
            mismatches < other.mismatches;
    }
    
    /**
     * Return true if this MismatchResult belongs before the other one when
     * sorting. Needed for completeness.
     */
    inline bool operator>(const MismatchResult& other) const {
        return result.getForwardStart() > other.result.getForwardStart() ||
            // Note that we ignore the reverse start on purpose.
            result.getEndOffset() > other.result.getEndOffset() ||
            startIndex > other.startIndex || endIndex > other.endIndex ||
            // We can also use the deque's built-in greater than operator.
            mismatches > other.mismatches;
    }
    
    /**
     * Return the number of characters currently searched in this result.
     */
    inline size_t getCharacters() const {
        return endIndex - startIndex;
    }
    
    /**
     * Return the number of mismatches in this result.
     */
    inline size_t getMismatches() const {
        return mismatches.size();
    }
    
    /**
     * Return true if the leftmost character searched is a match, and false
     * otherwise (or if there is no leftmost character).
     */
    inline bool leftIsMatch() const {
        return (getCharacters() > 0) && ((getMismatches() == 0) ||
            // Remember: new (lefter) characters go at the back.
            mismatches.back() != endIndex); 
    }
    
    /**
     * Return true if the rightmost character searched is a match, and false
     * otherwise (or if there is no rightmost character).
     */
    inline bool rightIsMatch() const {
        return (getCharacters() > 0) && ((getMismatches() == 0) ||
            // Remember: old (righter) characters come out the front.
            mismatches.front() != startIndex); 
    }
    
    /**
     * Return true if we have no results under the given mask.
     */
    inline bool isEmpty(const GenericBitVector* mask = NULL) const {
        return result.isEmpty(mask);
    }
    
    /**
     * Return the number of positions actually found, under the given mask.
     */
    inline size_t getLength(const GenericBitVector* mask = NULL) const {
        // Remember that this call can return a negative length if you get into
        // a strange state (TODO: when?).
        int64_t length = result.getLength(mask);
        
        return length < 0 ? 0 : length;
    }
    
    /**
     * Return the range number that the result belongs to, or -1 if there is no
     * such range.
     */
    inline int64_t range(const GenericBitVector& ranges, 
        const GenericBitVector* mask = NULL) const {
        
        // Delegate to the FMDPosition.
        return result.range(ranges, mask);
    }
    
    /**
     * Extend this result left with every possible character, and discard
     * extensions that are completely empty. Gives a maximum of four results.
     * (Some may still be empty under whatever BWT mask you happen to be using.)
     * If the caller wants to look at only exact matches or limit total mismatch
     * count, the caller has to drop results (or not bother extending those at
     * the mismatch limit already).
     */
    std::vector<MismatchResult> extendLeft(const FMDIndex& index, 
        char match) const;
        
    /**
     * Retract the rightmost character from this result, producing a new result.
     */
    MismatchResult retractRight(const FMDIndex& index) const;
    

protected:
    /**
     * What BWT positions are our results? Only the forward interval is used.
     */
    FMDPosition result;
    
    /**
     * Keep a queue of mismatch positions. Numbers in here are the number of the
     * extension at which the mismatch occurred, and are inserted in extension
     * order.
     */
    std::deque<size_t> mismatches = std::deque<size_t>();
    
    /**
     * Keep track of what the earliest extension that hasn't been retracted is.
     * This lets us know when to pick up a mismatch.
     */
    size_t startIndex = 0;
    
    /**
     * Keep track of what the latest extension is is, so we can add mismatches
     * there in our queue, to be picked up on the appropriate retract.
     */
    size_t endIndex = 0;
};

#endif
