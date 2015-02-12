#ifndef FMDPOSITION_HPP
#define FMDPOSITION_HPP

#include <iostream>
#include <algorithm>
#include <utility>

#include "GenericBitVector.hpp"
#include "Log.hpp"

/**
 * Represents the state (or result) of an FMD-index search, which is two ranges
 * (one for the forward sequence, and one for the reverse complement) of equal
 * length in the BWT. The ranges are stored as two start indices and a length.
 *
 * Range semantics are inclusive, so a length = 0 range holds 1 thing and its
 * reverse complement.
 */
class FMDPosition
{
public:
    /**
     * Create a new FMDPosition with undefined contents.
     */
    FMDPosition();
    
    /**
     * Create a new FMDPosition representing two ranges starting at the given
     * positions and extending for the given offset (where an offset of 0
     * indicates two 1-character ranges).
     */
    FMDPosition(int64_t forward_start, int64_t reverse_start, 
        int64_t end_offset);
    
    /**
     * Get the forward-strand start position of the range.
     */
    inline int64_t getForwardStart() const {
        return forward_start;
    }
    
    /**
     * Set the forward-strand start position.
     */
    inline void setForwardStart(int64_t value) {
        forward_start = value;
    }
    
    /**
     * Get the reverse-strand start position of the range.
     */
    inline int64_t getReverseStart() const {
        return reverse_start;
    }
    
    /**
     * Set the reverse-strand start position.
     */
    inline void setReverseStart(int64_t value) {
        reverse_start = value;
    }
    
    /**
     * Get the offsets of the range ends from the range starts. 0 means a
     * 1-element range on each strand.
     */
    inline int64_t getEndOffset() const {
        return end_offset;
    }
    
    /**
     * Set the interval's end offset. 0 means a 1-element range on each strand.
     */
    inline void setEndOffset(int64_t value) {
        end_offset = value;
    }
    
    
    
    /** 
     * Flip the FMDPosition around so the reverse complement interval is the
     * forward interval and visa versa.
     */
    FMDPosition flip() const;
    
    /**
     * Flip the FMDPosition around in place, modifying it so the reverse
     * complement interval is the forward interval and visa versa.
     */
    inline void flipInPlace() {
        std::swap(forward_start, reverse_start);
    }

    /**
     * Are two FMDPositions equal?
     */
    bool operator==(const FMDPosition& other) const;

    /**
     * Is an FMDPosition empty? If a mask is specified, only counts matches with
     * 1s in the mask.
     *
     * Does not depend on ranges.
     */
    inline bool isEmpty(const GenericBitVector* mask = NULL) const {
        // This is guaranteed to take the not-super-slow path, since no ranges
        // vector is specified.
        return getLength(mask) <= 0;
    }
    
    /**
     * Does an FMDPosition include exactly one masked-in position or range of
     * positions?
     *
     * If mask is null, assumes all positions are masked in.
     *
     * If ranges is null, assumes no positions are merged.
     */
    inline bool isUnique(const GenericBitVector* mask = NULL,
        const GenericBitVector* ranges = NULL) const {
        
        if(isEmpty(mask)) {
            // If it's empty, it's not unique.
            return false;
        }
        
        // Return true if we can identify a range number, false otherwise (in
        // which case we span multiple ranges).
        return range(*ranges, mask) != -1;
        
    }
    
    /**
     * Does an FMDPosition select multiple masked-in positions or merged ranges?
     *
     * If mask is null, assumes all positions are masked in.
     *
     * If ranges is null, assumes no positions are merged.
     */
    inline bool isAmbiguous(const GenericBitVector* mask = NULL,
        const GenericBitVector* ranges = NULL) const {
        
        // If it's not empty and it's not unique, it must have multiple things
        // in it.
        return !isEmpty(mask) && !isUnique(mask, ranges);
    }

    /**
     * Return the actual number of matches represented by an FMDPosition. If a
     * mask is specified, only counts matches with 1s in the mask.
     *
     * If ranges is not null, counts each range marked by a leading 1 in that
     * vector as a single position.
     */
    inline int64_t getLength(const GenericBitVector* mask = NULL,
        const GenericBitVector* ranges = NULL) const {
        
        if((mask == NULL && ranges == NULL) || end_offset == -1) {
            // Fast path: no mask or ranges, or an actually empty interval. Can
            // just look at our end offset.
            return end_offset + 1;
        } else if(ranges == NULL) {
            // Slow path: need to make rank queries, but we know the interval is
            // nonempty. Get the rank at the end of the region (inclusive), and
            // subtract the rank at the beginning of the region (exclusive). We
            // need a +1 since we actually measure 1 the inclusive rank of the
            // previous position and need to get rid of the extra 1.
            
            Log::trace() << "Mask rank at " << forward_start + end_offset <<
                ": " << mask->rank(forward_start + end_offset) << std::endl;
                
            Log::trace() << "Mask rank at least at " << forward_start <<
                ": " << mask->rank(forward_start, true) << std::endl;
            
            return mask->rank(forward_start + end_offset) + 1 - 
                mask->rank(forward_start, true);
        } else if(mask == NULL) {
            // Slow path: need to make rank queries in the ranges, but we don't
            // have to deal with a mask.
            
            // The answer is 1 more than the number of 1s between the start and
            // the past-the-end position.
            return ranges->rank(forward_start + end_offset) - 
                ranges->rank(forward_start) + 1;
        } else {
            // Slowest path: need to deal with both a mask and ranges.
            
            // This will just be the number of masked-in ranges.
            // TODO: are we basically replicating that method here?
            return this->ranges(*ranges, mask);
        }
    }
    
    /**
     * Get the BWT index of each selected, masked-in position.
     */
    inline std::vector<int64_t> getResults(
        const GenericBitVector* mask = NULL) const {
        // TODO: Implement based on rank and select. For now we just scan.
        
        std::vector<int64_t> toReturn;
        
        for(size_t i = getForwardStart();
            i <= getForwardStart() + getEndOffset(); i++) {
        
            // For each base in the range
            
            if(!mask || mask->isSet(i)) {
                // If there is no mask or we're in it, take this position.
                toReturn.push_back(i);
            }
        }
        
        return toReturn;
    }
    
    /**
     * Get a single BWT position that is selected and masked in. Assumes such a
     * position exists.
     */
    inline int64_t getResult(const GenericBitVector* mask = NULL) const {
        // TODO: Implement based on rank and select. For now we just scan.
        
        std::vector<int64_t> toReturn;
        
        for(size_t i = getForwardStart();
            i <= getForwardStart() + getEndOffset(); i++) {
        
            // For each base in the range
            
            if(!mask || mask->isSet(i)) {
                // If there is no mask or we're in it, take this position.
                return i;
            }
        }
        
        throw std::runtime_error("Attempted to get non-existent result.");
    }

    /**
     * Return the index of the range that the forward-strand interval of this
     * FMDPosition is contained in, or -1 if it is not contained in any such
     * interval.
     *
     * Note that empty intervals (where the end is before the start) may still
     * be contained in ranges.
     *
     * If a mask is specified, only counts matches with 1s in the mask.
     */
    int64_t range(const GenericBitVector& ranges, 
        const GenericBitVector* mask = NULL) const;

    /**
     * Return the number of ranges that the forward-strand interval of this
     * FMDPosition overlaps. If a mask is specified, only counts matches with 1s
     * in the mask.
     */
    int64_t ranges(const GenericBitVector& ranges, 
        const GenericBitVector* mask = NULL) const;
    
    /**
     * Provide pretty-printing for FMDPositions. See
     * <http://www.parashift.com/c++-faq/output-operator.html>
     */
    friend std::ostream& operator<< (std::ostream& o, 
        FMDPosition const& position);
    
protected:
    int64_t forward_start;
    int64_t reverse_start;
    // Offset 0 = only the entry at start/end. -1 = empty.
    int64_t end_offset;
  
};



const FMDPosition EMPTY_FMD_POSITION = FMDPosition(0, 0, -1);

#endif
