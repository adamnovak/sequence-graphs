#ifndef FMDPOSITION_HPP
#define FMDPOSITION_HPP

#include <iostream>
#include <algorithm>
#include <utility>

#include "BitVector.hpp"
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
     */
    inline bool isEmpty(BitVectorIterator* mask = NULL) const {
        return getLength(mask) <= 0;
    }

    /**
    * Return the actual number of matches represented by an FMDPosition. If a
    * mask is specified, only counts matches with 1s in the mask.
    */
    inline int64_t getLength(BitVectorIterator* mask = NULL) const
    {
        if(mask == NULL || end_offset == -1) {
            // Fast path: no mask or an actually empty interval. Can just look
            // at our end offset.
            return end_offset + 1;
        } else {
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
        }
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
    int64_t range(BitVectorIterator& ranges, BitVectorIterator* mask = NULL)
        const;

    /**
     * Return the number of ranges that the forward-strand interval of this
     * FMDPosition overlaps. If a mask is specified, only counts matches with 1s
     * in the mask.
     */
    int64_t ranges(BitVectorIterator& ranges,BitVectorIterator* mask = NULL)
        const;
    
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
