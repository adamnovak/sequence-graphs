#ifndef FMDPOSITION_HPP
#define FMDPOSITION_HPP

#include <iostream>

#include "RangeVector.hpp"

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
    FMDPosition();
    FMDPosition(usint forward_start, usint reverse_start, usint end_offset);
    
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
     * Are two FMDPositions equal?
     */
    bool operator==(const FMDPosition& other) const;

    /**
     * Is an FMDPosition empty?
     */
    inline bool isEmpty() const
    {
        return end_offset < 0;
    }

    /**
    * Return the actual number of matches represented by an FMDPosition.
    */
    inline int64_t getLength() const
    {
        return end_offset + 1;
    }

    /**
     * Return the index of the range that the forward-strand interval of this
     * FMDPosition is contained in, or -1 if it is not contained in any such
     * interval.
     *
     * Note that empty intervals (where the end is before the start) may still
     * be contained in ranges.
     */
    int64_t range(const RangeVector& ranges) const;

    /**
     * Return the number of ranges that the forward-strand interval of this
     * FMDPosition overlaps.
     */
    int64_t ranges(const RangeVector& ranges) const;
    
protected:
    int64_t forward_start;
    int64_t reverse_start;
    // Offset 0 = only the entry at start/end. -1 = empty.
    int64_t end_offset;
  
};

/**
 * Provide pretty-printing for FMDPositions. See
 * <http://www.parashift.com/c++-faq/output-operator.html>
 */
std::ostream& operator<< (std::ostream& o, FMDPosition const& position);

const FMDPosition EMPTY_FMD_POSITION = FMDPosition(0, 0, -1);

#endif
