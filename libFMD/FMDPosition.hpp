#ifndef FMDPOSITION_HPP
#define FMDPOSITION_HPP

#include <iostream>
#include <algorithm>
#include <utility>

// Forward declaration because we have a bit of a circular dependency with us
// and the index and the view.
class FMDPosition;

#include "FMDIndexView.hpp"
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
     * Is an FMDPosition empty under the specified view?
     */
    inline bool isEmpty(const FMDIndexView& view) {
        return getLength(view) <= 0;
    }

    /**
     * Does this FMDPosition indicate a unique merged position under the gievn
     * view? This relies on the view not having ranges defined with nothing
     * masked in, and thus can avoid needing to check if multiple merged ranges
     * belong to the same position.
     */
    inline bool isUnique(const FMDIndexView& view) {
    
        if(isEmpty(view)) {
            // If it's empty, it's not unique.
            return false;
        }
        
        // Return true if we can identify a range number, false otherwise (in
        // which case we span multiple ranges).
        return range(view) != -1;
    
    }
    
    /**
     * Does an FMDPosition select multiple masked-in positions or merged ranges
     * under the given view?
     */
    inline bool isAmbiguous(const FMDIndexView& view) const {
        
        // If it's not empty and it's not unique, it must have multiple things
        // in it.
        return !isEmpty(view) && !isUnique(view);
    }
    
    /**
     * Given that this FMDPosition is unique, get the TextPosition it uniquely
     * maps to.
     */
    inline TextPosition getTextPosition(const FMDIndexView& view) const {
        // Since we are unique, we know range is not -1. Grab it.
        int64_t rangeNumber = getRangeNumber(view);
        
        if(view.getPositions().count(rangeNumber)) == 0) {
            // This range has no assigned position, so it must belong to the
            // TextPosition you get if you locate its first BWT position (i.e.
            // the one at that 1).
            
            if(view.getRanges() == nullptr) {
                // But ranges aren't even merged, so the range number is just a
                // BWT index. We can just locate it.
                return view.getIndex().locate(rangeNumber);
            } else {
                // This is an actual range number. Find the 1 that begins this
                // range, and locate it, and use that TextPosition.
                return view.getIndex().locate(view.getRanges().select(
                    rangeNumber));
            }
            
        } else {
            // Go look up the right TextPosition for this range and use that.
            return view.getPositions()[rangeNumber];
        }
        
    }


    /**
     * Provide pretty-printing for FMDPositions. See
     * <http://www.parashift.com/c++-faq/output-operator.html>
     */
    friend std::ostream& operator<< (std::ostream& o, 
        FMDPosition const& position);

protected:
    /**
     * Return the index of the range that the forward-strand interval of this
     * FMDPosition is contained in, or -1 if it is not contained in any such
     * interval.
     */
    int64_t getRangeNumber(const FMDIndexView& view) const;

    
    int64_t forward_start;
    int64_t reverse_start;
    // Offset 0 = only the entry at start/end. -1 = empty.
    int64_t end_offset;
  
};



const FMDPosition EMPTY_FMD_POSITION = FMDPosition(0, 0, -1);

#endif
