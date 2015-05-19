#ifndef FMDPOSITION_HPP
#define FMDPOSITION_HPP

#include <iostream>
#include <algorithm>
#include <utility>
#include <set>

#include "GenericBitVector.hpp"
#include "TextPosition.hpp"
#include "Log.hpp"

// Forward declaration of FMDPosition
class FMDIndexView;

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
     * How many BWT indices are selected by this range (ignoring any sort of
     * view, mask, etc.)?
     *
     * Minimum length is 0.
     */
    inline size_t getLength() const {
        // Make sure we don't try and return a negative number as a size_t, by
        // maxing against 0 while signed.
        return std::max(end_offset + (int64_t) 1, (int64_t) 0);        
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
     * Is one FMDPosition less than another (so they can go in a set)?
     */
    bool operator<(const FMDPosition& other) const;


    /**
     * Provide pretty-printing for FMDPositions. See
     * <http://www.parashift.com/c++-faq/output-operator.html>
     */
    friend std::ostream& operator<< (std::ostream& o, 
        FMDPosition const& position);

    ////////////////////////////////////////////////////////////////////////////
    // Operations under a view.
    ////////////////////////////////////////////////////////////////////////////

    /**
     * Extend left with no mismatches.
     */
    void extendLeftOnly(const FMDIndexView& view, char character);
    
    /**
     * Retract to the given search string length, on the right.
     */
    void retractRightOnly(const FMDIndexView& view, size_t newLength);
    
    /**
     * Retract until more BWT positions are selected. They may not actually be
     * different merged ranges or masked in.
     */
    size_t retractRightOnly(const FMDIndexView& view);

    /**
     * Is nothing selected under the given view?
     */
    bool isEmpty(const FMDIndexView& view) const;
    
    /**
     * Is exactly one merged position selected under the given view?
     */
    bool isUnique(const FMDIndexView& view) const;
    
    /**
     * Is more than one merged position selected under the given view?
     */
    bool isAmbiguous(const FMDIndexView& view) const;
    
    /**
     * Get the unique TextPosition selected under the viven view. isUnique()
     * must be true.
     */
    TextPosition getTextPosition(const FMDIndexView& view) const;
    
    /**
     * Provides an overestimate of the number of ranges selected under the given
     * view.
     */
    size_t getApproximateNumberOfRanges(const FMDIndexView& view) const;
    
    /**
     * Provides an overestimate of the number of new ranges selected under the
     * given view, relative to the given old FMDPosition.
     */
    size_t getApproximateNumberOfNewRanges(const FMDIndexView& view,
        const FMDPosition& old);
        
    /**
     * Get the TextPositions selected under the given view.
     */
    std::set<TextPosition> getTextPositions(const FMDIndexView& view) const;
    
    /**
     * Get the new TextPositions selected under the given view, relative to the
     * given old FMDPosition.
     */
    std::set<TextPosition> getNewTextPositions(const FMDIndexView& view,
        const FMDPosition& old) const;

protected:
    int64_t forward_start;
    int64_t reverse_start;
    // Offset 0 = only the entry at start/end. -1 = empty.
    int64_t end_offset;
  
};



const FMDPosition EMPTY_FMD_POSITION = FMDPosition(0, 0, -1);

#endif
