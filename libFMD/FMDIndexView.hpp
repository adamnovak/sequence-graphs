#ifndef FMDINDEXVIEW_HPP
#define FMDINDEXVIEW_HPP

#include "FMDIndex.hpp"
#include "GenericBitVector.hpp"
#include "TextPosition.hpp"

#include <vector>

/**
 * Represents an FMDIndex taken together with a graph structure merging
 * positions, and a mask including some positions and excluding others.
 */
class FMDIndexView {
public:
    /**
     * Make a new FMDIndexView for the given FMDIndex, with the given mask, the
     * given bitvector of merged ranges, and the given map from range number to
     * assigned position.
     *
     * All pointers must be to objects that will outlive this FMDIndexView. No
     * ownership is taken.
     */
    FMDIndexView(const FMDIndex& index, const GenericBitVector* mask = nullptr,
        const GenericBitVector* ranges = nullptr,
        const std::map<size_t, TextPosition> positions =
        std::map<size_t, TextPosition>());
    
    /**
     * Get the FMDIndex this is a view of.
     */    
    const inline FMDIndex& getIndex() const {
        return index;
    }
    
    /**
     * Get the mask that this view uses.
     */
    const inline GenericBitVector* getMask() const {
        return mask;
    }
    
    /**
     * Get the ranges bitvector that this view uses to specifty the start points
     * of merged ranges.
     */
    const inline GenericBitVector* getRanges() const {
        return ranges;
    }
    
    /**
     * Get a reference to the range-number-to-position-and-orientation map that
     * this view uses.
     */
    const inline std::map<size_t, TextPosition>& getPositions() const {
        return positions;    
    }
    
    /***************************************************************************
     * Functions for working with FMDPositions
     **************************************************************************/
    
    /**
     * Is an FMDPosition empty under this view?
     */
    inline bool isEmpty(const FMDPosition& position) const {
        if(position.getEndOffset() < 0) {
            // It has no BWT positions in it
            return true;
        } else if(getMask() != nullptr) {
            // It has some BWT positions in it, and we have a mask
            
            // Where's the last position that's in the BWT range? TODO: we use
            // this pattern so much it should be a method on FMDPosition.
            size_t endInclusive = position.getForwardStart() + 
                position.getEndOffset();
                
            // Return true if we *don't* have any 1s in the mask between the
            // interval endpoints.
            return (getMask()->rank(endInclusive) + 1 - getMask()->rank(
                position.getForwardStart(), true) <= 0);
        } else {
            // It contains BWT positions and we have no mask, so they are all
            // masked in.
            return false;
        }
    }

    /**
     * Does this FMDPosition indicate a unique merged position under this view?
     * This relies on the view not having ranges defined with nothing masked in,
     * and thus can avoid needing to check if multiple merged ranges belong to
     * the same position.
     */
    inline bool isUnique(const FMDPosition& position) const {
    
        if(isEmpty(position)) {
            // If it's empty, it's not unique.
            return false;
        }
        
        // Return true if we can identify a range number, false otherwise (in
        // which case we span multiple ranges).
        return getRangeNumber(position) != -1;
    
    }
    
    /**
     * Does an FMDPosition select multiple masked-in positions or merged ranges
     * under this view?
     */
    inline bool isAmbiguous(const FMDPosition& position) const {
        
        // If it's not empty and it's not unique, it must have multiple things
        // in it.
        return !isEmpty(position) && !isUnique(position);
    }
    
    /**
     * Convert a range number (which may actually be a BWT position if the view
     * does not use a range vector) to a TextPosition.
     *
     * There must be at least one masked-in BWT position in the range.
     *
     * TODO: Make this protected? Would anyone else ever want it?
     */
    inline TextPosition rangeToTextPosition(size_t rangeNumber) const {
        
        if(getPositions().count(rangeNumber) == 0) {
            // This range has no assigned position, so it must belong to the
            // TextPosition you get if you locate its first masked in BWT
            // position.
            
            if(getRanges() == nullptr) {
                // But ranges aren't even merged, so the range number is just a
                // BWT index. We can just locate it.
                return getIndex().locate(rangeNumber);
            } else {
                // This is an actual range number. Find the 1 that begins this
                // range, get the first masked in position after that (if
                // applicable) and locate it, and use that TextPosition.
                
                // Where does this range start?
                auto bwtIndex = getRanges()->select(rangeNumber);
                
                if(getMask() != nullptr) {
                    // Make sure we select a masked-in position.
                    bwtIndex = getMask()->valueAfter(bwtIndex).first;
                    
                    // TODO: make sure this obeys all the preconditions on the
                    // bitvector.
                    
                    if(bwtIndex >= getRanges()->select(rangeNumber + 1)) {
                        // Complain if the next masked-in position is not in
                        // this range.
                        throw std::runtime_error("Tried to get a "
                            "TextPosition for a range with nothing in it!");
                    }
                }
                
                // Find the TextPosition for the selected BWT position.
                auto pos = getIndex().locate(bwtIndex);
                
                Log::debug() << "Range " << rangeNumber <<
                    " is TextPosition " << pos << std::endl;
                return pos;
            }
            
        } else {
            // Go look up the right TextPosition for this range and use that.
            return getPositions().at(rangeNumber);
        }
        
    }
    
    /**
     * Given that this FMDPosition is unique, get the TextPosition it uniquely
     * maps to.
     */
    inline TextPosition getTextPosition(const FMDPosition& position) const {
        
        // Since we are unique, we know range is not -1. Grab it and immediately
        // make it a TextPosition.
        return rangeToTextPosition(getRangeNumber(position));
        
    }
    
    /**
     * Get all of the TextPositions that this FMDPosition selects. Automatically
     * de-duplicates those that might be included through multiple ranges.
     */
    inline std::set<TextPosition> getTextPositions(
        const FMDPosition& position) const {
        
        // Get all the range numbers we have masked-in positions in.
        std::vector<size_t> rangeNumbers = getRangeNumbers(position);
        
        // We'll convert the range numbers to text positions and populate this.
        std::set<TextPosition> toReturn;
        
        for(const auto& rangeNumber : rangeNumbers) {
            // Map the range-number-to-text-position over the range numbers,
            // putting the results in the set.
            toReturn.insert(rangeToTextPosition(rangeNumber));
        }
        
        // Return the results.
        return toReturn;
    }
    
    /**
     * Find TextPositions for all the BWT positions which were not selected in
     * the old FMDPosition but which are selected in the wider one. Note that
     * some of the returned TextPositions may be ones that were already selected
     * in the old range, if new BWT positions merged into the same TextPositions
     * are selected.
     *
     * The wider FMDPosition must represent a range containing that of the old
     * FMDPosition.
     */
    inline std::set<TextPosition> getNewTextPositions(
        const FMDPosition& old, const FMDPosition& wider) const {
        
        // Get all the range numbers we found new stuff in.
        std::vector<size_t> rangeNumbers = getNewRangeNumbers(old, wider);
        
        // We'll convert the range numbers to text positions and populate this.
        std::set<TextPosition> toReturn;
        
        for(const auto& rangeNumber : rangeNumbers) {
            // Map the range-number-to-text-position over the range numbers,
            // putting the results in the set.
            toReturn.insert(rangeToTextPosition(rangeNumber));
        }
        
        // Return the results.
        return toReturn;
        
    }
    
    /**
     * Find (approximately) the number of merged ranges selected by an
     * FMDPosition. Provides an overestimate of the number of items in the set
     * getTextPositions() will return.
     */
    size_t getApproximateNumberOfRanges(const FMDPosition& position) const;
        
    /**
     * Find (approximately) the number of merged ranges selected by a new
     * FMDPosition over an old one. Provides an overestimate of the number of
     * items in the set getNewTextPositions() will return.
     */
    size_t getApproximateNumberOfNewRanges(const FMDPosition& old,
        const FMDPosition& wider) const;
    
protected:
    /**
     * Return the index of the range that the forward-strand interval of this
     * FMDPosition is contained in, or -1 if it is not contained in any such
     * interval.
     */
    int64_t getRangeNumber(const FMDPosition& position) const;
    
    /**
     * Return the indices of all of the ranges that the forward-strand interval
     * of this FMDPosition has masked-in positions in.
     */
    std::vector<size_t> getRangeNumbers(const FMDPosition& position) const;
    
    /**
     * Return the indices of all of the ranges that the forward-strand interval
     * of the wider FMDPosition has masked-in positions in, in the part of it
     * not overlapped by the old interval.
     *
     * The wider FMDPosition must represent a range containing that of the old
     * FMDPosition.
     */
    std::vector<size_t> getNewRangeNumbers(const FMDPosition& old,
        const FMDPosition& wider) const;

    /**
     * What FMDIndex are we a view of? It must of course outlive us.
     */
    const FMDIndex& index;
    
    /**
     * What positions from that index are included in the view? If this is null,
     * all positions are included. Otherwise only positions with 1s are
     * included.
     */
    const GenericBitVector* mask;
    
    /**
     * What are the boundaries of merged ranges that exist? If this is null,
     * each position is its own unique range, not merged with any other
     * position. Otherwise, 1s in this bitvector mark the starts of merged
     * ranges.
     *
     * This bitvector may onyl have a 1 for a range if at least one position in
     * that range is masked in under the mask. Otherwise, it would be possible
     * to have a BWT range that spans multiple ranges here, but which ought to
     * be unique because the intervening range had everything masked out.
     */
    const GenericBitVector* ranges;
    
    /**
     * What is the position and orientation that each merged range belongs to?
     * This map stores, by range number, a TextPosition representing the
     * position that each merged range belongs to, and the orientation of the
     * position that the range is associated with. If ranges is null, this will
     * be empty. If an entry does not exist for a range, that range belongs to
     * the position and orientation associated with its lowest BWT position.
     * (Thus, single-item ranges do not need entries here.)
     *
     * This is the structure that ties multiple merged ranges that ought to
     * belong to the same position together. It also ties forward and reverse-
     * complement ranges together.
     *
     * If ranges is null, this must be empty.
     */
    const std::map<size_t, TextPosition> positions;
};

#endif
