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

private:
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
