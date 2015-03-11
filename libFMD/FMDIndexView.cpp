#include "FMDIndexView.hpp"

FMDIndexView::FMDIndexView(const FMDIndex& index, const GenericBitVector* mask,
    const GenericBitVector* ranges,
    const std::map<size_t, TextPosition> positions): index(index), mask(mask),
    ranges(ranges), positions(std::move(positions)) {
    
    
    // Nothing to do!
}

int64_t FMDIndexView::getRangeNumber(const FMDPosition& position) const {

    // What's the first index in the BWT that the range would need to cover?
    size_t interval_start;
    // What's the last index in the BWT that the range would need to cover?
    size_t interval_end;

    if(getMask() != nullptr) {
        // Slow path. We need to only count positions with 1s in the mask.
        // Like this:
        //  Mask:                               1  1
        //  Ranges:            |            |            |         |
        //  Position Ends:   ^                                        ^
        //  Answer:                             This one
        // Or this:
        //  Mask:                               1  
        //  Ranges:            |            |            |         |
        //  Position Ends:   ^                                        ^
        //  Answer:                             This one
        // Or this:
        //  Mask:                               1  1
        //  Ranges:            |            |            |         |
        //  Position Ends:                    ^       ^                    
        //  Answer:                             This one
        
        // Eevn though we now have the constraint that each range contain at
        // least one 1 in the mask, this still works.
        
        // Algorithm:
        // Find the first 1 in the mask in our interval.
        // Find the last 1 in the mask in our interval.
        // See what range spans them, if any.
        
        // Find the first 1 at position >= our start. This holds (index, rank-1)
        auto firstOne = getMask()->valueAfter(position.getForwardStart());
        // TODO: What if there is no such 1?
        
        // Find the last 1 at position <= our end. This holds (index, rank-1),
        // or (size, items) if no such 1 exists.
        auto lastOne = getMask()->valueBefore(position.getForwardStart() +
            position.getEndOffset());
        if(lastOne.first > position.getForwardStart() +
            position.getEndOffset()) {
            
            // We didn't find such a 1. There must be no 1 in our interval and
            // we are thus empty.
            return -1;
        }
        
        // Fill in the specs for the interval we want.
        interval_start = firstOne.first;
        interval_end = lastOne.first;
    } else {
        // Fast path; no need to mess about with the mask. Just look at the
        // first and last indices in the BWT interval.
        interval_start = position.getForwardStart();
        interval_end = position.getForwardStart() + position.getEndOffset();   
    }

    if(interval_end < interval_start || position.getEndOffset() < 0) {
        // Our range is actually empty.
        // TODO: Do we need this check?
        return -1;
    }

    // If we get here we actually have a pair of endpoints (both inclusive), and
    // we want to know if they are in the same range.
    
    if(getRanges() != nullptr) {
        // Merged ranges are actually defined.
    
        // Look up the range that the forward starting position is in
        int64_t start_range = getRanges()->rank(interval_start);

        // And the range the forward end is in
        int64_t end_range = getRanges()->rank(interval_end);
        
        if(start_range == end_range) {
            // Both ends of the interval are in the same range.
            return start_range;
        } else {
            // There is no range spanning our forward-strand interval.
            return -1;
        }
    } else {
        // No positions are merged
        if(interval_end == interval_start) {
            // Luckily we only have one position selected, so return that.
            return interval_start;
        } else {
            // We have multiple positions selected, and they can't be merged
            return -1;
        }
    }
}
