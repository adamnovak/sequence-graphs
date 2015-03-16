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
        int64_t start_range = getRanges()->rank(interval_start) - 1;

        // And the range the forward end is in
        int64_t end_range = getRanges()->rank(interval_end) - 1;
        
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

std::vector<size_t> FMDIndexView::getRangeNumbers(
    const FMDPosition& position) const {
    
    // TODO: make this use callbacks or something. TODO: Be able to go from
    // TextPosition to reverse complement to range number somehow?

    // We'll populate this with the number of every range for which we overlap a
    // 1 in the mask, each of which will appear once.
    std::vector<size_t> toReturn;

    if(position.getEndOffset() < 0) {
        // Special case the empty FMDPosition so we don't have to deal with that
        // case later.
        
        // Just skip all this stuff and end up returning our empty vector.
        
    } else if(getMask() == nullptr && getRanges() == nullptr) {
        // No mask and no ranges to deal with. Return the position number for
        // each selected position.
        
        for(size_t i = position.getForwardStart(); 
            i <= position.getForwardStart() + position.getEndOffset(); i++) {
            // Make an entry for each of the BWT positions
            toReturn.push_back(i);
        }
        
    } else if(getMask() == nullptr && getRanges() != nullptr) {
        // Fast path. Since there's no mask, every range counts. But we do have
        // to go get the ranges instead of just using position numbers.
    
        // Look up the range that the starting position is in. We have to
        // subtract 1 because the 0th range begins with a 1.
        int64_t startRange = getRanges()->rank(position.getForwardStart()) - 1;

        // And the range the end is in
        int64_t endRange = getRanges()->rank(position.getForwardStart() +
            position.getEndOffset()) - 1;
            
        Log::debug() << "Looking for ranges between " <<
            position.getForwardStart() << " in range " << startRange << 
            " which starts at " << getRanges()->select(startRange) <<
            " and " << position.getForwardStart() + position.getEndOffset() <<
            " in range " << endRange << " which starts at " << 
            getRanges()->select(endRange) << std::endl;
    
        
        if(endRange < startRange) {
            throw std::runtime_error("End of ranges interval is before start");
        }
        
        for(size_t i = startRange; i <= endRange; i++) {
            // Make an entry for each of those ranges in our result vector.
            // TODO: make sure we don't have the end range as the max size_t.
            toReturn.push_back(i);
        }
    
    } else {
        // We have ranges and a mask. Not every range counts, only those that
        // overlap positions with 1s in the mask.
        
        // TODO: test this codepath!

        // Keep track of the left edge of the interval we still need to count
        // the ranges in.
        size_t left = position.getForwardStart();

        while(left <= position.getForwardStart() + position.getEndOffset()) {
            // Until the next range starts outside this FMDPosition...
        
            // Find the (index, rank) of the first 1 in the mask that's here or
            // beyond our current position.
            auto nextPosition = getMask()->valueAfter(left);
            
            if(nextPosition.first <= position.getForwardStart() +
                position.getEndOffset()) {
                
                // This masked-in position is within this FMDPosition. We should
                // count its range. Go look at what range that BWT position is
                // in and grab it. Make sure to make it 0-based.
                toReturn.push_back(getRanges()->rank(nextPosition.first) - 1);
                
                // Find the next 1 in the ranges vector after this masked-in
                // position (exclusive), indicating the start of the next range.
                auto nextRange = getRanges()->valueAfter(
                    nextPosition.first + 1);
                    
                // Look in or after that range for the next masked-in position.
                left = nextRange.first;
            } else {
                // The next masked-in position is out of us. We're done finding
                // ranges.
                break;
            }
            
        }
    }
    
    // Return the occupied ranges.
    return toReturn;
}

std::vector<size_t> FMDIndexView::getNewRangeNumbers(
    const FMDPosition& old, const FMDPosition& wider) const {
    
    // We can just make FMDPositions for the new ranges, keeping track of only
    // the forward bits.
    
    // Make sure to subtract 1 from the length, since a 0 offset means 1 base.
    FMDPosition newLeft(wider.getForwardStart(), 0, 
        old.getForwardStart() - wider.getForwardStart() - 1);
        
    FMDPosition newRight(old.getForwardStart() + old.getEndOffset() + 1, 0, 
        wider.getForwardStart() + wider.getEndOffset() - old.getForwardStart() -
        old.getEndOffset() - 1);
        
    // Get the answers on each side
    std::vector<size_t> leftRanges = getRangeNumbers(newLeft);
    
    std::vector<size_t> rightRanges = getRangeNumbers(newRight);
    
    // Add the new ranges on the right to the new ranges on the left. See
    // <http://stackoverflow.com/a/201729/402891>
    leftRanges.insert(leftRanges.end(), rightRanges.begin(), rightRanges.end());
        
    // Return them all.
    return leftRanges;
}
