#include "FMDIndexView.hpp"

FMDIndexView::FMDIndexView(const FMDIndex& index, const GenericBitVector* mask,
    const GenericBitVector* ranges,
    const std::map<size_t, TextPosition> positions): index(index), mask(mask),
    ranges(ranges), positions(std::move(positions)), invertedPositions() {
    
    if(ranges != nullptr) {
        // We have to deal with ranges.
    
        for(size_t i = 0; i < ranges->rank(ranges->getSize()); i++) {
            // For each 1 in the ranges bitvector
            
            if(positions.count(i) != 0) {
                // This particular range is marked as owned.
                
                // Put in an inverted entry from the owning TextPosition to this
                // range number.
                invertedPositions.emplace(positions.at(i), i);
                
                Log::trace() << "Range " << i << " explicitly owned by " <<
                    positions.at(i) << std::endl;
                
            } else {
                // We will use the first index's TextPosition as the owner.
                
                // TODO: maybe we can use an inverted sampled suffix array or
                // something here instead of the full inverted index.
                
                // Where does this range start?
                int64_t rangeStart = ranges->select(i);
                
                if(rangeStart >= index.getBWTLength()) {
                    // Sometimes we put trailing 1s. TODO: figure out where we
                    // do that and stop it.
                    
                    // We know there are no more ranges.
                    break;
                }
                
                // And what TextPosition is there?
                TextPosition owner = index.locate(rangeStart);
                
                invertedPositions.emplace(owner, i);
                
                Log::trace() << "Range " << i << " implicitly owned by " <<
                    owner << std::endl;
            }
        }
    
        Log::debug() << "Made " << invertedPositions.size() <<
            " inverted positions entries" << std::endl;
            
    }
}

int64_t FMDIndexView::getRangeNumber(size_t start, size_t length) const {

    // What's the first index in the BWT that the range would need to cover?
    size_t interval_start;
    // What's the last index in the BWT that the range would need to cover?
    size_t interval_end;

    if(getMask() != nullptr) {
        // We know every range has at least 1 masked-in position, but we don't
        // know if we cover it.
        
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
        auto firstOne = getMask()->valueAfter(start);
        // TODO: What if there is no such 1?
        
        // Find the last 1 at position <= our end. This holds (index, rank-1),
        // or (size, items) if no such 1 exists.
        auto lastOne = getMask()->valueBefore(start + length - 1);
        if(lastOne.first > start + length - 1) {
            
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
        interval_start = start;
        interval_end = start + length - 1;  
    }

    if(interval_end < interval_start || length == 0) {
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

std::pair<size_t, size_t> FMDIndexView::getRangeByNumber(
    size_t rangeNumber) const {
    
    // Select the 1 at the start of the range
    auto start = getRanges()->select(rangeNumber);
    // Select the 1 after the end of the range, and then move back 1 to be
    // inclusive.
    auto afterEnd = getRanges()->select(rangeNumber + 1);
    
    // Make a new range with and return it.
    return std::make_pair(start, afterEnd - start);
    
}

std::vector<size_t> FMDIndexView::getRangeNumbers(size_t start,
    size_t length) const {
    
    // TODO: make this use callbacks or something. TODO: Be able to go from
    // TextPosition to reverse complement to range number somehow?

    // We'll populate this with the number of every range for which we overlap a
    // 1 in the mask, each of which will appear once.
    std::vector<size_t> toReturn;

    if(length == 0) {
        // Special case the empty interval so we don't have to deal with that
        // case later.
        
        // Just skip all this stuff and end up returning our empty vector.
        
    } else if(getMask() == nullptr && getRanges() == nullptr) {
        // No mask and no ranges to deal with. Return the position number for
        // each selected position.
        
        for(size_t i = start; i < start + length; i++) {
            // Make an entry for each of the BWT positions
            toReturn.push_back(i);
        }
        
    } else if(getMask() == nullptr && getRanges() != nullptr) {
        // Fast path. Since there's no mask, every range counts. But we do have
        // to go get the ranges instead of just using position numbers.
    
        // Look up the range that the starting position is in. We have to
        // subtract 1 because the 0th range begins with a 1.
        int64_t startRange = getRanges()->rank(start) - 1;

        // And the range the end is in
        int64_t endRange = getRanges()->rank(start + length - 1) - 1;
            
        Log::trace() << "Looking for ranges between " << start << 
            " in range " << startRange << " which starts at " << 
            getRanges()->select(startRange) << " and " << start + length - 1 <<
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
    } else if(getMask() != nullptr && getRanges() == nullptr) {
        // We have a mask but no ranges vector.
        
        // Return the position number for each selected position that's masked
        // in.
        
        // Keep track of the left edge of the interval we still need to count
        // the masked-in BWT positions in.
        size_t left = start;

        while(left < start + length) {
            // Until the next masked-in position is outside this FMDPosition...
        
            // Find the (index, rank) of the first 1 in the mask that's here or
            // beyond our current position.
            auto nextPosition = getMask()->valueAfter(left);
            
            if(nextPosition.first < start + length) {
                
                // This masked-in position is within this FMDPosition. We should
                // count it.
                toReturn.push_back(nextPosition.first);
                
                // Look after there.
                left = nextPosition.first + 1;
            } else {
                // The next masked-in position is out of us. We're done finding
                // ranges.
                break;
            }
            
        }
    
    } else if(getMask() != nullptr && getRanges() != nullptr) {
        // We have ranges and a mask. Not every range counts, only those that
        // overlap positions with 1s in the mask.
        
        // TODO: test this codepath!

        // Keep track of the left edge of the interval we still need to count
        // the ranges in.
        size_t left = start;

        while(left < start + length) {
            // Until the next range starts outside this FMDPosition...
        
            // Find the (index, rank) of the first 1 in the mask that's here or
            // beyond our current position.
            auto nextPosition = getMask()->valueAfter(left);
            
            if(nextPosition.first < start + length) {
                
                // This masked-in position is within this FMDPosition. We should
                // count its range. Go look at what range that BWT position is
                // in and grab it. Make sure to make it 0-based.
                toReturn.push_back(getRanges()->rank(nextPosition.first) - 1);
                
                // Find the next 1 in the ranges vector after this masked-in
                // position (exclusive), indicating the start of the next range.
                auto nextRange = getRanges()->valueAfter(
                    nextPosition.first + 1);
                    
                if(nextRange.first <= left) {
                    // We wrapped around because we ran out of values, but
                    // landed in the selected range again.
                    break;
                }
                    
                // Look in or after that range for the next masked-in position.
                left = nextRange.first;
            } else {
                // The next masked-in position is out of us. We're done finding
                // ranges.
                break;
            }
            
        }
    } else {
        // We broke the law of the excluded middle; we checked all the cases but
        // none of them were true.
        throw std::runtime_error(
            "Impossible combination of range vector and mask vector");
    }
    
    // Return the occupied ranges.
    return toReturn;
}

std::vector<size_t> FMDIndexView::getNewRangeNumbers(size_t oldStart,
    size_t oldLength, size_t newStart, size_t newLength) const {
    
    // We can just make the new ranges, keeping track of only the forward
    // intervals.
    
    // What's the new range on the left?
    size_t newLeftStart = newStart;
    size_t newLeftLength = oldStart - newStart;

    // And on the right?
    size_t newRightStart = oldStart + oldLength;
    size_t newRightLength = newStart + newLength - oldStart - oldLength;
        
    // Get the answers on each side
    std::vector<size_t> ranges = getRangeNumbers(newLeftStart,
        newLeftLength);
    std::vector<size_t> rightRanges = getRangeNumbers(newRightStart,
        newRightLength);
    
    // Add the new ranges on the right to the new ranges on the left. See
    // <http://stackoverflow.com/a/201729/402891>
    ranges.insert(ranges.end(), rightRanges.begin(), rightRanges.end());
        
    // Return them all.
    return ranges;
}

std::vector<size_t> FMDIndexView::textPositionToRanges(
    const TextPosition& textPosition) const {
        
    // Find the bounding iterators of the range of range numbers belonging to
    // this position.
    auto bounds = invertedPositions.equal_range(textPosition);
    
    // We're going to fill up this vector with the results.
    std::vector<size_t> toReturn;
    
    for(auto kvIterator = bounds.first; kvIterator != bounds.second;
        ++kvIterator) {
    
        // For every TextPosition, range number pair that has the key we asked
        // for, put the range number in the vector.
        toReturn.push_back((*kvIterator).second);
    }
    
    // Return the result.
    return toReturn;
}

size_t FMDIndexView::getApproximateNumberOfRanges(
    size_t start, size_t length) const {

    if(getRanges() != nullptr) {
        // We have a ranges vector. Don't worry about the mask, just provide an
        // over-estimate on the number of ranges we may have selected. We know
        // at least one mask 1 is in every range, so this can't be off by too
        // much because of that.
        
        // Look up the range that the starting position is in. We have to
        // subtract 1 because the 0th range begins with a 1.
        int64_t startRange = getRanges()->rank(start) - 1;

        // And the range the end is in
        int64_t endRange = getRanges()->rank(start + length - 1) - 1;
            
        // Return the total number of ranges we overlap
        return endRange - startRange + 1;
        
    } else {
        // Ranges is null. Just count the positions, since each may be a
        // different merged range.
        
        if(getMask() != nullptr) {
            // Count only masked-in positions
            
            // How many masked-in positions exist before we start?
            int64_t startMask = getMask()->rank(start);
            
            // And how many exist before we end?            
            int64_t endMask = getMask()->rank(start + length - 1);
                
            // Cound how many we touch because of that. TODO: Unit test this!
            return endMask - startMask + 1;
        } else {
            // Count all the BWT positions.
            return length;
        }
    }
}

size_t FMDIndexView::getApproximateNumberOfNewRanges(size_t oldStart,
    size_t oldLength, size_t newStart, size_t newLength) const {

    // We can just make the new ranges, keeping track of only the forward
    // intervals.
    
    // What's the new range on the left?
    size_t newLeftStart = newStart;
    size_t newLeftLength = oldStart - newStart;

    // And on the right?
    size_t newRightStart = oldStart + oldLength;
    size_t newRightLength = newStart + newLength - oldStart - oldLength;
    
    // Sum up the number of ranges on both sides.
    return getApproximateNumberOfRanges(newLeftStart, newLeftLength) +
        getApproximateNumberOfRanges(newRightStart, newRightLength);

}







