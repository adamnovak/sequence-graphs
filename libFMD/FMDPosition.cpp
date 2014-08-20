#include <stdexcept>

#include "FMDPosition.hpp"

FMDPosition::FMDPosition(int64_t forward_start, int64_t reverse_start,
  int64_t end_offset): forward_start(forward_start), 
  reverse_start(reverse_start), end_offset(end_offset) {
}

FMDPosition::FMDPosition(): forward_start(0), reverse_start(0), end_offset(-1) {
}

FMDPosition FMDPosition::flip() const {
    // Swap the two intervals of the bi-interval
    return FMDPosition(reverse_start, forward_start, end_offset);
}

bool FMDPosition::operator==(const FMDPosition& other) const {
    // Compare all the fields.
    return
        forward_start == other.forward_start &&
        reverse_start == other.reverse_start &&
        end_offset == other.end_offset;
}

int64_t FMDPosition::range(const GenericBitVector& ranges, 
    const GenericBitVector* mask) const {

    // What's the first index in the BWT that the range would need to cover?
    size_t interval_start;
    // What's the last index in the BWT that the range would need to cover?
    size_t interval_end;

    if(mask != NULL) {
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
        
        // Algorithm:
        // Find the first 1 in the mask in our interval.
        // Find the last 1 in the mask in our interval.
        // See what range spans them, if any.
        
        // Find the first 1 at position >= our start. This holds (index, rank-1)
        // TODO: Convert to std::tie
        CSA::pair_type firstOne = mask->valueAfter(forward_start);
        // TODO: What if there is no such 1?
        
        // Find the last 1 at position <= our end. This holds (index, rank-1),
        // or (size, items) if no such 1 exists.
        CSA::pair_type lastOne = mask->valueBefore(forward_start + end_offset);
        if(lastOne.first > forward_start + end_offset) {
            // We didn't find such a 1. There must be no 1 in our interval and
            // we are thus empty.
            return -1;
        }
        
        // Fill in the specs for the interval we want.
        interval_start = firstOne.first;
        interval_end = lastOne.first;
    } else {
        // Fast path; no need to mess about with the mask.
        interval_start = forward_start;
        interval_end = forward_start + end_offset;   
    }

    if(interval_end < interval_start || end_offset < 0) {
        // Our range is actually empty.
        // TODO: Do we need this check?
        return -1;
    }

    // If we get here we actually have a range to try to look up

    // Look up the range that the forward starting position is in
    int64_t start_range = ranges.rank(interval_start);

    // And the range the forward end is in
    int64_t end_range = ranges.rank(interval_end);
    
    if(start_range == end_range) {
        // Both ends of the interval are in the same range.
        return start_range;
    } else {
        // There is no range spanning our forward-strand interval.
        return -1;
    }
}

int64_t FMDPosition::ranges(const GenericBitVector& ranges, 
    const GenericBitVector* mask) const {

    if(end_offset < 0) {
        // Special case the empty FMDPosition so we don't have to deal with that
        // case later.
        return 0;
    } else if(mask == NULL) {
        // Fast path. Since there's no mask, every range counts.
    
        // Look up the range that the starting position is in
        int64_t start_range = ranges.rank(forward_start);

        // And the range the end is in
        int64_t end_range = ranges.rank(forward_start + end_offset);
        
        if(end_range < start_range) {
            throw std::runtime_error("End of ranges interval is before start");
        }
        
        // Return the number of ranges we intersect (1s hit plus 1)
        return end_range - start_range + 1;
    
    } else {
        // Slow path accounting for mask. Not every range counts, only those
        // that overlap positions with 1s in the mask.

        // Keep track of the left edge of the interval we still need to count
        // the ranges in.
        size_t left = forward_start;

        // How many ranges have we found so far?
        int64_t rangeCount = 0;
        
        while(left <= forward_start + end_offset) {
            // Until the next range starts outside this FMDPosition...
        
            // Find the index and rank of the first 1 in the mask that's here or
            // beyond our current position.
            CSA::pair_type nextPosition = mask->valueAfter(left);
            
            if(nextPosition.first <= forward_start + end_offset) {
                // This masked-in position is within us. We should count its
                // range.
                rangeCount++;
                
                // Find the next 1 in the ranges vector after this masked-in
                // position (exclusive), indicating the start of the next range.
                CSA::pair_type nextRange = ranges.valueAfter(
                    nextPosition.first + 1);
                    
                // Look in or after that range for the next masked-in position.
                left = nextRange.first;
            } else {
                // The next masked-in position is out of us. We're done counting
                // ranges.
                break;
            }
            
        }
        
        // Return the count of occupied ranges.
        return rangeCount;
    }
}

std::ostream& operator<< (std::ostream& o, FMDPosition const& position) {
    // Report both the ranges that we represent.
    return o << position.forward_start << "-" << 
        (position.forward_start + position.end_offset) << "|" << 
        position.reverse_start << "-" << (position.reverse_start + 
        position.end_offset);
}
