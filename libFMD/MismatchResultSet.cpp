#include "MismatchResultSet.hpp"
#include "FMDIndex.hpp"

MismatchResultSet::MismatchResultSet(const FMDIndex& index): 
    index(index), results({MismatchResult(index)}) {
    // We started with a single result covering the whole index.
    
    // Nothing to do!
}

void MismatchResultSet::extendLeft(char base) {
    // Make a new set to hold all the new results after extension.
    std::set<MismatchResult> newResults;
    
    for(auto oldResult : results) {
        // For each current result
        for(auto newResult : oldResult.extendLeft(index, base)) {
            // For each new result you would get by extending it, put it in the
            // set of new results (if we don't have it already, which we
            // shouldn't, because these are all subranges of originally unique
            // ranges).
            newResults.insert(newResult);
        }
    }
    
    // Keep the new results. Swap to avoid copying the whole thing; the old set
    // can go out of scope here.
    // TODO: Is there a less clever way to do this now? Some kind of move?
    std::swap(results, newResults);
}

void MismatchResultSet::retractRight() {
    // Make a new set to hold all the new results after retracting.
    std::set<MismatchResult> newResults;
    
    for(auto oldResult : results) {
        // For each current result, retract it and throw that in the new set.
        // This will never make an empty result, but it will make duplicates,
        // which is why we're using the set in the first place.
        newResults.insert(oldResult.retractRight(index));
    }
    
    // Replace the old result set.
    std::swap(results, newResults);
}

void MismatchResultSet::filter(size_t z_max, const GenericBitVector* mask) {
    // Make a new set to hold all the new results after retracting.
    std::set<MismatchResult> newResults;
    
    // TODO: Don't have this method, add these filters when things get added.
    // This is going to be O(n log n) or something due to set being some fancy
    // tree thing and not O(1) access.
    // TODO: If I move over to an unordered_set would this be easier?
    
    for(auto result : results) {
        // For each current result
        
        if(result.getMismatches() <= z_max && !result.isEmpty(mask)) {
            // Keep the ones that are sufficiently good.
            newResults.insert(result);
        }
    }
    
    // Replace the old result set.
    std::swap(results, newResults);
}

int64_t MismatchResultSet::range(const GenericBitVector& ranges, 
    const GenericBitVector* mask) const {
    
    // What range has been found?
    int64_t range = -1;
    
    for(auto result : results) {
        // For each result, get its range.
        int64_t resultRange = result.range(ranges);
        
        if(range == -1) {
            // Adopt this as our range, since we have none.
            range = resultRange;
        }
        
        if(resultRange == -1 || range != resultRange) {
            // No range, or a contradictory one, for this result.
            return -1;
        }
    }
    
    // Everybody agrees on this. Will be -1 if no results.
    return range;

}

bool MismatchResultSet::isEmpty(const GenericBitVector* mask) const {
    for(auto result : results) {
        // Go through each result
        
        if(!result.isEmpty(mask)) {
            // And if it has anything under the mask, say we aren't empty
            return false;
        }
    }
    
    // Otherwise we are empty.
    return true;
}

int64_t MismatchResultSet::getLength(const GenericBitVector* mask) const {
    // How many total result positions are there?
    int64_t total = 0;
    
    for(auto result : results) {
        // For each result
    
        // The positions in different results will all be distinct, so we can
        // count all of them.
        total += result.getLength(mask);
    }
    
    // Return the total number of BWT positions in results in this result set,
    // under the mask.
    return total;
}

bool MismatchResultSet::isMapped(const GenericBitVector& ranges, 
    const GenericBitVector* mask) const {
    
    // We are mapped if we overlap exactly one range.
    return range(ranges, mask) != -1;    
}
