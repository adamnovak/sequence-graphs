#include "ResultRange.hpp"

ResultRange::ResultRange(): position(EMPTY_FMD_POSITION), searchStringStart(0),
    searchStringEnd(0), mismatches() {

    // Nothing to do
}

ResultRange::ResultRange(const FMDIndex& field, size_t queryLength): 
    position(field.getCoveringPosition()), searchStringStart(queryLength),
    searchStringEnd(queryLength), mismatches() {
    
    // Nothing to do
}

ResultRange ResultRange::extendLeftMatch(const FMDIndex& index,
    const std::string& query) const {
    
    // Copy ourselves
    ResultRange toReturn(*this);
    
    if(searchStringStart != 0) {
        // It hasn't hit the start of the query yet.
        
        // Include another base
        toReturn.searchStringStart--;
        
        // And extend with it (possibly making the range empty).
        index.extendLeftOnly(toReturn.position,
            query[toReturn.searchStringStart]);
    } else {
        // We can't extend with something that isn't there.
        toReturn.position = EMPTY_FMD_POSITION;
    }
    
    return toReturn;
}

std::array<ResultRange, 3> ResultRange::extendLeftMismatch(
    const FMDIndex& index, const std::string& query) const {

    // Make the array, copying ourselves 3 times.
    std::array<ResultRange, 3> toReturn = {*this, *this, *this};
    
    
    // What's the next result to fill in?
    size_t nextResult = 0;
    
    if(searchStringStart != 0) {
        // We haven't hit the start of the query yet.
        
        // Grab the character that is actually there.
        char matchChar = query[searchStringStart - 1];
        
        for(char base : BASES) {
            // Try all the bases
            
            if(base == matchChar) {
                // Skip the match base
                continue;
            }
            
            // So for all the mismatch bases, extend the next available result
            // range with this mismatch character.
            toReturn[nextResult].searchStringStart--;
            index.extendLeftOnly(toReturn[nextResult].position, base);
            // Make sure to log the mismatch.
            toReturn[nextResult].mismatches.push_back(
                toReturn[nextResult].searchStringStart);
            
            // And use up that result slot.
            nextResult++;
        }
        
    } else {
        // We can't extend with something that isn't there.
        
        for(ResultRange& range : toReturn) {
            // Mark each result empty.
            range.position = EMPTY_FMD_POSITION;
        }
    }
    
    // Return the array of results.
    return toReturn;
}
ResultRange ResultRange::retractRight(const FMDIndex& index) const {
    // Copy ourselves.
    ResultRange toReturn(*this);
    
    if(toReturn.mismatches.front() == toReturn.searchStringEnd) {
        // We're retracting off a mismatch. Remove it.
        toReturn.mismatches.pop_front();
    }
    
    // Throw a character out of the search string.
    toReturn.searchStringEnd--;
    
    // Retract by one character.
    index.retractRightOnly(toReturn.position, toReturn.getSearchStringLength());
    
    return toReturn;
}

bool ResultRange::isEmpty(GenericBitVector* mask) const {
    return position.isEmpty(mask);
}

size_t ResultRange::getLength(GenericBitVector* mask) const {
    return position.getLength(mask);
}

size_t ResultRange::getSearchStringLength() const {
    return searchStringEnd - searchStringStart;
}

size_t ResultRange::mismatchesUsed() const {
    return mismatches.size();
}

bool ResultRange::operator==(const ResultRange& other) const {
    
    
    return(position == other.position && 
        searchStringStart == other.searchStringStart && 
        searchStringEnd == other.searchStringEnd && 
        mismatches.size() == other.mismatches.size() && 
        std::equal(mismatches.begin(), mismatches.end(), 
        other.mismatches.begin()));
        
}











