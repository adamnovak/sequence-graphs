#include "CreditStrategy.hpp"

#include "util.hpp"
#include "Log.hpp"

CreditStrategy::CreditStrategy(const FMDIndexView& view, bool enabled): 
    view(view), enabled(enabled) {
    // Nothing to do!
}

void CreditStrategy::applyCredit(const std::string& query, 
    std::vector<Mapping>& toUpdate) const {

    if(!enabled) {
        // We belong to a mapping scheme and someone has disabled us (probably
        // in response to command-line options). Don't do anything.
        return;
    }

    // Scan the mappings from left to right to find pairs of mapped positions
    // bordering runs of unmapped positions.
    
    // What is the last mapped position you saw?
    size_t leftAnchor = (size_t) -1;
    
    // Was the last position you looked at mapped?
    bool lastWasMapped = true;
    
    for(size_t i = 0; i < toUpdate.size(); i++) {
        // Look at every base and observe where we change from mapped to
        // unmapped or visa versa.
        
        Log::debug() << "Checking position " << i << " (mapped: " <<
            toUpdate[i].isMapped() << ", last: " << lastWasMapped << ")" <<
            std::endl;
    
        if(!toUpdate[i].isMapped() && lastWasMapped && i != 0) {
            // We have just entered a run of unmapped bases. We have to set the
            // left anchor to the position before here.
            leftAnchor = i - 1;
            
            Log::debug() << "No longer mapped" << std::endl;
            
        } else if(toUpdate[i].isMapped() && !lastWasMapped &&
            leftAnchor != (size_t) -1) {
            
            // We have just finished a run of unmapped bases. We need to run the
            // searches.
            
            Log::debug() << "Newly mapped" << std::endl;
            
            // We should use this base as the right anchor of the credit.
            size_t rightAnchor = i;
            
            // Actually go and apply the credit between those mapped bases.
            applyCreditBetween(query, toUpdate, leftAnchor, rightAnchor);
        }
        
        // Remember the state of the last base we saw.        
        lastWasMapped = toUpdate[i].isMapped();
    }
}

void CreditStrategy::applyCreditBetween(const std::string& query, 
    std::vector<Mapping>& toUpdate, size_t leftAnchor,
    size_t rightAnchor) const {
    
    Log::info() << "Applying credit between " << leftAnchor << " and " << 
        rightAnchor << std::endl;
    
    // How many bases do we want to give credit to?
    size_t unmappedRegionSize = rightAnchor - leftAnchor - 1; 
    
    // Make a vector of unmapped mappings for credit from the left.
    std::vector<Mapping> leftMappings(unmappedRegionSize, Mapping());
    // And from the right
    std::vector<Mapping> rightMappings(unmappedRegionSize, Mapping());
    
    // Run a search in from each end, to the depth of the other end or until you
    // find too many mismatches.
    
    // Do the right-side search
    breadthFirstSearch(query, rightAnchor, toUpdate[rightAnchor].getLocation(), 
        unmappedRegionSize, [&](size_t index, TextPosition mappedTo) {
        
        // Save every mapping as coming from the right-side search. Make sure to
        // convert from global to bubble-local coordinates.
        rightMappings[index - (leftAnchor + 1)] = Mapping(mappedTo);
    
    });
    
    // Flip the query around. TODO: do this efficiently by using a special DNA
    // string type or soemthing.
    std::string rcQuery = reverseComplement(query);
    
    // Get the left starting position, flipped around.
    TextPosition leftFlipped = toUpdate[leftAnchor].getLocation();
    // TODO: We need the contig length to flip it, and it is hard to get.
    leftFlipped.flip(view.getIndex().getContigLength(
        leftFlipped.getContigNumber()));
    
    // Do the left-side search
    breadthFirstSearch(rcQuery, rcQuery.size() - leftAnchor - 1, leftFlipped, 
        unmappedRegionSize, [&](size_t index, TextPosition mappedTo) {
        
        // Flip around again
        index = rcQuery.size() - index - 1;
        mappedTo.flip(view.getIndex().getContigLength(
            mappedTo.getContigNumber()));
        
        // Save every mapping as coming from the left-side search. Make sure to
        // convert from global to bubble-local coordinates.
        leftMappings[index - (leftAnchor + 1)] = Mapping(mappedTo);
    
    });
    
    // Make sure they agree and apply them
    for(size_t i = 0; i < leftMappings.size(); i++) {
        // For each position in the region we were trying to give credit to
        
        if(leftMappings[i].isMapped() && !rightMappings[i].isMapped()) {
            // Apply the mapping from the left, since it's the only one
            toUpdate[leftAnchor + i + 1] = leftMappings[i];
        } else if(rightMappings[i].isMapped() && !leftMappings[i].isMapped()) {
            // Apply the mapping from the right, since it's the only one.
            toUpdate[leftAnchor + i + 1] = rightMappings[i];
        } else if (leftMappings[i] == rightMappings[i] &&
            leftMappings[i].isMapped()) {
            // Both mapped but they agree, so apply both.
            toUpdate[leftAnchor + i + 1] = leftMappings[i];
        }
    }
    
    // Now we've updated the vector of mappings.
}

void CreditStrategy::breadthFirstSearch(const std::string& query,
    size_t queryStart, const TextPosition& referenceStart, size_t maxDepth, 
    std::function<void(size_t, TextPosition)> callback) const {

    // See where the reference starts and turn that into range numbers
    auto ranges = view.textPositionToRanges(referenceStart);
    
    // We're going to turn them all into FMDPositions
    std::vector<FMDPosition> startPositions;
    
    for(auto range : ranges) {
        // Each range gets expanded into its covering FMDPosition.
        std::pair<size_t, size_t> startAndLength = view.getRangeByNumber(range);
        startPositions.emplace_back(startAndLength.first, 0,
            ((int64_t) startAndLength.second) - 1);
        
        Log::debug() << "Got range " << startAndLength.first << "-" <<
            startAndLength.first + startAndLength.second << " as #" << range <<
            std::endl;
    }

    // Make the FMDPositionGroup holding them all
    FMDPositionGroup search(startPositions);
    
    Log::trace() << "Applying credit from " << queryStart - 1 << " to " << 
        queryStart - maxDepth << std::endl;
    
    for(size_t queryIndex = queryStart - 1; 
        (queryIndex >= queryStart - maxDepth) && 
        (queryIndex != (size_t) -1) && !search.isEmpty(view); 
        queryIndex--) {
        
        // While it is not empty, and we haven't hit our depth limit (or run off
        // the left edge of the string)
    
        // Try extending with a character. Remember if we found that character or not.
        bool correctCharacter = search.extendGreedy(view, query[queryIndex],
            maxMismatches);
        
        if(correctCharacter && search.isUnique(view)) {
            // If we have a unique result, and we actually found the base we
            // were looking for, issue a callback.
            callback(queryIndex, search.getTextPosition(view));
        }
        
    }
}
