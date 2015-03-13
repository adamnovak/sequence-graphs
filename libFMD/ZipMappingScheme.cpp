#include "ZipMappingScheme.hpp"
#include "Log.hpp"
#include "util.hpp"

#include <vector>
#include <algorithm>
#include <cstdlib>

#include <boost/range/adaptor/reversed.hpp>

std::vector<FMDPosition> ZipMappingScheme::findRightContexts(
    const std::string& query) const {
 
    // This will hold the right context of every position.
    std::vector<FMDPosition> toReturn(query.size());
 
    // Start with everything selected.
    FMDPosition results = view.getIndex().getCoveringPosition();
    
    // How many characters are currently searched?
    size_t patternLength = 0;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--) {
        // For each position in the query from right to left, we're going to
        // inchworm along and get the FMDPosition that is extended out right as
        // far as possible while still having results.
        
        // We're going to extend left with this new base.
        FMDPosition extended = results;
        view.getIndex().extendLeftOnly(extended, query[i]);
        
        while(view.isEmpty(extended)) {
            // If you couldn't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPosition retracted = results;
            // Make sure to drop characters from the total pattern length.
            view.getIndex().retractRightOnly(retracted, --patternLength);
            
            // Try extending again
            extended = retracted;
            view.getIndex().extendLeftOnly(extended, query[i]);
            
            // Say that last step we came from retracted.
            results = retracted;
        }
    
        // We successfully extended (if only with the base itself).
        results = extended;
        // Increment the pattern length since we did actually extedn by 1. 
        patternLength++;
        
        // Save the search results.
        toReturn[i] = results;
    }

    // Now we have inchwormed all the way from right to left, retracting only
    // when necessary. Return all the results.    
    return toReturn;
}

void ZipMappingScheme::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // Get the right contexts
    Log::info() << "Looking for right contexts" << std::endl;
    std::vector<FMDPosition> rightContexts = findRightContexts(query);
    
    // And the left contexts (which is the reverse of the contexts for the
    // reverse complement.
    Log::info() << "Looking for left contexts" << std::endl;
    std::vector<FMDPosition> leftContexts = findRightContexts(
        reverseComplement(query));
    std::reverse(leftContexts.begin(), leftContexts.end());
    
    // For each pair, figure out if the forward and reverse FMDPositions select
    // one single consistent TextPosition.
    for(size_t i = 0; i < query.size(); i++) {
    
        Log::debug() << "Base " << i << " = " << query[i] << " selects " <<
            rightContexts[i] << " and " << leftContexts[i] << std::endl;
    
        // Find the forward-orientation positions selected by the right context.
        std::set<TextPosition> forwardPositions = view.getTextPositions(
            rightContexts[i]);
            
        // And the reverse-orientation positions selected by the left contexts.
        std::set<TextPosition> reversePositions = view.getTextPositions(
            leftContexts[i]);
            
        // We're going to keep a set of TextPositions that appeared forward in
        // the forward set and backwards in the reverse set here.
        std::set<TextPosition> intersection;
            
        for(auto reversePosition : reversePositions) {
            // Flip our temporary TextPosition in place. TODO: Can we make this
            // use pattern shorter? Maybe pass the index or have the index do it
            // or something?
            reversePosition.flip(view.getIndex().getContigLength(
                reversePosition.getContigNumber()));
            
            if(forwardPositions.count(reversePosition)) {
                intersection.insert(reversePosition);
                Log::debug() << "\tTextPosition " << reversePosition <<
                    " is shared." << std::endl;
            }
            
            // TODO: Fail mapping fast if we find a nonempty intersection.
        }
        
        // TODO: If we can't find anything, try retracting a few bases on one or
        // both sides until we get a shared result.
        
        if(intersection.size() == 1) {
            // We map! Call the callback with the only text position.
            Log::debug() << "Index " << i << " maps." << std::endl;
            callback(i, *(intersection.begin()));
        } else if(intersection.size() > 1) {
            // Too many results
            Log::debug() << "Index " << i << " is ambiguous." << std::endl;
        } else {
            // Too few results
            Log::debug() << "Index " << i << " not in reference." << std::endl;
        }
        
    }
    
    // Now we have gone through every position in the query and found the text
    // positions it matches for both forward and reverse complement.
    
}

