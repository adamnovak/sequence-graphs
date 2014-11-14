#include "MismatchResult.hpp"
#include "util.hpp"

std::vector<MismatchResult> MismatchResult::extendLeft(const FMDIndex& index, 
    char match) const {
    
    // Make a vector to hold our result.
    std::vector<MismatchResult> toReturn;
    
    for(auto base : BASES) {
        // For each base we can try extending with...
        
        // Copy ourselves
        MismatchResult extended(*this);
        
        // Extend the copy on that base.
        index.extendLeftOnly(extended.result, base);
        
        if(extended.result.isEmpty()) {
            // Skip any extensions that are unconditionally empty.
            continue;
        }
        
        // Note that we searched another character.
        extended.endIndex++;
        
        if(base != match) {
            // Note that there is a mismatch here.
            extended.mismatches.push_back(extended.endIndex);
        }
        
        // Keep this extension around.
        toReturn.push_back(extended);
    }
    
    // Give back the four extensions.
    return toReturn;

}
    
MismatchResult MismatchResult::retractRight(const FMDIndex& index) const {
    // Copy ourselves
    MismatchResult toReturn(*this);
    
    if(toReturn.rightIsMismatch()) {
        // We're about to retract a mismatched character. Remove the record of
        // its mismatch.
        toReturn.mismatches.pop_front();
    }
    
    // Say we have searched one fewwe character.
    toReturn.startIndex++;
    
    // Retract the result to match our recorded search length.
    index.retractRightOnly(toReturn.result, toReturn.getCharacters());
}
