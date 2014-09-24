#include "CreditFilter.hpp"
#include <algorithm>

CreditFilter::CreditFilter(FMDIndex& index): index(index) {
    // Nothing to do!
}

std::vector<Mapping> CreditFilter::apply(
    std::vector<Mapping> leftMappings, std::vector<Mapping> rightMappings) {

    // First disambiguate everything
    std::vector<Mapping> disambiguated = disambiguate.apply(leftMappings,
        rightMappings);

    // Find the left and right sentinels (must be mapeped on the left/right
    // respectively, and not have the other side disagree). Maybe this can be
    // solved by disambiguation.
    
    int64_t leftSentinel = -1;
    int64_t rightSentinel = -1;
    
    for(size_t i = 0; i < disambiguated.size(); i++) {
        // Look in from the left for the left sentinel
        if(leftMappings[i].isMapped() && 
            disambiguated[i].isMapped()) {
            // This is the first mapped thing we have found mapped on the left.
            leftSentinel = i;
            break;
        }
    }
    
    for(size_t i = disambiguated.size() - 1; i != (size_t) -1; i--) {
        // Look in from the right for the right sentinel
        if(rightMappings[i].isMapped() && 
            disambiguated[i].isMapped()) {
            // This is the first mapped thing we have found mapped on the right.
            rightSentinel = i;
            break;
        }
    }
    
    if(leftSentinel == -1 || rigthSentinel == -1) {
        // We didn't find one of the sentinels, so we can't do credit between
        // them. Just disambiguate without applying credit.
        return std::move(disambiguated);
    }
    
    // Find the max left and right contexts we need to worry about checking
    // consistency from.
    size_t maxLeftContext = 0;
    size_t maxRightContext = 0;
    
    for(size_t i = 0; i < leftMappings.size(); i++) {
        // Take the max left and right context when we find them.
        maxLeftContext = std::max(maxLeftContext, leftMappings[i].geContext());
        maxRightContext = std::max(maxRightContext, 
            rightMappings[i].geContext());
    }
    
    // This is going to hold our output
    std::vector<Mapping> toReturn;
    
    for(size_t i = 0; i <= leftSentinel; i++) {
        // For each base before or at the left sentinel, disambiguate normally.
        toReturn.push_back(disambiguated[i]);
    }

    for(size_t i = leftSentinel + 1; i < rightSentinel; i++) {    
        // For each base between the sentinels
    
        if(disambiguated[i].isMapped()) {
            // Look at how it has been mapped
            // If it is mapped on one or more sides, disambiguate normally.
            toReturn.push_back(disambiguated[i]);   
        } else {
        
            // Set this to true if you find a base that implies a position for
            // this one by its right context.
            bool rightFound = false;
            // This is where that base would place this one
            TextPosition rightCreditPosition;
            // Set this if right contexts all place this base in one spot
            bool rightConsistent = true;
        
            for(size_t j = i - 1; j != (size_t) -1 && j >= i - maxRightContext;
                j--) {
                // Look left from here until you find a base that maps and
                // places us.
                
                if(!rightMappings[j].isMapped()) {
                    // This base never mapped, so it can't give credit.
                    continue;
                }
                
                if(rightMappings[j].getContext() - 1 < i - j) {
                    // This base's context didn't reach all the way out, after
                    // accounting for the fact that it includes the base itself.
                    continue;
                }
                
                // OK, we imply some mapping. What is it?
                // Grab the mapped base and go forwards on the forward strand
                // (since right contexts reach forward).
                TextPosition implied = rightMappings[j].getLocation();
                implied.addOffset(i - j);
                
                if(!rightFound) {
                    // This is the first one. Make sure all the others match.
                    rightFound = true;
                    rightCreditPosition = implied;
                } else if (rightCreditPosition != implied) {
                    // There are two or more implied locations from right
                    // contexts.
                    rightConsistent = false;
                }
                
            }
            
            // And the same thing for left contexts
            
            // Set this to true if you find a base that implies a position for
            // this one by its left context.
            bool leftFound = false;
            // This is where that base would place this one
            TextPosition leftCreditPosition;
            // Set this if left contexts all place this base in one spot
            bool leftConsistent = true;
        
            for(size_t j = i + 1; j < i + maxLeftContext; j++) {
                // Look right from here until you find a base that maps and
                // places us.
                
                if(!leftMappings[j].isMapped()) {
                    // This base never mapped, so it can't give credit.
                    continue;
                }
                
                if(leftMappings[j].getContext() - 1 < j - i) {
                    // This base's context didn't reach all the way out, after
                    // accounting for the fact that it includes the base itself.
                    continue;
                }
                
                // OK, we imply some mapping. What is it?
                // Grab the mapped base and go backward on the forward strand
                // (since left contexts reach backward).
                TextPosition implied = leftMappings[j].getLocation();
                implied.addOffset((int64_t) i - (int64_t) j);
                
                if(!leftFound) {
                    // This is the first one. Make sure all the others match.
                    leftFound = true;
                    leftCreditPosition = implied;
                } else if (leftCreditPosition != implied) {
                    // There are two or more implied locations from left
                    // contexts.
                    leftConsistent = false;
                }
                
            }
            
            // TODO: somehow call into the disambiguation filter here?
        
            if(leftFound && leftConsistent) {
                if(rightFound && rightConsistent) {
                    // We have credit from both the left and the right. Do they
                    // agree?
                    
                    // Flip the text position we got from the right using the
                    // index.
                    TextPosition rightFlipped = rightCreditPosition.flip(
                        index.getContigLength(index.getContigNumber(
                        rightCreditPosition)));
                        
                    if(leftCreditPosition == rightFlipped) {
                        // They agree! Record a successful mapping.
                        toReturn.push_back(Mapping(leftCreditPosition));
                    } else {
                        // They disagree. Fail the mapping
                        toReturn.push_back(Mapping());
                    }
                } else {
                    // We have credit only from the left.
                    toReturn.push_back(Mapping(leftCreditPosition));
                }
            } else if(rightFound && rightConsistent) {
                // We have credit only from the right
                
                // We need to flip it into left semantics.
                TextPosition rightFlipped = rightCreditPosition.flip(
                    index.getContigLength(index.getContig(
                    rightCreditPosition)));
            
                toReturn.push_back(Mapping(rightFlipped));
            } else {
                // We have credit from nowhere
                toReturn.push_back(Mapping());
            }
        }
        
        
        
    }
    
    for(size_t i = rightSentinel; i < disambiguated.size(); i++) {
        // Then for each base at or after the right sentinel, disambiguate
        // normally.
        toReturn.push_back(disambiguated[i]);
    }
    
    // Return the credit-mapped stuff
    return std::move(toReturn);
    
}
