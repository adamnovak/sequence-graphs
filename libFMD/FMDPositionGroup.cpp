#include "FMDPositionGroup.hpp"

#include "util.hpp"

FMDPositionGroup::FMDPositionGroup(const FMDIndexView& view): positions(),
    view(view) {
    // Nothing to do!
}

FMDPositionGroup::FMDPositionGroup(const FMDIndexView& view,
    const std::vector<FMDPosition>& startPositions): FMDPositionGroup(view) {
    
    for(const auto& position : startPositions) {
        // We need to tag all the positions as having 0 mismatches
        positions.emplace_back(position, 0);
    }
}

void FMDPositionGroup::extendGreedy(char correctCharacter,
    size_t maxMismatches) {

    // This will hold all the extensions we find that aren't empty
    std::vector<std::pair<FMDPosition, size_t>> nonemptyExtensions;
    
    for(const auto& positionAndMismatches : positions) {
        // For each existing FMDPosition
        
        // Pull it out
        FMDPosition toExtend = positionAndMismatches.first;
        
        // Extend it
        view.getIndex().extendLeftOnly(toExtend, correctCharacter);
        
        if(!view.isEmpty(toExtend)) {
            // If we got anything, keep the range (and don't increment the
            // mismatches)
            nonemptyExtensions.emplace_back(toExtend,
                positionAndMismatches.second);
        }
    }
    
    if(nonemptyExtensions.size() == 0) {
        // If they are all empty, extend with all the mismatch characters and
        // charge a mismatch each.
    
        for(char base : BASES) {
            if(base == correctCharacter) {
                // Don't even try the right base, we just did.
                continue;
            }
            
            for(const auto& positionAndMismatches : positions) {
                // Try each position we have as a starting point.
                
                if(positionAndMismatches.second >= maxMismatches) {
                    // We can't afford another mismatch on this one.
                    continue;
                }
                
                // Pull it out
                FMDPosition toExtend = positionAndMismatches.first;
                
                // Extend it
                view.getIndex().extendLeftOnly(toExtend, base);
                
                if(!view.isEmpty(toExtend)) {
                    // If we got anything, keep the range (and charge a
                    // mismatch)
                    nonemptyExtensions.emplace_back(toExtend,
                        positionAndMismatches.second + 1);
                }
            }
            
        }
        
        Log::debug() << "Extended with all but " << correctCharacter <<
            std::endl;
        
    } else {
        Log::debug() << "Extended with " << correctCharacter << std::endl;
    }
    
    // Replace our FMDPositions with the new extended ones.
    positions = std::move(nonemptyExtensions);
    
    Log::debug() << "Have " << positions.size() << " ranges" << std::endl;

}


bool FMDPositionGroup::isEmpty() const {
    for(const auto& positionAndMismatches : positions) {
        // For each interval we contain
        if(!view.isEmpty(positionAndMismatches.first)) {
            // If any is nonempty, we are nonempty
            return false;
        }
    }
    
    // If none exists or they are all empty, we are empty.
    return true;
}

bool FMDPositionGroup::isUnique() const {
    for(const auto& positionAndMismatches : positions) {
        // For each interval we contain
        if(view.isAmbiguous(positionAndMismatches.first)) {
            // If any is ambiguous, we are non-unique
            return false;
        }
    }
    
    // If they are all unique, we need to make sure they are all unique to the
    // same place.
    
    // This is that place
    TextPosition uniquePlace;
    
    // This is a flag that we only set once we actually find a nonempty
    // FMDPosition.
    bool found = false;
    
    for(const auto& positionAndMismatches : positions) {
        // For each position
        auto position = positionAndMismatches.first;
        
        if(view.isEmpty(position)) {
            continue;
        }
        
        if(found) {
            // We already have a candidate unique place, so make sure this one
            // matches.
            if(uniquePlace != view.getTextPosition(position)) {
                return false;
            }
        } else {
            // This is the first one we found
            
            // Now we know it's nonempty and unique
            found = true;
            // Figure out where to
            uniquePlace = view.getTextPosition(position);
        }
    }
    
    // If we found something we're unique. Otherwise we're empty.
    return found;
}

TextPosition FMDPositionGroup::getTextPosition() const {
    // TODO: Unify with above function, since they differ only in what they
    // return.

    for(const auto& positionAndMismatches : positions) {
        // For each interval we contain
        if(view.isAmbiguous(positionAndMismatches.first)) {
            // If any is ambiguous, we are non-unique
            throw std::runtime_error(
                "No TextPosition for ambiguous FMDPositionGroup");
        }
    }
    
    // If they are all unique, we need to make sure they are all unique to the
    // same place.
    
    // This is that place
    TextPosition uniquePlace;
    
    // This is a flag that we only set once we actually find a nonempty
    // FMDPosition.
    bool found = false;
    
    for(const auto& positionAndMismatches : positions) {
        // For each position
        auto position = positionAndMismatches.first;
        
        if(view.isEmpty(position)) {
            continue;
        }
        
        if(found) {
            // We already have a candidate unique place, so make sure this one
            // matches.
            if(uniquePlace != view.getTextPosition(position)) {
                 throw std::runtime_error(
                    "No TextPosition for ambiguous FMDPositionGroup");
            }
        } else {
            // This is the first one we found
            
            // Now we know it's nonempty and unique
            found = true;
            // Figure out where to
            uniquePlace = view.getTextPosition(position);
        }
    }
    
    // If we found something we're unique. Otherwise we're empty.
    return uniquePlace;
}















