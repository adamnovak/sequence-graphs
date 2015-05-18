#include "FMDPositionGroup.hpp"

#include "util.hpp"

FMDPositionGroup::FMDPositionGroup(const FMDIndexView& view): positions(),
    view(view) {
    // Nothing to do!
}

FMDPositionGroup::FMDPositionGroup(const FMDIndexView& view,
    const std::vector<FMDPosition>& startPositions): FMDPositionGroup(view) {
    
    for(const auto& position : startPositions) {
        // We need to tag all the positions with empty annotations
        positions.emplace_back(position);
    }
}

bool FMDPositionGroup::extendGreedy(char correctCharacter,
    size_t maxMismatches) {

    // This will hold all the extensions we find that aren't empty
    decltype(positions) nonemptyExtensions;
    
    // Set this to true if we found an exact match, false if we had to use a
    // mismatch here (and thus, even if we are unique, we shouldn't produce a
    // mapping).
    bool exactMatch;
    
    for(const auto& annotated : positions) {
        // For each existing FMDPosition
        
        // Pull it out
        FMDPosition toExtend = annotated.position;
        
        // Extend it
        view.getIndex().extendLeftOnly(toExtend, correctCharacter);
        
        if(!view.isEmpty(toExtend)) {
            // If we got anything, keep the range (and don't increment the
            // mismatches)
            nonemptyExtensions.emplace_back(toExtend, annotated, false);
        }
    }
    
    if(nonemptyExtensions.size() == 0) {
        // If they are all empty, extend with all the mismatch characters and
        // charge a mismatch each.
    
        // We have to use mismatches, if we can find anything at all.
        exactMatch = false;
    
        for(char base : BASES) {
            if(base == correctCharacter) {
                // Don't even try the right base, we just did.
                continue;
            }
            
            for(const auto& annotated : positions) {
                // Try each position we have as a starting point.
                
                if(annotated.mismatches >= maxMismatches) {
                    // We can't afford another mismatch on this one.
                    continue;
                }
                
                // Pull it out
                FMDPosition toExtend = annotated.position;
                
                // Extend it
                view.getIndex().extendLeftOnly(toExtend, base);
                
                if(!view.isEmpty(toExtend)) {
                    // If we got anything, keep the range (and charge a
                    // mismatch)
                    nonemptyExtensions.emplace_back(toExtend, annotated, true);
                }
            }
            
        }
        
        Log::debug() << "Extended with all but " << correctCharacter <<
            std::endl;
        
    } else {
        // We managed to extend with an exact match.
        exactMatch = true;
        Log::debug() << "Extended with " << correctCharacter << std::endl;
    }
    
    // Replace our FMDPositions with the new extended ones.
    positions = std::move(nonemptyExtensions);
    
    Log::debug() << "Have " << positions.size() << " ranges" << std::endl;
    
    // Let the caller know if we found the base they wanted or not.
    return exactMatch;

}


bool FMDPositionGroup::isEmpty() const {
    for(const auto& annotated : positions) {
        // For each interval we contain
        if(!view.isEmpty(annotated.position)) {
            // If any is nonempty, we are nonempty
            return false;
        }
    }
    
    // If none exists or they are all empty, we are empty.
    return true;
}

bool FMDPositionGroup::isUnique() const {
    for(const auto& annotated : positions) {
        // For each interval we contain
        if(view.isAmbiguous(annotated.position)) {
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
    
    for(const auto& annotated : positions) {
        // For each position
        auto position = annotated.position;
        
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

    for(const auto& annotated : positions) {
        // For each interval we contain
        if(view.isAmbiguous(annotated.position)) {
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
    
    for(const auto& annotated : positions) {
        // For each position
        auto position = annotated.position;
        
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















