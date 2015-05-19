#include "FMDPositionGroup.hpp"

#include "util.hpp"

FMDPositionGroup::FMDPositionGroup(): positions() {
    // Nothing to do!
}

FMDPositionGroup::FMDPositionGroup(
    const std::vector<FMDPosition>& startPositions): FMDPositionGroup() {
    
    for(const auto& position : startPositions) {
        // We need to tag all the positions with empty annotations
        positions.emplace(position);
    }
}

bool FMDPositionGroup::extendGreedy(const FMDIndexView& view, 
    char correctCharacter, size_t maxMismatches) {

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
        
        if(!toExtend.isEmpty(view)) {
            // If we got anything, keep the range (and don't increment the
            // mismatches)
            nonemptyExtensions.emplace(toExtend, annotated, false);
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
                
                if(!toExtend.isEmpty(view)) {
                    // If we got anything, keep the range (and charge a
                    // mismatch)
                    nonemptyExtensions.emplace(toExtend, annotated, true);
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


void FMDPositionGroup::retractOne(const FMDIndexView& view) {
    // Make a new collection to hold everything after we retract.
    decltype(positions) retracted;

    for(auto& annotated : positions) {
        if(annotated.searchedCharacters == 0) {
            // We can't retract if it's 0 length (or we don't know how long it
            // is.
            throw std::runtime_error("Cannot retract 0-length search.");
        }
        
        // Make a copy of the range to retract.
        FMDPosition toRetract = annotated.position;
        
        // Retract exactly one base (assuming we started with a covering
        // FMDPosition).
        view.getIndex().retractRightOnly(toRetract,
            annotated.searchedCharacters - 1);
            
        // Call the retraction constructor and stick in this new retracted
        // range.
        retracted.emplace(toRetract, annotated, (size_t) 1);
        
    }
    
    // Replace our old positions with the new retracted ones.
    positions = std::move(retracted);
    
}

bool FMDPositionGroup::isEmpty(const FMDIndexView& view) const {
    for(const auto& annotated : positions) {
        // For each interval we contain
        if(!annotated.position.isEmpty(view)) {
            // If any is nonempty, we are nonempty
            return false;
        }
    }
    
    // If none exists or they are all empty, we are empty.
    return true;
}

bool FMDPositionGroup::isUnique(const FMDIndexView& view) const {
    for(const auto& annotated : positions) {
        // For each interval we contain
        if(annotated.position.isAmbiguous(view)) {
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
        
        if(position.isEmpty(view)) {
            continue;
        }
        
        if(found) {
            // We already have a candidate unique place, so make sure this one
            // matches.
            if(uniquePlace != position.getTextPosition(view)) {
                return false;
            }
        } else {
            // This is the first one we found
            
            // Now we know it's nonempty and unique
            found = true;
            // Figure out where to
            uniquePlace = position.getTextPosition(view);
        }
    }
    
    // If we found something we're unique. Otherwise we're empty.
    return found;
}

TextPosition FMDPositionGroup::getTextPosition(const FMDIndexView& view) const {
    // TODO: Unify with above function, since they differ only in what they
    // return.

    for(const auto& annotated : positions) {
        // For each interval we contain
        if(annotated.position.isAmbiguous(view)) {
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
        
        if(position.isEmpty(view)) {
            continue;
        }
        
        if(found) {
            // We already have a candidate unique place, so make sure this one
            // matches.
            if(uniquePlace != position.getTextPosition(view)) {
                 throw std::runtime_error(
                    "No TextPosition for ambiguous FMDPositionGroup");
            }
        } else {
            // This is the first one we found
            
            // Now we know it's nonempty and unique
            found = true;
            // Figure out where to
            uniquePlace = position.getTextPosition(view);
        }
    }
    
    // If we found something we're unique. Otherwise we're empty.
    return uniquePlace;
}















