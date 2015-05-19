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

FMDPositionGroup::FMDPositionGroup(
    const FMDIndexView& view): FMDPositionGroup() {
    
    // Start with just one FMDPosition, covering the whole thing.
    positions.emplace(view);
    
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
                toExtend.extendLeftOnly(view, base);
                
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

bool FMDPositionGroup::extendFull(const FMDIndexView& view,
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
    
    // If we found anything, we have an exact match. But keep going with the
    // inexact ones regardless.
    exactMatch = (nonemptyExtensions.size() != 0);
    
        
    
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
            toExtend.extendLeftOnly(view, base);
            
            if(!toExtend.isEmpty(view)) {
                // If we got anything, keep the range (and charge a
                // mismatch)
                nonemptyExtensions.emplace(toExtend, annotated, true);
            }
        }
        
    }
    
    Log::debug() << "Extended with all but " << correctCharacter <<
        std::endl;
    
    // Replace our FMDPositions with the new extended ones.
    positions = std::move(nonemptyExtensions);
    
    Log::debug() << "Have " << positions.size() << " ranges" << std::endl;
    
    // Let the caller know if we found the base they wanted or not.
    return exactMatch;
}


void FMDPositionGroup::extendLeftOnly(const FMDIndexView& view,
    char character) {
    
    // This will hold all the extensions we find that aren't empty
    decltype(positions) nonemptyExtensions;
    
    for(const auto& annotated : positions) {
        // For each existing FMDPosition
        
        // Pull it out
        FMDPosition toExtend = annotated.position;
        
        // Extend it
        toExtend.extendLeftOnly(view, character);
        
        if(!toExtend.isEmpty(view)) {
            // If we got anything, keep the range (and don't increment the
            // mismatches)
            nonemptyExtensions.emplace(toExtend, annotated, false);
        }
    }
    
    // Replace our FMDPositions with the new extended ones.
    positions = std::move(nonemptyExtensions);
    
}

void FMDPositionGroup::retractRightOnly(const FMDIndexView& view,
    size_t newLength) {
    
    // This will hold all the retractions we find
    decltype(positions) retractions;
    
    for(const auto& annotated : positions) {
        // For each existing FMDPosition
        
        // Pull it out
        FMDPosition toRetract = annotated.position;
        
        // Retract it
        toRetract.retractRightOnly(view, newLength);
        
        // Call the retraction constructor with the number of bases we
        // retracted, and stick in this new retracted range.
        retractions.emplace(toRetract, annotated, 
            annotated.searchedCharacters - newLength);
    }
    
    // Replace our FMDPositions with the new extended ones.
    positions = std::move(retractions);
}

size_t FMDPositionGroup::retractRightOnly(const FMDIndexView& view) {

    // What length do we need to retract to to get something new?
    size_t newLength = 0;
    
    for(const auto& annotated : positions) {
        // For each existing FMDPosition
        
        // Pull it out
        FMDPosition toRetract = annotated.position;
        
        // See how far it needs to retract, and keep it if it stays longer.
        newLength = std::max(newLength, toRetract.retractRightOnly(view));
    }

    if(newLength == 0) {
        // We couldn't get a number out of anything. Maybe we ran out of
        // FMDPositions?
        throw std::runtime_error("No retraction possible");
    }

    // Retract to the longest length that lets us cover anything more.
    // TODO: We're half as fast as we could be, since we repeat a lot of work.
    retractRightOnly(view, newLength);
    
    // Report how logn we retracted to.
    return newLength;
}

size_t FMDPositionGroup::mismatchesUsed() const {
    // What't she max mismatches we've observed used?
    size_t used = 0;
    
    for(const auto& annotated : positions) {
        // Look through all the FMDPositions and get the max number of
        // mismatches.
        used = std::max(used, annotated.mismatches);
    }
    
    return used;
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

size_t FMDPositionGroup::getApproximateNumberOfRanges(
    const FMDIndexView& view) const {
    
    size_t total = 0;
    
    for(const auto& annotated : positions) {
        // Just sum up over all the intervals we contain.
        total += annotated.position.getApproximateNumberOfRanges(view);
    }
    
    return total;
}

size_t FMDPositionGroup::getApproximateNumberOfNewRanges(
    const FMDIndexView& view, const FMDPositionGroup& old) {

    // How many have we found across all our positions in the group?
    size_t total = 0;

    // We know intervals won't overlap in either group.
    
    // So we can just get the parent starting at or after the start of the
    // widened interval, using map::lowert_bound.

    // Index the parents
    std::map<size_t, FMDPosition> oldPositions;
    
    for(auto parentAnnotated : old.positions) {
        // For each position in the old group, pull it out.
        const auto& parent = parentAnnotated.position;
        
        // Put it in the index. Since the FMDPositions won't overlap, we will
        // never overwrite another one like this.
        oldPositions[parent.getForwardStart()] = parent;
    }
    
    for(auto newAnnotated : positions) {
        // For each new position, find the old position that started at or after
        // this one's start.
        auto parentIterator = oldPositions.lower_bound(
            newAnnotated.position.getForwardStart());
            
        if(parentIterator == oldPositions.end()) {
            // Complain if we can't find where a range came from.
            throw std::runtime_error(
                "Could not find parent for expanded interval");
        }
        
        // Go work out how many new ranges this expanded range has on top of
        // this parent. TODO: there may be some inaccuracy if it had multiple
        // parents. Handle that case and provide a tighter estimate.
        total += newAnnotated.position.getApproximateNumberOfNewRanges(view,
            (*parentIterator).second);
    }
    
    return total;


}
    
std::set<TextPosition> FMDPositionGroup::getTextPositions(
    const FMDIndexView& view) const {

    // Just loop over all the positions in the group and collect in here.
    std::set<TextPosition> toReturn;
    
    for(const auto& annotated : positions) {
        // Grab the resulkts from everything in the group.
        std::set<TextPosition> found = annotated.position.getTextPositions(
            view);
            
        // And stick them in the set. See
        // <http://bytes.com/topic/c/answers/718481-better-prettier-way-copy-
        // one-set-another> if you are like me and don't know how to STL at all.
        std::copy(found.begin(), found.end(), std::inserter(toReturn,
            toReturn.begin()));
    }
    
    return toReturn;
}

std::set<TextPosition> FMDPositionGroup::getNewTextPositions(
    const FMDIndexView& view, const FMDPositionGroup& old) const {
    
    // We know intervals won't overlap in either group.
    
    // So we can just get the parent starting at or after the start of the
    // widened interval, using map::lowert_bound.
    
    // Just loop over all the positions in the group and collect in here.
    std::set<TextPosition> toReturn;

    // Index the parents
    std::map<size_t, FMDPosition> oldPositions;
    
    for(auto parentAnnotated : old.positions) {
        // For each position in the old group, pull it out.
        const auto& parent = parentAnnotated.position;
        
        // Put it in the index. Since the FMDPositions won't overlap, we will
        // never overwrite another one like this.
        oldPositions[parent.getForwardStart()] = parent;
    }
    
    for(auto newAnnotated : positions) {
        // For each new position, find the old position that started at or after
        // this one's start.
        auto parentIterator = oldPositions.lower_bound(
            newAnnotated.position.getForwardStart());
            
        if(parentIterator == oldPositions.end()) {
            // Complain if we can't find where a range came from.
            throw std::runtime_error(
                "Could not find parent for expanded interval");
        }
        
        // Go work what new ranges this expanded range has on top of
        // this parent. TODO: there may be some inaccuracy if it had multiple
        // parents. Handle that case and provide a tighter estimate.
        std::set<TextPosition> found = 
            newAnnotated.position.getNewTextPositions(view,
                (*parentIterator).second);
            
        // And stick them in the set. See
        // <http://bytes.com/topic/c/answers/718481-better-prettier-way-copy-
        // one-set-another> if you are like me and don't know how to STL at all.
        std::copy(found.begin(), found.end(), std::inserter(toReturn,
            toReturn.begin()));
    }
    
    return toReturn;
    
}















