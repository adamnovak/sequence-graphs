#include "ZipMappingScheme.hpp"
#include "Log.hpp"
#include "util.hpp"

#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

#include <boost/range/adaptor/reversed.hpp>

// Only template specializations can really be here.

template<>
bool ZipMappingScheme<FMDPosition>::canExtendThrough(FMDPosition context,
    const std::string& opposingQuery) const {
    
    Log::debug() << "Trying to extend " << context << " through " << 
        opposingQuery.size() << " opposing context" << std::endl;
    stats.add("extendThroughAttampts", 1);
        
    // We're going to retract it until it's no longer unique, then go
    // back and retract it one less.
    FMDPosition noLongerUnique = context;
    size_t nonUniqueLength = noLongerUnique.retractRightOnly(view);
    
    while(noLongerUnique.isUnique(view)) {
        // It may take multiple retracts because of masks and stuff.
        nonUniqueLength = noLongerUnique.retractRightOnly(view);
    }
    
    // Go back and retract one less (to a length one longer).
    FMDPosition barelyUnique = context;
    barelyUnique.retractRightOnly(view, nonUniqueLength + 1);

    for(size_t i = opposingQuery.size() - 2; i != (size_t) -1 &&
        !barelyUnique.isEmpty(view); i--) {
        // Now go extend through with the opposing string, from right to left.
        // We skip the rightmost character (the one we are actually in the
        // process of mapping) because that's already searched, and we keep
        // going until we run out of results or we make it all the way through.
        
        Log::trace() << "Extending with " << opposingQuery[i] << std::endl;
        
        barelyUnique.extendLeftOnly(view, opposingQuery[i]);
    }
    
    if(!barelyUnique.isEmpty(view)) {
        // We extended through!
        stats.add("extendThroughSuccesses", 1);
        return true;
    } else {
        // We didn't get any results upon extending through
        Log::debug() << "Extension failed" << std::endl;
        return false;
    }
}

template<>
bool ZipMappingScheme<FMDPositionGroup>::canExtendThrough(
    FMDPositionGroup context, const std::string& opposingQuery) const {
    
    // Now we need to make sure to allow the right number of mismatches.
    
    Log::debug() << "Trying to extend " << context << " through " << 
        opposingQuery.size() << " opposing context" << std::endl;
    stats.add("extendThroughAttampts", 1);
        
    // We're going to retract it until it's no longer unique, then go
    // back and retract it one less.
    FMDPositionGroup noLongerUnique = context;
    size_t nonUniqueLength = noLongerUnique.retractRightOnly(view);
    
    while(noLongerUnique.isUnique(view)) {
        // It may take multiple retracts because of masks and stuff.
        nonUniqueLength = noLongerUnique.retractRightOnly(view);
    }
    
    // Go back and retract one less (to a length one longer).
    FMDPositionGroup barelyUnique = context;
    barelyUnique.retractRightOnly(view, nonUniqueLength + 1);

    for(size_t i = opposingQuery.size() - 2; i != (size_t) -1 &&
        !barelyUnique.isEmpty(view); i--) {
        // Now go extend through with the opposing string, from right to left.
        // We skip the rightmost character (the one we are actually in the
        // process of mapping) because that's already searched, and we keep
        // going until we run out of results or we make it all the way through.
        
        Log::trace() << "Extending with " << opposingQuery[i] <<
            " up to " << mismatchTolerance << " mismatches" << std::endl;
        
        // Do the extension allowing for mismatches
        barelyUnique.extendFull(view, opposingQuery[i], mismatchTolerance);
    }
    
    if(!barelyUnique.isEmpty(view)) {
        // We extended through!
        stats.add("extendThroughSuccesses", 1);
        return true;
    } else {
        // We didn't get any results upon extending through
        Log::debug() << "Extension failed" << std::endl;
        return false;
    }
}

template<>
std::pair<bool, std::set<TextPosition>>
    ZipMappingScheme<FMDPosition>::exploreRetraction(
    const typename DPTable::DPTask& task, DPTable& table,
    const std::string& query, size_t queryBase) const {
    
    // Get the left and right SideRetractionEntries. Also computes any
    // uncomputed retractions and enumerates any sufficiently small sets.
    const typename DPTable::SideRetractionEntry& left = table.getRetraction(
        false, task.leftIndex, view, maxRangeCount);
    const typename DPTable::SideRetractionEntry& right = table.getRetraction(
        true, task.rightIndex, view, maxRangeCount);
    
    Log::debug() << "Exploring " << left.selection << " (" <<
        left.contextLength << ") left, " << right.selection << " (" <<
        right.contextLength << ") right" << std::endl;
    
    // How many bases are searched total (accounting for the 1-base central
    // overlap)?
    size_t totalContext = left.contextLength + right.contextLength - 1;    
    if(totalContext < minContextLength) {
        // If we have too little total context, we would ignore any results we
        // found, so say we have no results but should still map.
        Log::debug() << "Skipping task due to too little context." <<
            std::endl;
            
        return {true, {}};
    }
          
    // We keep track of whether we tried to extend the left through the right,
    // or the right through the left. If we did, we know we will try the same
    // things if we retract on the opposite side, and we will queue up those
    // retractions even if we can't afford to correctly populate their sets.
    bool triedLeftExtend = false;
    bool triedRightExtend = false;
          
    
    if(left.selection.isUnique(view) && right.contextLength < maxExtendThrough) {
        
        // We will try extending the left through the right.
        triedLeftExtend = true;
        
        if(canExtendThrough(left.selection, reverseComplement(query.substr(
            queryBase, right.contextLength)))) {
            // We tried to extend the left context through the (reverse
            // complemented) right context and succeeded.
            
            Log::debug() << "Successfully extended left through right" <<
                std::endl;
            
            // We can now cheat and return the one-element set of whatever the
            // left context selects, flipped.
            TextPosition matchedTo = left.selection.getTextPosition(view);
            matchedTo.flip(view.getIndex().getContigLength(
                matchedTo.getContigNumber()));
            return {true, {matchedTo}};
        } else {
            Log::debug() << "Failed to extend left though right" << std::endl;
        }
        
        // If we don't succeed, we still need to check for set overlap.
    }
    if(right.selection.isUnique(view) && left.contextLength < maxExtendThrough) {
        // TODO: Don't try this extension if we failed the other direction.
        
        // We will try extending the right through the left.
        triedRightExtend = true;
        
        if(canExtendThrough(right.selection, query.substr(
            queryBase - (left.contextLength - 1), left.contextLength))) {
            // We tried to extend the right context through the left context and
            // succeeded
            
            Log::debug() << "Successfully extended right through left" <<
                std::endl;
            
            // We can now cheat and return the one-element set of whatever the
            // right context selects.
            TextPosition matchedTo = right.selection.getTextPosition(view);
            return {true, {matchedTo}};
        } else {
            Log::debug() << "Failed to extend right though left" << std::endl;
        }
        
        // If we don't succeed, we still need to check for set overlap.
    }
    
    // If we can't or won't try to extend one through the other, can we
    // successfully bang the sets together?
    
    // Can we do the set banging? Or would one of the sets have been too big?
    bool canTouchResults = left.setsValid && right.setsValid;
    
    if(canTouchResults) {
        // We have properly populated sets for both sides.
        
        // Note that root tasks will have equal selected and newlySelected sets
        // on each side, so we can handle them fairly simply.
        
        // If we retracted on the right, use the old positions from the left.
        // Otherwise (if we retracted on the left or are a root), use the old
        // positions from the right
        const std::set<TextPosition>& oldPositions = task.retractedRight ? 
            left.selected : right.selected;
            
        // If we retracted on the right, bang the new right positions against
        // the old positions. Otherwise bang the new left positions.
        const std::set<TextPosition>& newPositions = task.retractedRight ?
            right.newlySelected : left.newlySelected;
            
        Log::debug() << "Banging " << newPositions.size() <<
            " new positions against " << oldPositions.size() << " old ones" <<
            std::endl;
        
        // This is going to hold the TextPositiuons in both sets, in right
        // orientation.
        std::set<TextPosition> shared;
        
        // TODO: scan the smaller set always?
        for(auto result : newPositions) {
                
            // Flip each result around and see if it is in the opposing set
            TextPosition flipped = result;
            flipped.flip(view.getIndex().getContigLength(
                flipped.getContigNumber()));
                
            if(oldPositions.count(flipped)) {
                // We found it!
                if(task.retractedRight) {
                    // We're going through right contexts.
                    shared.insert(result);
                } else {
                    // We're going through left contexts, so insert the flipped
                    // version.
                    shared.insert(flipped);
                }
                
                if(shared.size() > 1) {
                    
                    Log::debug() << "Found multiple shared new results" <<
                        std::endl;
                
                    // We can short circuit now because we found multiple
                    // shared TextPositions.
                    return {true, shared};
                }
            }
        }
        
        // If we get here, we went through all the new results.
        if(shared.size() > 0) {
            // We found results, but didn't terminate the loop early due to
            // having 2 or more. So there must be just 1.
            
            Log::debug() << "Found " << shared.size() <<
                " shared new results" << std::endl;
        
            // If we found any overlap at all, say we found something.
            return {true, shared};
        }
    
        Log::debug() << "Could look at result sets but found no overlap" <<
            std::endl;
    } else {
        Log::debug() << "Result set(s) would be too big!" << std::endl;
    }
    
    // Make a task that we can retract to make child tasks, if needed.
    typename DPTable::DPTask toRetract = task;
    
    if(left.contextLength > 1 && (canTouchResults || triedRightExtend)) {
        // We have stuff on the left, and either we're operating normally or we
        // just tried to extend the right through the left and failed.
        
        // Retract on the left and enqueue that child task.
        Log::debug() << "Queueing retraction on the left" << std::endl;
        table.taskQueue.push(toRetract.retract(false));
    }
    
    if(((task.isRightEdge && canTouchResults) || triedLeftExtend) &&
        right.contextLength > 1) {
        
        // Either we're on the right edge (haven't retracted left yet) and we
        // can look at our results correctly, or we just tried and failed
        // extending the left through the right (and want to try with a shorter
        // right). TODO: Am I going to duplicate work here? I think I really do
        // need to approach some retractions from this side to make things work.
        
        Log::debug() << "Queueing retraction on the right" << std::endl;
        // Also retract on the right, if the task we just did never retracted on
        // the left.
        table.taskQueue.push(toRetract.retract(true));
    }
    
    Log::debug() << "Found no shared results" << std::endl;
    
    // Report whether we got too many results to think about, and say we found
    // no overlaps.
    return {canTouchResults, {}};
        
}


template<>
std::pair<bool, std::set<TextPosition>>
    ZipMappingScheme<FMDPositionGroup>::exploreRetraction(
    const typename DPTable::DPTask& task, DPTable& table,
    const std::string& query, size_t queryBase) const {
    
    // Get the left and right SideRetractionEntries. Also computes any
    // uncomputed retractions and enumerates any sufficiently small sets.
    const typename DPTable::SideRetractionEntry& left = table.getRetraction(
        false, task.leftIndex, view, maxRangeCount);
    const typename DPTable::SideRetractionEntry& right = table.getRetraction(
        true, task.rightIndex, view, maxRangeCount);
        
    Log::debug() << "Exploring " << left.selection << " (" <<
        left.contextLength << ") left, " << right.selection << " (" <<
        right.contextLength << ") right" << std::endl;
    
    // We need to look at the searches in both directions at this level of
    // retraction, and determine if they can agree. They can only agree if they
    // overlap on some BWT position, and they don't have too many mismatches. If
    // they overlap on a certain BWT position, each direction includes that
    // position with a certain number of mismatches. However, that can change
    // with retraction, and a position could be included with fewer mismatches.
    // So we can't just do what we have been doing and only look at the stuff we
    // gained by retraction; we also have to look at the stuff we had already
    // but which is now available with fewer mismatches.
    
    // And we have to know that at this smaller number of mismatches, the
    // searches don't overlap at two places.
    
    // How many bases are searched total (accounting for the 1-base central
    // overlap)?
    size_t totalContext = left.contextLength + right.contextLength - 1;    
    if(totalContext < minContextLength) {
        // If we have too little total context, we would ignore any results we
        // found, so say we have no results but should still map.
        Log::debug() << "Skipping task due to too little context." <<
            std::endl;
            
        return {true, {}};
    }
          
    // We keep track of whether we tried to extend the left through the right,
    // or the right through the left. If we did, we know we will try the same
    // things if we retract on the opposite side, and we will queue up those
    // retractions even if we can't afford to correctly populate their sets.
    bool triedLeftExtend = false;
    bool triedRightExtend = false;
          
    
    if(left.selection.isUnique(view) && right.contextLength < maxExtendThrough) {
        
        // We will try extending the left through the right.
        triedLeftExtend = true;
        
        if(canExtendThrough(left.selection, reverseComplement(query.substr(
            queryBase, right.contextLength)))) {
            // We tried to extend the left context through the (reverse
            // complemented) right context and succeeded. This accounts for
            // total mismatches.
            
            Log::debug() << "Successfully extended left through right" <<
                std::endl;
            
            // We can now cheat and return the one-element set of whatever the
            // left context selects, flipped.
            TextPosition matchedTo = left.selection.getTextPosition(view);
            matchedTo.flip(view.getIndex().getContigLength(
                matchedTo.getContigNumber()));
            return {true, {matchedTo}};
        } else {
            Log::debug() << "Failed to extend left though right" << std::endl;
        }
        
        // If we don't succeed, we still need to check for set overlap.
    }
    if(right.selection.isUnique(view) && left.contextLength < maxExtendThrough) {
        // TODO: Don't try this extension if we failed the other direction.
        
        // We will try extending the right through the left.
        triedRightExtend = true;
        
        if(canExtendThrough(right.selection, query.substr(
            queryBase - (left.contextLength - 1), left.contextLength))) {
            // We tried to extend the right context through the left context and
            // succeeded. This accounts for total mismatches.
            
            Log::debug() << "Successfully extended right through left" <<
                std::endl;
            
            // We can now cheat and return the one-element set of whatever the
            // right context selects.
            TextPosition matchedTo = right.selection.getTextPosition(view);
            return {true, {matchedTo}};
        } else {
            Log::debug() << "Failed to extend right though left" << std::endl;
        }
        
        // If we don't succeed, we still need to check for set overlap.
    }
    
    // If we can't or won't try to extend one through the other, can we
    // successfully bang the sets together?
    
    // Can we do the set banging? Or would one of the sets have been too big?
    bool canTouchResults = left.setsValid && right.setsValid;
    
    if(canTouchResults) {
        // We have properly populated sets for both sides.
        
        // How many mismatches are used by this combined context? We know there
        // won't be a mismatch at the base in question, so we don't have to
        // worry about double-counting there.
        size_t mismatchesUsed = left.selection.mismatchesUsed() +
            right.selection.mismatchesUsed();
        
        if(mismatchesUsed <= mismatchTolerance) {
        
            // How many mismatches did we use in our parent?
            // TODO: Just remember this.
            size_t lastMismatchesUsed = (size_t) -1;
                
            if(task.retractedRight && task.rightIndex > 0) {
                // Our parent had this left and a different right. Figure out
                // how many mismatches it used.
                lastMismatchesUsed = left.selection.mismatchesUsed() + 
                    table.getRetraction(true, task.rightIndex - 1, view,
                    maxRangeCount).selection.mismatchesUsed();
            } else if(!task.retractedRight && task.leftIndex > 0) {
                // Our parent had this right and a different left. Figure out
                // how many mismatches it used.
                lastMismatchesUsed = right.selection.mismatchesUsed() + 
                    table.getRetraction(false, task.leftIndex - 1, view,
                    maxRangeCount).selection.mismatchesUsed();
            }
        
            // This is going to hold the TextPositions in both sets, in right
            // orientation.
            std::set<TextPosition> shared;
        
            if(lastMismatchesUsed > mismatchTolerance) {
                // This is the first retraction that has few enough mismatches.
        
                if(lastMismatchesUsed == -1) {
                    Log::debug() << "Only " << mismatchesUsed << 
                        " mismatches; comparing entire result sets." <<
                        std::endl;
                } else {
                    Log::debug() << "Only " << mismatchesUsed << 
                        " mismatches, down from from " << lastMismatchesUsed << 
                        " mismatches; comparing entire result sets." <<
                        std::endl;
                }
            
                // We can accept overlaps because not too many mismatches are
                // used.
                
                // TODO: For now we can just bang the selected sets against each
                // other, and ignore the newly selected sets. In the future we
                // need to be able to work out if we just dropped enough
                // mismatches.
                
                // TODO: scan the smaller set always?
                for(auto result : left.selected) {
                        
                    // Flip each result around and see if it is in the opposing
                    // set
                    TextPosition flipped = result;
                    flipped.flip(view.getIndex().getContigLength(
                        flipped.getContigNumber()));
                        
                    if(right.selected.count(flipped)) {
                        // We found it!
                        
                        // We're going through left contexts, so insert the
                        // flipped version.
                        shared.insert(flipped);
                        
                        if(shared.size() > 1) {
                            
                            Log::debug() <<
                                "Found multiple shared new results" <<
                                std::endl;
                        
                            // We can short circuit now because we found
                            // multiple shared TextPositions.
                            return {true, shared};
                        }
                    }
                }
            } else {
                // We have few enough mismatches, but we did before, too. Only
                // check the new results.
                
                // TODO: share code
                
                // Note that root tasks will have equal selected and newlySelected sets
                // on each side, so we can handle them fairly simply.
                
                // If we retracted on the right, use the old positions from the
                // left. Otherwise (if we retracted on the left or are a root),
                // use the old positions from the right
                const std::set<TextPosition>& oldPositions = 
                    task.retractedRight ? left.selected : right.selected;
                    
                // If we retracted on the right, bang the new right positions
                // against the old positions. Otherwise bang the new left
                // positions.
                const std::set<TextPosition>& newPositions =
                    task.retractedRight ? right.newlySelected :
                    left.newlySelected;
                    
                Log::debug() << "Comparing " << newPositions.size() <<
                    " new positions against " << oldPositions.size() << 
                    " old ones" << std::endl;
                
                // TODO: scan the smaller set always?
                for(auto result : newPositions) {
                        
                    // Flip each result around and see if it is in the opposing
                    // set
                    TextPosition flipped = result;
                    flipped.flip(view.getIndex().getContigLength(
                        flipped.getContigNumber()));
                        
                    if(oldPositions.count(flipped)) {
                        // We found it!
                        if(task.retractedRight) {
                            // We're going through right contexts.
                            shared.insert(result);
                        } else {
                            // We're going through left contexts, so insert the
                            // flipped version.
                            shared.insert(flipped);
                        }
                        
                        if(shared.size() > 1) {
                            
                            Log::debug() <<
                                "Found multiple shared new results" <<
                                std::endl;
                        
                            // We can short circuit now because we found
                            // multiple shared TextPositions.
                            return {true, shared};
                        }
                    }
                }
            }
            
            // If we get here, we went through all the results.
            if(shared.size() > 0) {
                // We found results, but didn't terminate the loop early due to
                // having 2 or more. So there must be just 1.
                
                Log::debug() << "Found " << shared.size() <<
                    " shared results" << std::endl;
            
                // If we found any overlap at all, say we found something.
                return {true, shared};
            }
        
            Log::debug() << "Found no overlap" << std::endl;
        } else {
            // We used too many mismatches to accept an overlap. Try after
            // retracting.
            Log::debug() << "Skipping for too many mismatches: " <<
                mismatchesUsed << std::endl;
                
            // We need to retract more.
        }
        
    } else {
        Log::debug() << "Result set(s) would be too big!" << std::endl;
    }
    
    // Make a task that we can retract to make child tasks, if needed.
    typename DPTable::DPTask toRetract = task;
    
    if(left.contextLength > 1 && (canTouchResults || triedRightExtend)) {
        // We have stuff on the left, and either we're operating normally or we
        // just tried to extend the right through the left and failed.
        
        // Retract on the left and enqueue that child task.
        Log::debug() << "Queueing retraction on the left" << std::endl;
        table.taskQueue.push(toRetract.retract(false));
    }
    
    if(((task.isRightEdge && canTouchResults) || triedLeftExtend) &&
        right.contextLength > 1) {
        
        // Either we're on the right edge (haven't retracted left yet) and we
        // can look at our results correctly, or we just tried and failed
        // extending the left through the right (and want to try with a shorter
        // right). TODO: Am I going to duplicate work here? I think I really do
        // need to approach some retractions from this side to make things work.
        
        Log::debug() << "Queueing retraction on the right" << std::endl;
        // Also retract on the right, if the task we just did never retracted on
        // the left.
        table.taskQueue.push(toRetract.retract(true));
    }
    
    Log::debug() << "Found no shared results" << std::endl;
    
    // Report whether we got too many results to think about, and say we found
    // no overlaps.
    return {canTouchResults, {}};
        
}

template<>
std::vector<std::pair<FMDPosition, size_t>>
    ZipMappingScheme<FMDPosition>::findRightContexts(const std::string& query,
    bool reverse) const {
 
    // This will hold the right context of every position.
    std::vector<std::pair<FMDPosition, size_t>> toReturn(query.size());
 
    // Start with everything selected.
    FMDPosition results(view);
    
    // How many characters are currently searched?
    size_t patternLength = 0;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--) {
        // For each position in the query from right to left, we're going to
        // inchworm along and get the search that is extended out right as
        // far as possible while still having results.
        
        // We're going to extend left with this new base.
        FMDPosition extended = results;
        extended.extendLeftOnly(view, query[i]);
        
        while(extended.isEmpty(view)) {
            // If you couldn't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPosition retracted = results;
            // Make sure to drop characters from the total pattern length.
            retracted.retractRightOnly(view, --patternLength);
            
            // Try extending again
            extended = retracted;
            extended.extendLeftOnly(view, query[i]);
            
            // Say that last step we came from retracted.
            results = retracted;
        }
    
        // We successfully extended (if only with the base itself).
        results = extended;
        // Increment the pattern length since we did actually extedn by 1. 
        patternLength++;
        
        Log::trace() << "Index " << i << " has " << query[i] << " + " << 
            patternLength - 1 << " selecting " << results << std::endl;
        
        // Save the search results to the appropriate location, depending on if
        // we want to reverse the results or not.
        toReturn[reverse ? toReturn.size() - i - 1 : i] = std::make_pair(
            results, patternLength);
    }

    // Now we have inchwormed all the way from right to left, retracting only
    // when necessary. Return all the results.    
    return toReturn;
}

template<>
std::vector<std::pair<FMDPositionGroup, size_t>>
    ZipMappingScheme<FMDPositionGroup>::findRightContexts(
    const std::string& query, bool reverse) const {
 
    // This will hold the right context of every position.
    std::vector<std::pair<FMDPositionGroup, size_t>> toReturn(query.size());
 
    // Start with everything selected.
    FMDPositionGroup results(view);
    
    // How many characters are currently searched?
    size_t patternLength = 0;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--) {
        // For each position in the query from right to left, we're going to
        // inchworm along and get the search that is extended out right as
        // far as possible while still having results.
        
        // We're going to extend left with this new base.
        FMDPositionGroup extended = results;
        extended.extendFull(view, query[i], mismatchTolerance);
        
        while(extended.isEmpty(view)) {
            // If you couldn't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPositionGroup retracted = results;
            // Make sure to drop characters from the total pattern length.
            retracted.retractRightOnly(view, --patternLength);
            
            // Try extending again
            extended = retracted;
            extended.extendFull(view, query[i], mismatchTolerance);
            
            // Say that last step we came from retracted.
            results = retracted;
        }
        
        // This is what we want to use on the next base. But for this base we
        // have an additional constraint of no mismatch here.
        FMDPositionGroup carryForwardResults = extended;
        size_t carryForwardPatternLength = patternLength + 1;
        
        // Apply that constraint
        extended.dropMismatchesHere();
        
        while(extended.isEmpty(view)) {
            // If you couldn't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPositionGroup retracted = results;
            // Make sure to drop characters from the total pattern length.
            retracted.retractRightOnly(view, --patternLength);
            
            // Try extending again
            extended = retracted;
            extended.extendFull(view, query[i], mismatchTolerance);
            extended.dropMismatchesHere();
            
            // Say that last step we came from retracted.
            results = retracted;
        }
    
        // We successfully extended (if only with the base itself).
        results = extended;
        // Increment the pattern length since we did actually extend by 1. 
        patternLength++;
        
        Log::trace() << "Index " << i << " has " << query[i] << " + " << 
            patternLength - 1 << " selecting " << results << std::endl;
        
        // Save the search results to the appropriate location, depending on if
        // we want to reverse the results or not.
        toReturn[reverse ? toReturn.size() - i - 1 : i] = std::make_pair(
            results, patternLength);
            
        // Now wind back to before applying the constraint. TODO: this is
        // awkward.
        results = carryForwardResults;
        patternLength = carryForwardPatternLength;
    }

    // Now we have inchwormed all the way from right to left, retracting only
    // when necessary. Return all the results.    
    return toReturn;
}




