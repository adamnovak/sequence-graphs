#include "ZipMappingScheme.hpp"
#include "Log.hpp"
#include "util.hpp"

#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

#include <boost/range/adaptor/reversed.hpp>

std::vector<std::pair<FMDPosition, size_t>> ZipMappingScheme::findRightContexts(
    const std::string& query) const {
 
    // This will hold the right context of every position.
    std::vector<std::pair<FMDPosition, size_t>> toReturn(query.size());
 
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
        
        Log::trace() << "Index " << i << " has " << query[i] << " + " << 
            patternLength - 1 << " selecting " << results << std::endl;
        
        // Save the search results.
        toReturn[i] = std::make_pair(results, patternLength);
    }

    // Now we have inchwormed all the way from right to left, retracting only
    // when necessary. Return all the results.    
    return toReturn;
}

void ZipMappingScheme::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // Get the right contexts
    Log::info() << "Looking for right contexts" << std::endl;
    auto rightContexts = findRightContexts(query);
    
    // And the left contexts (which is the reverse of the contexts for the
    // reverse complement.
    Log::info() << "Looking for left contexts" << std::endl;
    auto leftContexts = findRightContexts(reverseComplement(query));
    std::reverse(leftContexts.begin(), leftContexts.end());
    
    // For each pair, figure out if the forward and reverse FMDPositions select
    // one single consistent TextPosition.
    for(size_t i = 0; i < query.size(); i++) {
    
        Log::debug() << "Base " << i << " = " << query[i] << " (+" << 
            leftContexts[i].second << "|+" << rightContexts[i].second << 
            ") selects " << leftContexts[i].first << " and " << 
            rightContexts[i].first << std::endl;
            
        // Go look at these two contexts in opposite directions, and all their
        // (reasonable to think about) retractions, and see whether this base
        // belongs to 0, 1, or multiple TextPositions.
        std::set<TextPosition> intersection = exploreRetractions(
            leftContexts[i].first, leftContexts[i].second,
            rightContexts[i].first, rightContexts[i].second, query, i);
        
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

bool ZipMappingScheme::canExtendThrough(FMDPosition context,
    const std::string& opposingQuery) const {
    
    Log::debug() << "Trying to extend " << context << " through " << 
        opposingQuery.size() << " opposing context" << std::endl;
        
    // We're going to retract it until it's no longer unique, then go
    // back and retract it one less.
    FMDPosition noLongerUnique = context;
    size_t nonUniqueLength = view.getIndex().retractRightOnly(noLongerUnique);
    
    // Go back and retract one less (to a length one longer).
    FMDPosition barelyUnique = context;
    view.getIndex().retractRightOnly(barelyUnique, nonUniqueLength + 1);
    
    for(size_t i = 1; i < opposingQuery.size() && !view.isEmpty(barelyUnique);
        i++) {
        // Now go extend through with the opposing string. We skip the first
        // character (the one we are actually in the process of mapping) because
        // that's already searched, and we keep going until we run out of
        // results or we make it all the way through.
        
        Log::debug() << "Extending with " << opposingQuery[i] << std::endl;
        
        view.getIndex().extendLeftOnly(barelyUnique, opposingQuery[i]);
    }
    
    // Let the caller know whether we got all the way through (in which case the
    // LR context exists) or not.
    return !view.isEmpty(barelyUnique);
}

std::pair<bool, std::set<TextPosition>> ZipMappingScheme::exploreRetraction(
        const DPTask& task, std::queue<DPTask>& taskQueue, 
        const std::string& query, size_t queryBase) const {
    
    Log::debug() << "Exploring " << task.left << " (" << task.leftContext <<
        ") left, " << task.right << " (" << task.rightContext << ") right" <<
        std::endl;
    
    // How many bases are searched total (accounting for the 1-base central
    // overlap)?
    size_t totalContext = task.leftContext + task.rightContext - 1;    
    if(totalContext < minContextLength) {
        // If we have too little total context, we would ignore any results we
        // found, so say we have no results but should still map.
        Log::debug() << "Skipping task due to too little context." <<
            std::endl;
            
        return {true, {}};
    }
          
    
    if(view.isUnique(task.left) && view.isAmbiguous(task.right) && 
        task.rightContext < maxExtendThrough && canExtendThrough(task.left, 
        query.substr(queryBase, task.rightContext))) {
        // We tried to extend the left context through the right context and
        // succeeded. 
        
        Log::debug() << "Successfully extended left through right" << std::endl;
        
        // We can now cheat and return the one-element set of whatever the left
        // context selects, flipped.
        TextPosition matchedTo = view.getTextPosition(task.left);
        matchedTo.flip(view.getIndex().getContigLength(
            matchedTo.getContigNumber()));
        return {true, {matchedTo}};
    }
    
    if(view.isUnique(task.right) && view.isAmbiguous(task.left) && 
        task.leftContext < maxExtendThrough && canExtendThrough(task.right, 
        query.substr(queryBase - (task.leftContext - 1), task.leftContext))) {
        // We tried to extend the right context through the left context and
        // succeeded
        
        Log::debug() << "Successfully extended right through left" << std::endl;
        
        // We can now cheat and return the one-element set of whatever the right
        // context selects.
        TextPosition matchedTo = view.getTextPosition(task.right);
        return {true, {matchedTo}};
    }
    
    // If not, can we successfully bang the sets together?
    
    // What's the set of newly found results?
    std::set<TextPosition> newResults;
    
    // What's the set of old results they need to be merged into?
    std::set<TextPosition> oldResults;
    
    // What's the set of results we knew about on the opposing side?
    std::set<TextPosition> opposingResults;
    
    if(task.left == task.lastLeft && task.right == task.lastRight) {
        // We're a root task, we need the entirety of both sets.
        
        if(view.getApproximateNumberOfRanges(task.left) > maxRangeCount ||
            view.getApproximateNumberOfRanges(task.right) > maxRangeCount) {
            
            Log::debug() << "Can't hold all the ranges at the root" <<
                std::endl;
            
            // We can't actually evaluate this node because one side has too
            // many things.
            return {false, {}};
            
        }
        
        // We can touch all the results
        newResults = view.getTextPositions(task.right);
        oldResults = std::set<TextPosition>();
        opposingResults = view.getTextPositions(task.left);
    } else if(task.left == task.lastLeft) {
        // We retracted on the right and left the left alone
        
        if(view.getApproximateNumberOfNewRanges(task.lastRight, task.right) >
            maxRangeCount) {
            
            Log::debug() << "Can't hold all the ranges on the right" <<
                std::endl;
            
            // We found too many new things, we have to give up.
            return {false, {}};
        }
        
        // Our new results are just the new stuff we found
        newResults = view.getNewTextPositions(task.lastRight, task.right);
        // We pull our old already-known results from the task. TODO: Can we
        // move the set here?
        oldResults = task.lastRightPositions;
        opposingResults = task.lastLeftPositions;
        
    } else {
        // We retracted on the left and left the right alone
        
        if(view.getApproximateNumberOfNewRanges(task.lastLeft, task.left) >
            maxRangeCount) {
            
            Log::debug() << "Can't hold all the ranges on the left" <<
                std::endl;
            
            // We found too many new things, we have to give up.
            return {false, {}};
        }
        
        // Our new results are just the new stuff we found
        newResults = view.getNewTextPositions(task.lastLeft, task.left);
        // We pull our old already-known results from the task. TODO: Can we
        // move the set here?
        oldResults = task.lastLeftPositions;
        opposingResults = task.lastRightPositions;
    }
    
    // What TextPositions are shared? Always holds things from the right context
    // perspective.
    std::set<TextPosition> shared;
    
    for(auto result : newResults) {
        
        // Flip each result around and see if it is in the opposing set
        TextPosition flipped = result;
        flipped.flip(view.getIndex().getContigLength(
            flipped.getContigNumber()));
            
        if(opposingResults.count(flipped)) {
            // We found it!
            if(task.left == task.lastLeft) {
                // We didn't retract on the left, so we're going through right
                // contexts.
                shared.insert(result);
            } else {
                // We're going through left contexts, and we need to save
                // overlaps in right context space.
                shared.insert(flipped);
            }
            
            if(shared.size() > 1) {
                
                Log::debug() << "Found multiple shared new results" <<
                    std::endl;
            
                // We can short circuit now because we found multiple shared
                // TextPositions.
                return {true, shared};
            }
        }
        
        // If we're still going, make sure to add this new results to the
        // set with the old results from its side.
        oldResults.insert(result);
    }
    
    // If we get here, we went through all the new results.
    if(shared.size() > 0) {
    
        Log::debug() << "Found single shared new result" << std::endl;
    
        // If we found any overlap at all, say we found something.
        return {true, shared};
    }
    
    // Otherwise we didn't find anything and have to continue with our
    // descendants.
    
    // Make a task that reflects our new sets.
    DPTask toRetract = task;
    if(task.left == task.lastLeft) {
        // We didn't retract on the left. We can use this code path to fill in
        // the setd for roots and right retractions.
        toRetract.lastLeftPositions = opposingResults;
        // The old results have now been updated to include the new results
        toRetract.lastRightPositions = oldResults;
    } else {
        // We did retract on the left, flip these around.
        toRetract.lastLeftPositions = oldResults;
        toRetract.lastRightPositions = opposingResults;
    }
    
    if(task.leftContext > 1) {
        // Retract on the left and enqueue that child task.
        Log::debug() << "Queueing retraction on the left" << std::endl;
        taskQueue.push(toRetract.retract(view.getIndex(), false));
    }
    
    if(task.isRightEdge && task.rightContext > 1) {
        Log::debug() << "Queueing retraction on the right" << std::endl;
        // Also retract on the right, if the task we just did never retracted on
        // the left.
        taskQueue.push(toRetract.retract(view.getIndex(), true));
    }
    
    Log::debug() << "Found no shared results" << std::endl;
    
    // Say we found nothing in common but we don't want to kill the mapping.
    return {true, {}};
        
}

std::set<TextPosition> ZipMappingScheme::exploreRetractions(
    const FMDPosition& left, size_t patternLengthLeft, const FMDPosition& right,
    size_t patternLengthRight, const std::string& query,
    size_t queryBase) const {

    // List the unique TextPosition we found, if we found one. Holds more than
    // one TextPosition (though not necessarily all of them) if we're ambiguous,
    // and none if we have no results. Holds TextPositions for right contexts.
    std::set<TextPosition> toReturn;

    // We have this DP space we have to traverse defined in 2d by how much
    // context we have left on either side. Call this l and r. If we find we
    // have any overlap at (l, r), we don't have to explore (l1<=l, r1<=r)
    // because we know it will not have any new unique matches.
    
    // I can do my DP by always considering retracting on the left from a state,
    // and only considering retracting on the right on the very edge of the
    // space (i.e. if the left is still full-length).
    
    // Make a queue of DP tasks, and put in the root task.
    std::queue<DPTask> taskQueue({DPTask(left, patternLengthLeft, right,
        patternLengthRight)});
        
    // Make a map of the minimum r we will ever have to explore for any l this
    // size or smaller. We can use lower_bound for the lookups so we only have
    // to insert l,r pairs where we have overlap.
    std::map<size_t, size_t> minRightContext;
        
    // A function to see if we need to text a certain combination of left and
    // right context lengths, or if it's a more general context of one we
    // already had overlapping results for
    auto needToTest = [&](const DPTask& task) {
        if(task.leftContext == 0 || task.rightContext == 0) {
            // Way too short
            return false;
        }
        
        // Find the min right context for things with left
        // contexts greater than or equal to this one.
        auto minRightIterator = minRightContext.lower_bound(task.leftContext);
        
        if(minRightIterator != minRightContext.end() &&
            task.rightContext <= (*minRightIterator).second) {
            // We've retracted far enough that we don't need to process this
            // retraction, since it's a less general context than one we already
            // found results for (on the left it's a lower bound, and we just
            // saw it's a lower bound on the right).
            
            Log::debug() << "Skipping " << task.leftContext << ", " <<
                task.rightContext << " as it is covered by " <<
                (*minRightIterator).first << ", " <<
                (*minRightIterator).second << std::endl;
            
            return false;
        }
        
        // We do need to process this retraction.
        return true;
    };
        
    while(taskQueue.size() > 0) {
    
        if(!needToTest(taskQueue.front())) {
            // Skip queued tasks that have becomne redundant.
            taskQueue.pop();
            continue;
        }
    
        // Run the first task in the queue, possibly adding more, and getting
        // some results.
        auto flagAndSet = exploreRetraction(taskQueue.front(), taskQueue, 
            query, queryBase);
            
        // Drop taht task we just did.
        taskQueue.pop();
        
        if(!flagAndSet.first) {
            // We encountered something too hard to do, and have to abort
            // mapping.
            Log::debug() << "Aborting mapping base " << queryBase <<
                " because it is too hard." << std::endl;
            // Return an empty set.
            return std::set<TextPosition>();
        }
        
        for(auto position : flagAndSet.second) {
            // Put in all the places we found matchings to to.
            toReturn.insert(position);
        }
        
        if(toReturn.size() > 1) {
            // We're already ambiguous. Short circuit.
            Log::debug() << "Already ambiguous, not retracting any more" <<
                std::endl;
            return toReturn;
        }
        
    }
    
    Log::debug() << "Found " << toReturn.size() << " locations" << std::endl;
    
    // If we get here, we've finished all the DP jobs. Return all the
    // TextPositions we found where they were the unique overlap at some place
    // we explored in the DP table.
    return toReturn;

} 

