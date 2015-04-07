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
    
    // Create vectors of the minimum unique left and right context lengths, or
    // (size_t) -1 if no such lenght exists.
    std::vector<size_t> minUniqueLeftContexts;
    std::vector<size_t> minUniqueRightContexts;
    
    for(const auto& rangeAndLength : leftContexts) {
        // Find the min unique lengths for the left contexts.
        
        FMDPosition position = rangeAndLength.first;
        size_t length = rangeAndLength.second;
        
        // What's the min unique length for this base? We haven't found any yet.
        size_t minUniqueLength = (size_t) -1;
        
        while(view.isUnique(position)) {
            // Retract to a point where we may not be unique
            length = view.getIndex().retractRightOnly(position);
            // We know we are unique at 1 more character than that
            minUniqueLength = length + 1;
        }
        
        minUniqueLeftContexts.push_back(minUniqueLength);
        
        Log::trace() << "Min left context " <<
            minUniqueLeftContexts.size() - 1 << " is " << minUniqueLength <<
            " vs " << length << " selecting " <<
            view.getTextPositions(position).size() << std::endl;
        
    }
    
    for(const auto& rangeAndLength : rightContexts) {
        // And again for the right contexts. TODO: make a function to not repeat
        // all this code.
        
        FMDPosition position = rangeAndLength.first;
        size_t length = rangeAndLength.second;
        
        // What's the min unique length for this base? We haven't found any yet.
        size_t minUniqueLength = (size_t) -1;
        
        while(view.isUnique(position)) {
            // Retract to a point where we may not be unique
            length = view.getIndex().retractRightOnly(position);
            // We know we are unique at 1 more character than that
            minUniqueLength = length + 1;
        }
        
        minUniqueRightContexts.push_back(minUniqueLength);
    }
    
    // We're going to store up all the potential mappings, and then filter them
    // down. This stores true and the mapped-to TextPosition if a base would map
    // before filtering, and false and an undefined TextPosition otherwise.
    std::vector<Mapping> mappings;
    
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
            mappings.push_back(Mapping(*(intersection.begin())));
        } else if(intersection.size() > 1) {
            // Too many results
            Log::debug() << "Index " << i << " is ambiguous." << std::endl;
            mappings.push_back(Mapping());
        } else {
            // Too few results
            Log::debug() << "Index " << i << " not in reference." << std::endl;
            mappings.push_back(Mapping());
        }
        
    }
    
    for(size_t i = 0; i < mappings.size(); i++) {
        // Now we have to filter the mappings
        
        // Grab the mapping for this query index.
        const Mapping& mapping = mappings[i];
        
        if(!mapping.isMapped()) {
            // Skip unmapped query bases
            continue;
        }
        
        // We're going to find all the one-sided minimal unique matches
        // contained in our maximal match, and put them in this, and then do
        // activity selection to count how many are non-overlapping.
        std::vector<std::pair<size_t, size_t>> uniqueContexts;
        
        // What range of bases are within our left and right MUMs? TODO:
        // substitute range actually used to map on instead?
        size_t rangeStart = i - leftContexts[i].second + 1;
        size_t rangeEnd = i + rightContexts[i].second - 1; // Inclusive
        
        // Scan left and right to find the one-sided MUSes we overlap.
        for(size_t j = rangeStart; j <= rangeEnd; j++) {
            // For every base in the range
            
            if(minUniqueLeftContexts[j] != (size_t) -1) {
                // If it has a left context
                if(j - minUniqueLeftContexts[j] + 1 >= rangeStart) {
                    // And that left context is contained in our range
                
                    // Save this one-sided MUS
                    uniqueContexts.emplace_back(
                        j - minUniqueLeftContexts[j] + 1, j);
                
                }
            }
            
            if(minUniqueRightContexts[j] != (size_t) -1) {
                // If it has a right context
                if(j + minUniqueRightContexts[j] - 1 <= rangeEnd) {
                    // And that right context is contained in our range
                    
                    // Save this one-sided MUS
                    uniqueContexts.emplace_back(j,
                        j + minUniqueRightContexts[j] - 1);
                }
            }
        }
        
        Log::debug() << "Query position " << i << " contains " << 
            uniqueContexts.size() << " unique strings in its range from " <<
            rangeStart << " to " << rangeEnd << std::endl;
            
        // Do activity selection
        size_t nonOverlapping = selectActivities(uniqueContexts);
        
        Log::debug() << nonOverlapping << "/" << uniqueContexts.size() <<
            " are nonoverlapping" << std::endl;
        
        if(nonOverlapping < minUniqueStrings) {
            Log::debug() <<
                "Dropping mapping due to having too few unique strings." <<
                std::endl;
        } else {
            // Report all the mappings that pass.
            callback(i, mappings[i].getLocation());
        }
        
    }
    
}

bool ZipMappingScheme::canExtendThrough(FMDPosition context,
    const std::string& opposingQuery) const {
    
    Log::debug() << "Trying to extend " << context << " through " << 
        opposingQuery.size() << " opposing context" << std::endl;
        
    // We're going to retract it until it's no longer unique, then go
    // back and retract it one less.
    FMDPosition noLongerUnique = context;
    size_t nonUniqueLength = view.getIndex().retractRightOnly(noLongerUnique);
    
    while(view.isUnique(noLongerUnique)) {
        // It may take multiple retracts because of masks and stuff.
        nonUniqueLength = view.getIndex().retractRightOnly(noLongerUnique);
    }
    
    // Go back and retract one less (to a length one longer).
    FMDPosition barelyUnique = context;
    view.getIndex().retractRightOnly(barelyUnique, nonUniqueLength + 1);

    for(size_t i = opposingQuery.size() - 2; i != (size_t) -1 &&
        !view.isEmpty(barelyUnique); i--) {
        // Now go extend through with the opposing string, from right to left.
        // We skip the rightmost character (the one we are actually in the
        // process of mapping) because that's already searched, and we keep
        // going until we run out of results or we make it all the way through.
        
        Log::trace() << "Extending with " << opposingQuery[i] << std::endl;
        
        view.getIndex().extendLeftOnly(barelyUnique, opposingQuery[i]);
    }
    
    if(view.isEmpty(barelyUnique)) {
        Log::debug() << "Extension failed" << std::endl;
    }
    
    // Let the caller know whether we got all the way through (in which case the
    // LR context exists) or not.
    return !view.isEmpty(barelyUnique);
}

size_t ZipMappingScheme::selectActivities(
    std::vector<std::pair<size_t, size_t>> ranges) const {
    
    // First we have to sort the ranges by end, ascending.
    std::sort(ranges.begin(), ranges.end(), [](std::pair<size_t, size_t> a,
        std::pair<size_t, size_t> b) -> bool {
        
        // We'll give true if the first range ends before the second.
        // We also return true if there's a tie but a starts first.
        return (a.second < b.second) || (a.second == b.second &&
            a.first < b.first);
        
    });
    
    // How many non-overlapping things have we found?
    size_t found = 0;
    // What's the minimum start time we can accept?
    size_t nextFree = 0;
    
    for(const auto& range : ranges) {
        Log::trace() << "Have " << range.first << " - " << range.second <<
            std::endl;
        // For each range (starting with those that end soonest)...
        if(range.first >= nextFree) {
            // Take it if it starts late enough
            
            Log::trace() << "Taking " << range.first << " - " << range.second <<
                std::endl;
            
            nextFree = range.second + 1;
            found++;
        }
    }

    return found;
}

std::pair<bool, std::set<TextPosition>> ZipMappingScheme::exploreRetraction(
        const DPTable::DPTask& task, DPTable& table, const std::string& query,
        size_t queryBase) const {
    
    // Get the left and right SideRetractionEntries. Also computes any
    // uncomputed retractions and enumerates any sufficiently small sets.
    const DPTable::SideRetractionEntry& left = table.getRetraction(false,
        task.leftIndex, view, maxRangeCount);
    const DPTable::SideRetractionEntry& right = table.getRetraction(true,
        task.rightIndex, view, maxRangeCount);
    
    Log::debug() << "Exploring " << left.position << " (" <<
        left.contextLength << ") left, " << right.position << " (" <<
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
          
    
    if(view.isUnique(left.position) && right.contextLength < maxExtendThrough) {
        
        // We will try extending the left through the right.
        triedLeftExtend = true;
        
        if(canExtendThrough(left.position, reverseComplement(query.substr(
            queryBase, right.contextLength)))) {
            // We tried to extend the left context through the (reverse
            // complemented) right context and succeeded.
            
            Log::debug() << "Successfully extended left through right" <<
                std::endl;
            
            // We can now cheat and return the one-element set of whatever the
            // left context selects, flipped.
            TextPosition matchedTo = view.getTextPosition(left.position);
            matchedTo.flip(view.getIndex().getContigLength(
                matchedTo.getContigNumber()));
            return {true, {matchedTo}};
        } else {
            Log::debug() << "Failed to extend left though right" << std::endl;
        }
        
        // If we don't succeed, we still need to check for set overlap.
    }
    if(view.isUnique(right.position) && left.contextLength < maxExtendThrough) {
        // TODO: Don't try this extension if we failed the other direction.
        
        // We will try extending the right through the left.
        triedRightExtend = true;
        
        if(canExtendThrough(right.position, query.substr(
            queryBase - (left.contextLength - 1), left.contextLength))) {
            // We tried to extend the right context through the left context and
            // succeeded
            
            Log::debug() << "Successfully extended right through left" <<
                std::endl;
            
            // We can now cheat and return the one-element set of whatever the
            // right context selects.
            TextPosition matchedTo = view.getTextPosition(right.position);
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
    DPTable::DPTask toRetract = task;
    
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

std::set<TextPosition> ZipMappingScheme::exploreRetractions(
    const FMDPosition& left, size_t patternLengthLeft, const FMDPosition& right,
    size_t patternLengthRight, const std::string& query,
    size_t queryBase) const {

    // List the unique TextPosition we found, if we found one. Holds more than
    // one TextPosition (though not necessarily all of them) if we're ambiguous,
    // and none if we have no results. Holds TextPositions for right contexts.
    std::set<TextPosition> toReturn;

    // Make a DP table for this base. TODO: make DPTable remember the extra
    // parameters.
    DPTable table(left, patternLengthLeft, right, patternLengthRight,
        view, maxRangeCount);
        
    // I can do my DP by always considering retracting on the left from a state,
    // and only considering retracting on the right on the very edge of the
    // space (i.e. if the left is still full-length).
        
    // Make a map of the minimum right context (r) we will ever have to explore
    // for any left context (l) this size or smaller. We can use lower_bound for
    // the lookups so we only have to insert l,r pairs where we have overlap.
    std::map<size_t, size_t> minRightContext;
        
    // A function to see if we need to text a certain combination of left and
    // right context lengths, or if it's a more general context of one we
    // already had overlapping results for
    auto needToTest = [&](const DPTable::DPTask& task) {
        // Work out the context lengths you would have on the left and right for
        // this task (and probably actually compute the retracted sets too).
        const DPTable::SideRetractionEntry& left = table.getRetraction(false,
            task.leftIndex, view, maxRangeCount);
        const DPTable::SideRetractionEntry& right = table.getRetraction(true,
            task.rightIndex, view, maxRangeCount);
        
        if(left.contextLength == 0 || right.contextLength == 0) {
            // Way too short
            return false;
        }
        
        // Find the min right context for things with left contexts greater than
        // or equal to this one.
        auto minRightIterator = minRightContext.lower_bound(left.contextLength);
        
        if(minRightIterator != minRightContext.end() &&
            right.contextLength <= (*minRightIterator).second) {
            // We've retracted far enough that we don't need to process this
            // retraction, since it's a less general context than one we already
            // found results for (on the left it's a lower bound, and we just
            // saw it's a lower bound on the right).
            
            Log::debug() << "Skipping " << left.contextLength << ", " <<
                right.contextLength << " as it is covered by " <<
                (*minRightIterator).first << ", " <<
                (*minRightIterator).second << std::endl;
            
            return false;
        }
        
        // We do need to process this retraction.
        return true;
    };
        
    while(table.taskQueue.size() > 0) {
    
        if(!needToTest(table.taskQueue.front())) {
            // Skip queued tasks that have becomne redundant.
            table.taskQueue.pop();
            continue;
        }
    
        // Run the first task in the queue, possibly adding more, and getting
        // some results.
        auto flagAndSet = exploreRetraction(table.taskQueue.front(), table, 
            query, queryBase);
            
        // Drop taht task we just did.
        table.taskQueue.pop();
        
        if(!flagAndSet.first) {
            // We encountered something too hard to do.
            
            if(giveUpIfHard) {
                // We have to abort mapping.
                Log::debug() << "Aborting mapping base " << queryBase <<
                    " because it is too hard." << std::endl;
                // Return an empty set.
                return std::set<TextPosition>();
            } else {
                Log::debug() << "Ignoring hard task" << std::endl;
            }
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
        
        if(!useRetraction) {
            // Only explore the very first task, which is the root and not
            // actually retracting on either side.
            break;
        }
        
    }
    
    Log::debug() << "Found " << toReturn.size() << " locations" << std::endl;
    
    // If we get here, we've finished all the DP jobs. Return all the
    // TextPositions we found where they were the unique overlap at some place
    // we explored in the DP table.
    return toReturn;

} 

