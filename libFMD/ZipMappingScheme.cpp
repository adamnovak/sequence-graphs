#include "ZipMappingScheme.hpp"
#include "Log.hpp"
#include "util.hpp"

#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

#include <boost/range/adaptor/reversed.hpp>

template<>
bool ZipMappingScheme<FMDPosition>::canExtendThrough(FMDPosition context,
    const std::string& opposingQuery) const {
    
    Log::debug() << "Trying to extend " << context << " through " << 
        opposingQuery.size() << " opposing context" << std::endl;
    stats.add("extendThroughAttampts", 1);
        
    // We're going to retract it until it's no longer unique, then go
    // back and retract it one less.
    FMDPosition noLongerUnique = context;
    size_t nonUniqueLength = view.getIndex().retractRightOnly(noLongerUnique);
    
    while(noLongerUnique.isUnique(view)) {
        // It may take multiple retracts because of masks and stuff.
        nonUniqueLength = view.getIndex().retractRightOnly(noLongerUnique);
    }
    
    // Go back and retract one less (to a length one longer).
    FMDPosition barelyUnique = context;
    view.getIndex().retractRightOnly(barelyUnique, nonUniqueLength + 1);

    for(size_t i = opposingQuery.size() - 2; i != (size_t) -1 &&
        !barelyUnique.isEmpty(view); i--) {
        // Now go extend through with the opposing string, from right to left.
        // We skip the rightmost character (the one we are actually in the
        // process of mapping) because that's already searched, and we keep
        // going until we run out of results or we make it all the way through.
        
        Log::trace() << "Extending with " << opposingQuery[i] << std::endl;
        
        view.getIndex().extendLeftOnly(barelyUnique, opposingQuery[i]);
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
size_t ZipMappingScheme<FMDPosition>::selectActivities(
    std::vector<std::pair<size_t, size_t>> ranges) const {
    
    stats.add("activitySelectionRuns", 1);
    
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
            
            stats.add("activitiesSelected", 1);
        }
    }

    return found;
}

template<>
std::vector<std::pair<FMDPosition, size_t>>
    ZipMappingScheme<FMDPosition>::findRightContexts(const std::string& query,
    bool reverse) const {
 
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
        
        while(extended.isEmpty(view)) {
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
std::vector<size_t> ZipMappingScheme<FMDPosition>::createUniqueContextIndex(
    const std::vector<std::pair<FMDPosition, size_t>>& leftContexts,
    const std::vector<std::pair<FMDPosition, size_t>>& rightContexts) const {
    
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
        
        while(position.isUnique(view)) {
            // Retract to a point where we may not be unique
            length = view.getIndex().retractRightOnly(position);
            // We know we are unique at 1 more character than that
            minUniqueLength = length + 1;
        }
        
        if(minUniqueLength == 0) {
            throw std::runtime_error("Too much retracting!");
        }
        
        minUniqueLeftContexts.push_back(minUniqueLength);
        
        Log::trace() << "Min left context " <<
            minUniqueLeftContexts.size() - 1 << " is " << minUniqueLength <<
            " vs " << length << " selecting " <<
            position.getTextPositions(view).size() << std::endl;
        
    }
    
    for(const auto& rangeAndLength : rightContexts) {
        // And again for the right contexts. TODO: make a function to not repeat
        // all this code.
        
        FMDPosition position = rangeAndLength.first;
        size_t length = rangeAndLength.second;
        
        // What's the min unique length for this base? We haven't found any yet.
        size_t minUniqueLength = (size_t) -1;
        
        while(position.isUnique(view)) {
            // Retract to a point where we may not be unique
            length = view.getIndex().retractRightOnly(position);
            // We know we are unique at 1 more character than that
            minUniqueLength = length + 1;
        }
        
        minUniqueRightContexts.push_back(minUniqueLength);
    }
    
    // Now we're going to make a sorted list of (start, end) pairs for these
    // one-sided MUSes
    std::vector<std::pair<size_t, size_t>> uniqueContexts;
    
    for(size_t i = 0; i < minUniqueLeftContexts.size(); i++) {
        // For every base
        
        if(minUniqueLeftContexts[i] != (size_t) -1) {
            // If it has a left context,  Save this one-sided MUS
            uniqueContexts.emplace_back(
                i - minUniqueLeftContexts[i] + 1, i);
      
            if(minUniqueLeftContexts[i] == 0) {
                throw std::runtime_error("Left context runs backwards!");
            }
        }
        
        if(minUniqueRightContexts[i] != (size_t) -1) {
            // If it has a right context, save this one-sided MUS
            uniqueContexts.emplace_back(i,
                i + minUniqueRightContexts[i] - 1);
                
            if(minUniqueRightContexts[i] == 0) {
                throw std::runtime_error("Right context runs backwards!");
            }
        }
    }
    
    // The "first" range is the one that ends earliest.
    auto order = ([](std::pair<size_t, size_t> a, 
        std::pair<size_t, size_t> b) -> bool {
        
        // We'll give true if the first range ends before the second.
        // We also return true if there's a tie but a starts first.
        return (a.second < b.second) || (a.second == b.second &&
            a.first < b.first);
        
    });
    
    // Sort the unique contexts by start, descending, so we can pop stuff off
    // the front when we reach its start position.
    std::sort(uniqueContexts.begin(), uniqueContexts.end(),
        [](std::pair<size_t, size_t> a, std::pair<size_t, size_t> b) -> bool {
            // Define a total ordering
            return (a.first > b.first) || (a.first == b.first &&
                a.second > b.second);
    });
    
    // Go through the unique contexts from right to left, and save for each base
    // the end position of the earliest-ending unique context that starts at or
    // after it.
    std::vector<size_t> endOfBestActivity(leftContexts.size());
    
    // We're going to iterate through the unique contexts in order of decreasing
    // start index.
    auto nextActivity = uniqueContexts.begin();
    
    // We'll keep an iterator to the one that starts here or to the right and
    // ends soonest, but for now that's nothing.
    auto bestActivity = uniqueContexts.end();
    
    // We go from right to left. We have a current best (i.e. earliest-ending)
    // interval. When we reach the start point of an interval, we see if that
    // interval ends before our current best one. If so we replace the best
    // interval. If not we keep it and ignore the new interval.
    
    for(size_t i = leftContexts.size() - 1; i != (size_t) -1; i--) {
        // For each base from right to left
        
        while(nextActivity != uniqueContexts.end() && 
            (*nextActivity).first >= i) {
            // For every activity starting at or after here
            
            Log::trace() << "Considering activity " << (*nextActivity).first <<
                "-" << (*nextActivity).second << " which starts at or after " <<
                i << std::endl;
            
            // This activity starts here or later
            if(bestActivity == uniqueContexts.end() ||
                (*nextActivity).second < (*bestActivity).second) {
                
                // This is a sooner-ending activity than any we have so far.
                // Take it.
                bestActivity = nextActivity;
                
                Log::trace() << "It is the soonest ending" << std::endl;
                
            }
            // We accepted or rejected this activity, so go on to the next
            // lefter one to consider.
            ++nextActivity;
        }
        
        // Now bestActivity points to the soonest-ending activity that starts
        // here or later.
        if(bestActivity == uniqueContexts.end()) {
            // There is no such activity. Use (size_t) -1 as a no-such-value
            // value.
            endOfBestActivity[i] = (size_t) -1;
            
            Log::trace() << "No soonest ending activity starting at or after " <<
                i << std::endl;
        } else {
            // We found such an activity, so record its endpoint. We can go one
            // right from there and look up the end of the non-overlapping
            // activity that we would chain with.
            endOfBestActivity[i] = (*bestActivity).second;
            
            Log::trace() << "Soonest ending activity starting at or after " <<
                i << " ends at " << endOfBestActivity[i] << std::endl;
        }
    }
    
    // Giev back the index vector we made, which simplifies activity selection
    // problems immensely.
    return endOfBestActivity;
}

template<>
size_t ZipMappingScheme<FMDPosition>::selectActivitiesIndexed(size_t start,
    size_t end, const std::vector<size_t>& index, size_t threshold) const {
    
    // What base are we at?
    size_t position = start;
    // How many non-overlapping things have we found?
    size_t found = 0;

    while(found < threshold && position < index.size()) {
        // Look for activities until we find enough. A threshold of (size_t) -1
        // wil have us look forever. Also make sure we don't sneak out of the
        // index with that +1 later.
    
        // Jump to the end of the next activity we want to do.
        position = index[position];
        
        if(position > end) {
            // This activity runs too long and we can't do it. This also handles
            // (size_t) -1, the no-more-activities value, which is larger than
            // any other size_t.
            break;
        }
        
        // Otherwise we found an activity that we can do.
        found++;
        
        // Look for activities that start after this one ends.
        position++;
    }
    
    // Say how many activities were found.
    return found;
}

// Forward-declare specialization for mutual recursion of specializations.
template<>
Mapping ZipMappingScheme<FMDPosition>::exploreRetractions(
    const FMDPosition& left, size_t patternLengthLeft, const FMDPosition& right,
    size_t patternLengthRight, const std::string& query,
    size_t queryBase) const;

template<>
std::pair<bool, std::set<TextPosition>>
    ZipMappingScheme<FMDPosition>::exploreRetraction(
    const typename DPTable::DPTask& task, DPTable& table, const std::string& query,
    size_t queryBase) const {
    
    // Get the left and right SideRetractionEntries. Also computes any
    // uncomputed retractions and enumerates any sufficiently small sets.
    const DPTable::SideRetractionEntry& left = table.getRetraction(false,
        task.leftIndex, view, maxRangeCount);
    const DPTable::SideRetractionEntry& right = table.getRetraction(true,
        task.rightIndex, view, maxRangeCount);
    
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


template<>
Mapping ZipMappingScheme<FMDPosition>::exploreRetractions(
    const FMDPosition& left, size_t patternLengthLeft, const FMDPosition& right,
    size_t patternLengthRight, const std::string& query,
    size_t queryBase) const {

    stats.add("basesAttempted", 1);

    // List the unique TextPosition we found, if we found one. Holds more than
    // one TextPosition (though not necessarily all of them) if we're ambiguous,
    // and none if we have no results. Holds TextPositions for right contexts.
    std::set<TextPosition> found;
    
    // What are the max left and right contexts that were used to map the base
    // anywhere? Note that they may not be part of the same string that maps the
    // base there, but they would have to be if we accounted for crossover.
    size_t maxLeftContext = 0;
    size_t maxRightContext = 0;

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
            
            // Note that there was a covered retraction.
            stats.add("retractionCovered", 1);
            
            return false;
        }
        
        // We do need to process this retraction.
        return true;
    };
        
    while(table.taskQueue.size() > 0) {
    
        // Grab the task in the table
        const DPTable::DPTask& task = table.taskQueue.front();
    
        if(!needToTest(task)) {
            // Skip queued tasks that have becomne redundant.
            table.taskQueue.pop();
            continue;
        }
    
        // Grab the left and right retraction entries.
        // TODO: Stop repeating this code somehow.
        const DPTable::SideRetractionEntry& left = table.getRetraction(false,
            task.leftIndex, view, maxRangeCount);
        const DPTable::SideRetractionEntry& right = table.getRetraction(true,
            task.rightIndex, view, maxRangeCount);
            
        // Which we needed to pull out the context lengths on which we may be
        // mapping.
        size_t leftContext = left.contextLength;
        size_t rightContext = right.contextLength;
        
        // Which we in turn need to calculate the max left and right context
        // lengths observed for the base.
        maxLeftContext = std::max(maxLeftContext, leftContext);
        maxRightContext = std::max(maxRightContext, rightContext);
    
        // Run the first task in the queue, possibly adding more, and getting
        // some results.
        auto flagAndSet = exploreRetraction(task, table, 
            query, queryBase);
            
        // Drop that task we just did. We can't use task anymore now, or left or
        // right.
        table.taskQueue.pop();
        
        if(!flagAndSet.first) {
            // We encountered something too hard to do.
            stats.add("tooHardRetraction", 1);
            
            if(giveUpIfHard) {
                // We have to abort mapping.
                Log::debug() << "Aborting mapping base " << queryBase <<
                    " because it is too hard." << std::endl;
                // Return an empty mapping.
                return Mapping();
            } else {
                Log::debug() << "Ignoring hard task" << std::endl;
            }
        }
        
        for(auto position : flagAndSet.second) {
            // Put in all the places we found matchings to to.
            found.insert(position);
        }
        
        if(found.size() > 1) {
            // We're already ambiguous. Short circuit.
            Log::debug() << "Already ambiguous, not retracting any more" <<
                std::endl;
            stats.add("ambiguous", 1);
            return Mapping();
        }
        
        if(!useRetraction) {
            // Only explore the very first task, which is the root and not
            // actually retracting on either side.
            break;
        }
        
    }
    
    Log::debug() << "Found " << found.size() << " locations" << std::endl;
    Log::debug() << "Used" << maxLeftContext << ", " << maxRightContext <<
        " context" << std::endl;
        
    if(found.size() == 1) {
        // We mapped to one place, on these contexts.
        stats.add("unambiguous", 1);
        return Mapping(*(found.begin()), maxLeftContext, maxRightContext);
    } else {
        // We mapped to nowhere, because we're ambiguous. TODO: log
        // ambiguousness.
        stats.add("ambiguous", 1);
        return Mapping();
    }
} 

template<>
void ZipMappingScheme<FMDPosition>::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // Get the right contexts
    Log::info() << "Looking for right contexts..." << std::endl << std::flush;
    auto rightContexts = findRightContexts(query, false);
    
    // And the left contexts (which is the reverse of the contexts for the
    // reverse complement, which we can produce in backwards order already).
    Log::info() << "Flipping query..." << std::endl << std::flush;
    std::string complement = reverseComplement(query);
    
    Log::info() << "Looking for left contexts..." << std::endl << std::flush;
    auto leftContexts = findRightContexts(complement, true);
    
    // This is going to hold, for each query base, the endpoint of the soonest-
    // ending minimally unique context that starts at or after that base.
    Log::debug() << "Creating unique context index" << std::endl << std::flush;
    auto uniqueContextIndex = createUniqueContextIndex(leftContexts,
        rightContexts);
    
    
    Log::info() << "Exploring retractions..." << std::endl;
    
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
        Mapping mapping = exploreRetractions(
            leftContexts[i].first, leftContexts[i].second,
            rightContexts[i].first, rightContexts[i].second, query, i);
        
        // TODO: If we can't find anything, try retracting a few bases on one or
        // both sides until we get a shared result.
        
        if(mapping.isMapped()) {
            // We map!
            Log::debug() << "Index " << i << " maps on " << 
                mapping.getLeftMaxContext() << " left and " << 
                mapping.getRightMaxContext() << " right context" << std::endl;
        } else {
            // Too few results until we retracted back to too many
            Log::debug() << "Index " << i << " is not mapped." << std::endl;
       }
       // Save the mapping
       mappings.push_back(mapping);
        
    }
    
    Log::info() << "Applying filter..." << std::endl << std::flush;
    
    // We're going to filter everything first and then do the callbacks.
    std::vector<Mapping> filtered;
    
    for(size_t i = 0; i < mappings.size(); i++) {
        // Now we have to filter the mappings
        
        // Grab the mapping for this query index.
        const Mapping& mapping = mappings[i];
        
        if(!mapping.isMapped()) {
            // Skip unmapped query bases
            filtered.push_back(mapping);
            continue;
        }
        
        // What range of bases are within our left and right MUMs?
        size_t rangeStart = i - mapping.getLeftMaxContext() + 1;
        size_t rangeEnd = i + mapping.getRightMaxContext() - 1; // Inclusive
        
        // Do activity selection
        size_t nonOverlapping = selectActivitiesIndexed(rangeStart, rangeEnd,
            uniqueContextIndex, minUniqueStrings);
            
        if(nonOverlapping < minUniqueStrings) {
            Log::debug() <<
                "Dropping mapping due to having too few unique strings." <<
                std::endl;
            // Report a non-mapping Mapping
            stats.add("filterFail", 1);
            filtered.push_back(Mapping());
        } else {
            // Report all the mappings that pass.
            stats.add("filterPass", 1);
            filtered.push_back(mapping);
        }
        
    }
    
    if(credit.enabled) {
        // If credit isn't disabled, run it. The cannonical check for enabled-
        // ness is in the CreditStrategy itself, but no point running it if it's
        // going to do nothing.
        Log::info() << "Applying credit..." << std::endl;
        credit.applyCredit(query, filtered);
    }
    
    Log::info() << "Sending callbacks..." << std::endl << std::flush;
    for(size_t i = 0; i < filtered.size(); i++) {
        if(filtered[i].isMapped()) {
            // Send a callback for everything that passed the filter.
            callback(i, filtered[i].getLocation());
        }
    }
    
}

