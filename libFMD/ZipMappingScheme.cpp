#include "ZipMappingScheme.hpp"
#include "Log.hpp"
#include "util.hpp"

#include <vector>
#include <map>
#include <queue>
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
            rightContexts[i].first, rightContexts[i].second);
        
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

std::set<TextPosition> ZipMappingScheme::exploreRetractions(
    const FMDPosition& left, size_t patternLengthLeft,
    const FMDPosition& right, size_t patternLengthRight) const {

    // List the unique TextPosition we found, if we found one. Hosld more than
    // one TextPosition (though not necessarily all of them) if we're ambiguous,
    // and none if we have no results. Holds TextPositions for right contexts.
    std::set<TextPosition> toReturn;

    // We have this DP space we have to traverse defined in 2d by how much we
    // retract on either side. Call this l and r. If we find we have any overlap
    // at (l, r), we don't have to explore (l1>=l, r1>=r) because we know it
    // will not have any new unique matches.
    
    // I can do my DP by always considering retracting on the left from a state,
    // and only considering retracting on the right on the very edge of the
    // space (i.e. if the left is still full-length).
    
    if(patternLengthLeft + patternLengthRight - 1 < minContextLength) {
        // There isn't enough context for this base even with nothing retracted.
        return toReturn;
    }
    
    // First check the maximum extensions we have been given.
    
    if(left.getEndOffset() < 5) {
        for(size_t i = 0; i <= left.getEndOffset(); i++) {
            // Announce each actual BWT index we have selected on the left.            
            Log::debug() << "Left has selected BWT " <<
                left.getForwardStart() + i << " = " <<
                view.getIndex().locate(left.getForwardStart() + i) << std::endl;
        }
    }
    
    // Find the reverse-orientation positions selected by the right context.
    std::set<TextPosition> rightPositions = view.getTextPositions(right);
        
    // And the forward-orientation positions selected by the left contexts.
    std::set<TextPosition> leftPositions = view.getTextPositions(left);
    
    Log::debug() << leftPositions.size() << " positions on the left, " <<
        rightPositions.size() << " on the right" << std::endl;
        
    if(leftPositions.size() < 5) {
        Log::debug() << "Left:" << std::endl;
        for(auto position : leftPositions) {
            position.flip(view.getIndex().getContigLength(
                position.getContigNumber()));
            Log::debug() << "\t" << position << std::endl;
        }
    }
    
    if(rightPositions.size() < 5) {
        Log::debug() << "Right:" << std::endl;
        for(auto position : rightPositions) {
            Log::debug() << "\t" << position << std::endl;
        }
    }
    
    for(auto leftPosition : leftPositions) {
        // Flip the position from the reverse complement around so it's on
        // the right strand. TODO: this pattern should be simpler.
        leftPosition.flip(view.getIndex().getContigLength(
            leftPosition.getContigNumber()));
    
        if(rightPositions.count(leftPosition)) {
            // We found an intersection before even retracting.
            toReturn.insert(leftPosition);
            Log::debug() << "TextPosition " << leftPosition <<
                " is shared at top level." << std::endl;
               
            if(toReturn.size() > 1) { 
                // Fail mapping fast if we find more than one at this top level.
                return toReturn;
            }
        }        
    }
    
    if(toReturn.size() > 0) {
        // We found results without even doing any DP. It may be unique or not,
        // but either way we can just return the whole intersection we found.
        return toReturn;
    }
    
    if(!useRetraction) {
        // We let the user disable our retraction operations.
        Log::debug() << "Skipping retraction according to configuration." <<
            std::endl;
        return toReturn;
    }
    
    // If we get here, we had no intersection with the longest contexts on each
    // side, so we need to do the DP and try various retractions.
    
    // Define a type for an FMDPosition retracted a certain distance from the
    // original, with a set of selected TextPositions.
    using Retraction = std::tuple<FMDPosition, size_t, std::set<TextPosition>>;
    // Define a dynamic programming state as a pair of retractions: left and
    // right.
    using DPState = std::pair<Retraction, Retraction>;
    
    // Make a queue of the retractions to process.
    std::queue<DPState> todo;
    
    // Start by handling no retraction on either side.
    todo.push({std::make_tuple(left, 0, leftPositions), std::make_tuple(right,
        0, rightPositions)});
    
    // Make a map of the maximum r we will ever have to explore for any given l.
    // We can use upper_bound for the lookups so we only have to insert l,r
    // pairs where we have overlap.
    std::map<size_t, size_t> maxRightRetract;
    
    // A function to see if we need to text a certain combination of left and
    // right retractions, or if it's a more general context of one we already
    // had overlapping results for (or out of range for the pattern lengths)
    auto needToTest = [&](size_t leftRetraction, size_t rightRetraction) {
        if(leftRetraction >= patternLengthLeft ||
            rightRetraction >= patternLengthRight) {
            // We would retract all of the characters (or more) on one or both
            // sides. Since we need to keep a nonempty search string, we
            // shouldn't expore out to here.
            
            Log::debug() << leftRetraction << ", " <<
                rightRetraction << " is too far out of bounds " <<
                patternLengthLeft << ", " << patternLengthRight << std::endl;
            
            return false;
        }
        
        if(patternLengthLeft - leftRetraction + patternLengthRight - 
            rightRetraction - 1 < minContextLength) {
            
            // Apply our min context length. Subtract 1 because each direction
            // counts the present base.
            Log::debug() << "Context length would be too low" << std::endl;
            
            return false;
        }
        
        // Find the max right retraction distance for things with left
        // retractions equal to this one.
        auto maxRightIterator = maxRightRetract.find(leftRetraction);
        if(maxRightIterator == maxRightRetract.end()) {
            // Find the max right retraction distance for things with left
            // retractions less than this one.
            maxRightIterator = maxRightRetract.upper_bound(leftRetraction);
        }
        
        if(maxRightIterator != maxRightRetract.end() &&
            (*maxRightIterator).second <= rightRetraction) {
            // We've retracted far enough that we don't need to process this
            // retraction, since it's a less general context than one we already
            // found results for.
            
            Log::debug() << "Skipping " << leftRetraction << ", " <<
                rightRetraction << " as it is covered by " <<
                (*maxRightIterator).first << ", " <<
                (*maxRightIterator).second << std::endl;
            
            return false;
        }
        
        // We do need to process this retraction.
        return true;
    };
    
    Log::debug() << "Trying retractions..." << std::endl;
    
    while(todo.size() > 0) {
        // Process things off the queue until we run out.
        auto state = todo.front();
        todo.pop();
        
        // Pull out the result ranges, the retraction lengths, and the sets of
        // selected TextPositions on either side.
        FMDPosition leftResults = std::get<0>(state.first);
        size_t leftRetraction = std::get<1>(state.first);
        const std::set<TextPosition>& leftSet = std::get<2>(state.first);
        
        FMDPosition rightResults = std::get<0>(state.second);
        size_t rightRetraction = std::get<1>(state.second);
        const std::set<TextPosition>& rightSet = std::get<2>(state.first);
        
        // See if the spot down and to the left of us needs doing
        if(needToTest(leftRetraction + 1, rightRetraction)) {
        
            Log::debug() << "Considering retraction " << leftRetraction + 1 <<
                ", " << rightRetraction << std::endl;
        
            // Retract on the left. 
            FMDPosition leftRetracted = leftResults;
            view.getIndex().retractRightOnly(leftRetracted,
                patternLengthLeft - leftRetraction - 1);
                
            // See what new stuff we select on the left.
            std::set<TextPosition> newlySelected = view.getNewTextPositions(
                leftResults, leftRetracted);
        
            // Make a set of overlapping positions
            std::set<TextPosition> overlaps;
        
            for(auto newResult : newlySelected) {
                // Flip each result to look it up in the opposing set.
                newResult.flip(view.getIndex().getContigLength(
                    newResult.getContigNumber()));
                    
                if(rightSet.count(newResult)) {
                    // We have an overlap!
                    overlaps.insert(newResult);
                    
                    if(overlaps.size() > 1) { 
                        // Fail fast if we have more than one.
                        break;
                    }
                }
            }
        
        
            if(overlaps.size() > 0) {
                // We have overlap between the new stuff we select on the left
                // and the old stuff we had selected on the right.
                
                Log::debug() << overlaps.size() <<
                    " overlaps found for left child." << std::endl;
                
                // If so, we don't need to look at its descendants.
                maxRightRetract[leftRetraction + 1] = rightRetraction;
                
                if(overlaps.size() == 1) {
                    // If it's unique, we can spit out a unique match here.
                    toReturn.insert(*(overlaps.begin()));
                    
                    if(toReturn.size() > 1) { 
                        // Fail mapping fast if we have multiple unique matches
                        // for different retractions.
                        return toReturn;
                    }
                }
            } else {
                // We do need to look at the descendants
                
                Log::debug() << "No results, queueing left descendant" <<
                    std::endl;
            
                // Expand the left set with the new stuff. This is apparently
                // hard. See <http://choorucode.com/2010/07/16/c-stl-inserting-
                // vector-into-set/>
                std::set<TextPosition> newLeftSet = leftSet;
                std::copy(newlySelected.begin(), newlySelected.end(),
                    std::inserter(newLeftSet, newLeftSet.end()));
                
                // Queue up a DP job to expand off of there.
                todo.push({std::make_tuple(leftRetracted, leftRetraction + 1,
                    newLeftSet), std::make_tuple(rightResults, rightRetraction,
                    rightSet)});
            }
        } else {
            Log::debug() << "Don't need to test left descendant"  << std::endl;
        }
        
        if(leftRetraction == 0 &&
            needToTest(leftRetraction, rightRetraction + 1)) {
        
            Log::debug() << "Considering retraction " << leftRetraction <<
                ", " << rightRetraction + 1 << std::endl;
        
            // We are on the right edge, and the spot down and to the right of
            // us needs doing
            
            // Retract on the right. 
            FMDPosition rightRetracted = rightResults;
            view.getIndex().retractRightOnly(rightRetracted,
                patternLengthRight - rightRetraction - 1);
                
            // See what new stuff we select on the right.
            std::set<TextPosition> newlySelected = view.getNewTextPositions(
                rightResults, rightRetracted);
        
            // Make a set of overlapping positions
            std::set<TextPosition> overlaps;
        
            for(auto newResult : newlySelected) {
                // Flip each result to look it up in the opposing set.
                newResult.flip(view.getIndex().getContigLength(
                    newResult.getContigNumber()));
                    
                if(leftSet.count(newResult)) {
                    // We have an overlap!
                    overlaps.insert(newResult);
                    
                    if(overlaps.size() > 1) { 
                        // Fail fast if we have more than one.
                        break;
                    }
                    
                }
            }
        
        
            if(overlaps.size() > 0) {
                // We have overlap between the new stuff we select on the right
                // and the old stuff we had selected on the left.
                
                // If so, we don't need to look at its descendants. TODO: we
                // already handle this by not making a new job in the queue.
                // Simplify the needToTest system to work with our top-down
                // organization of job children.
                maxRightRetract[leftRetraction] = rightRetraction + 1;
                
                Log::debug() << overlaps.size() <<
                    " overlaps found for right child." << std::endl;
                
                if(overlaps.size() == 1) {
                    // If it's unique, we can spit out a unique match here. But
                    // we need to flip it around.
                    auto overlap = *(overlaps.begin());
                    overlap.flip(view.getIndex().getContigLength(
                        overlap.getContigNumber()));
                    
                    toReturn.insert(overlap);
                    
                    if(toReturn.size() > 1) { 
                        // Fail mapping fast if we have multiple unique matches
                        // for different retractions.
                        return toReturn;
                    }
                }
            } else {
                // We do need to look at the descendants
            
                Log::debug() << "No results, queueing right descendant" <<
                    std::endl;
            
                // Expand the right set with the new stuff.
                // TODO: Somehow avoid a copy?
                std::set<TextPosition> newRightSet = rightSet;
                std::copy(newlySelected.begin(), newlySelected.end(),
                    std::inserter(newRightSet, newRightSet.end()));
                
                // Queue up a DP job to expand off of there.
                todo.push({std::make_tuple(leftResults, leftRetraction,
                    leftSet), std::make_tuple(rightRetracted,
                    rightRetraction + 1, newRightSet)});
            }
        } else {
            Log::debug() << "Don't need to test right descendant"  << std::endl;
        }
        
    }
    
    Log::debug() << "Found " << toReturn.size() << " locations" << std::endl;
    
    // If we get here, we've finished all the DP jobs. Return all the
    // TextPositions we found where they were the unique overlap at some place
    // we explored in the DP table.
    return toReturn;

} 

