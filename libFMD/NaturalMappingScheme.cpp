#include "NaturalMappingScheme.hpp"
#include "Log.hpp"

#include <vector>
#include <algorithm>
#include <cstdlib>

std::vector<NaturalMappingScheme::Matching>
    NaturalMappingScheme::findMaxMatchings(const std::string& query) const {
 
    // What matchings have we found?
    std::vector<Matching> toReturn;
 
    // Start with everything selected.
    FMDPosition results = index.getCoveringPosition();
    
    // How many characters are currently searched?
    size_t patternLength = 0;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--) {
        // For each position in the query from right to left, we're going to
        // consider any maximal unique matches with left endpoints here.
        
        // We're going to extend backward with this new base.
        FMDPosition extended = results;
        index.extendLeftOnly(extended, query[i]);
        
        if(results.getLength(mask) == 1 && extended.isEmpty(mask)) {
            // We are already a maximal unique match and can't extend any more.
            // Report ourselves.
            // Note that we can only ever do this once per left endpoint.
            
            // Make sure we fix the endpoint as i + 1, since we moved i left
            // already and we want to talk about where it was last loop.
            toReturn.push_back(Matching(i + 1,  index.locate(results.getResult(
                mask)), patternLength));
        }
        
        while(extended.isEmpty(mask)) {
            // If you can't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPosition retracted = results;
            // Make sure to drop characters from the total pattern length.
            index.retractRightOnly(retracted, --patternLength);
            
            // Try extending again
            extended = retracted;
            index.extendLeftOnly(extended, query[i]);
            
            // Say that last step we came from retracted.
            results = retracted;
        }
    
        // We successfully extended.
        patternLength++;
        results = extended;
    }
    
    if(results.getLength(mask) == 1) {
        // We are a maximal unique match butted up against the left edge. Report
        // it.
        toReturn.push_back(Matching(0,  index.locate(results.getResult(
            mask)), patternLength));
    }
    
    // This works: if there were a longer unique match on the right, we would
    // not have retracted. And we explicitly check to see if there is a longer
    // unique match on the left.
    
    return toReturn;
}

std::vector<NaturalMappingScheme::Matching>
    NaturalMappingScheme::findMinMatchings(const std::string& query) const {

    // What matchings have we found?
    std::vector<Matching> toReturn;
 
    // Start with everything selected.
    FMDPosition results = index.getCoveringPosition();
    
    // How many characters are currently searched?
    size_t patternLength = 0;

    // Flag that says whether we need to retract at least once before we can
    // find another minimal unique match, since we just found one ending at a
    // certain place. TODO: Find a cleaner way to do this.
    bool mustRetract = false;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--) {
        // For each position in the query from right to left, we're going to
        // consider any minimal unique matches with left endpoints here.
        
        // Retract on the right until we can successfully extend on the left
        // without running out of results.
        
        // We're going to extend backward with this new base.
        FMDPosition extended = results;
        index.extendLeftOnly(extended, query[i]);
        
        while(extended.isEmpty(mask)) {
            // If you can't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPosition retracted = results;
            // Make sure to drop characters from the total pattern length.
            index.retractRightOnly(retracted, --patternLength);
            mustRetract = false;
            
            // Try extending again
            extended = retracted;
            index.extendLeftOnly(extended, query[i]);
            
            // Say that last step we came from retracted.
            results = retracted;
        }
        
        // Extend on the left.
        results = extended;
        patternLength++;
        
        // Retract on the right until the next retraction would make us not
        // unique, and report a minimal unique match starting at this position.
        FMDPosition retracted = results;
        index.retractRightOnly(retracted, patternLength - 1);
        
        while(retracted.getLength(mask) == 1) {
            // Retract until we would no longer be unique. Make sure to drop
            // characters from the total pattern length.
            results = retracted;
            index.retractRightOnly(retracted, (--patternLength) - 1);
            mustRetract = false;
        }
        
        if(results.getLength(mask) == 1 && retracted.getLength(mask) > 1 &&
            !mustRetract) {
            
            // We found a minimally unique match starting at this position and
            // ending patternLength right from here.
            toReturn.push_back(Matching(i,  index.locate(results.getResult(
                mask)), patternLength));
                
            // We can't find another minimal match until we move the right
            // endpoint.
            mustRetract = true;
        }
    }
    
    // When we get here, we already know we retracted as much as we could for
    // the leftmost extension, so there are no more results to report.
    
    // This works: if there were a shorter unique match on the right starting at
    // this position, we would have retracted to find it. And if there were a
    // shorter unique match on the left ending at this position, we would have
    // already reported it and not reported any more until we retracted.
    
    // Results will also come out in descending order by left endpoint.
    
    return toReturn;
}

std::vector<Mapping> NaturalMappingScheme::naturalMap(
    const std::string& query) const {
    
    // We need to find all the positions that map under the natural mapping
    // scheme without credit.
    
    // A base maps if:
    // It can't get to either end of the query without being unique.
    // All maximal unique matchings it is part of agree on where to place it.
    
    // So we can start at the end of the leftmost minimal unique match, and go
    // through the beginning of the rightmost (inclusive). Since anything
    // outside those can get to the end of the query without being unique.
    
    // Then we can look at all the maximal unique matches (divided into synteny
    // blocks and filtered) and apply them to actual query bases.
    
    // Then we can apply credit in another function.
    
    // How many bases are unmapped due to conflict?
    size_t unmappedByConflict = 0;
    
    // Get the min and max matchings.
    std::vector<Matching> minMatchings = findMinMatchings(query);
    std::vector<Matching> maxMatchings = findMaxMatchings(query);
    
    for(Matching m : minMatchings) { 
        Log::trace() << "Min matching: " << m.start << " - " <<
            m.start + m.length << " (+" << m.length << ")" << std::endl;
    }
    
    // We keep a map of SyntenyBlocks by (text, query offset) pairs.
    std::map<std::pair<size_t, int64_t>, SyntenyBlock> syntenyBlocks;
    
    // We need to build a graph from maximal Matchings to the previous maximal
    // Matchings they might be in runs with.
    std::vector<std::vector<size_t>> matchingGraph(maxMatchings.size());
    
    // We need to work out which min matchings belong to each max matching. This
    // holds how many min matchings have been used up and were contained in
    // previous (i.e. further right) max matchings. Every min matching is
    // contained within exactly one max matching (but may overlap others).
    size_t minMatchingsUsed = 0;
    
    for(size_t matchingNumber = 0; matchingNumber < maxMatchings.size();
        matchingNumber++) {
        
        // Go through all max matchings, and reference each.
        Matching& matching = maxMatchings[matchingNumber];
    
        if(matchingNumber != 0 && matchingNumber != maxMatchings.size() - 1 && 
            matching.length < ignoreMatchesBelow) {
            // This matching is not the first or the last, and it is too short
            // even to prevent mapping on bases it overlaps. Skip it entirely.
            // Pretend it isn't even there.
            Log::debug() << "\tSkipping matching entirely due to length " <<
                matching.length << " < " << ignoreMatchesBelow << std::endl;
            continue;
        }
        
        Log::debug() << "Max matching " << matching.start << " - " <<
            matching.start + matching.length << " (+" << matching.length <<
            ") @ " << matching.location << std::endl;
    
        // For each maximal unique match from right to left, see how many
        // minimal unique matchings it contains.
        
        while(minMatchings[minMatchingsUsed].start + 
            minMatchings[minMatchingsUsed].length > matching.start +
            matching.length) {
            
            Log::trace() << "\tDiscard min matching " << minMatchingsUsed <<
                ": " << minMatchings[minMatchingsUsed].start << " - " <<
                minMatchings[minMatchingsUsed].start +
                minMatchings[minMatchingsUsed].length << " (+ " <<
                minMatchings[minMatchingsUsed].length << ") @ " <<
                minMatchings[minMatchingsUsed].location << std::endl;
            
            // Throw away min matchings while they end further right than our
            // end.
            minMatchingsUsed++;
        }
        
        // How many min matchings have we grabbed that are contained in this max
        // matching? There will be at least one.
        size_t minMatchingsTaken = 0;
        
        // How many of those are non-overlapping? We know this maximal unique
        // match is at least this Hamming distance away from any other places in
        // the reference.
        size_t nonOverlapping = 0;
        
        // What was the left endpoint of the last non-overlapping minimal unique
        // matching we counted? We need this in order to implement the optimal
        // greedy algorithm for the "Activity Selection Problem"
        // <http://en.wikipedia.org/wiki/Activity_selection_problem>. Basically,
        // when reading from right to left, we can get the maximum number of
        // non-overlapping intervals by taking the non-overlapping interval with
        // the largest left endpoint (which, conveniently, we will encounter
        // first).
        size_t previousLeftEndpoint = query.size();
        
        while(minMatchingsUsed + matching.minMatchings < minMatchings.size() && 
            minMatchings[minMatchingsUsed + matching.minMatchings].start >=
            matching.start) {
            
            // What's the next min matching?
            const Matching& minMatching = minMatchings[minMatchingsUsed +
                matching.minMatchings];
            
            Log::trace() << "\tContains min matching " << 
                minMatchingsUsed + matching.minMatchings << ": " <<
                minMatching.start << " - " <<
                minMatching.start + minMatching.length << " @ " <<
                minMatching.location << std::endl;
            
            if(minMatching.start + minMatching.length <= previousLeftEndpoint) {
                // This min matching does not overlap the last non-overlapping
                // matching we collected. (Remember, it's start and length, so
                // start is inclusive and start+length is exclusive.) Since
                // we're taking min matchings in order of descending left
                // endpoint, it also must have the greatest left endpoint. So we
                // need to count it as non- overlapping and look only to the
                // left of it.
                matching.nonOverlapping++;
                previousLeftEndpoint = minMatching.start;
                
                Log::trace() <<
                    "\t\tMatching is on maximal non-overlapping path" <<
                    std::endl;
            }
            
            // Take matchings while they start further right than our start,
            // until we run out.
            matching.minMatchings++;
            
            // Make sure to count their length
            matching.minMatchingLength += minMatching.length;
        }
        
        // matching contains minMatchings[minMatchingsUsed] to
        // minMatchings[minMatchings + matching.minMatchings], inclusive on the
        // low end.
        
        Log::debug() << "\tTook " << matching.minMatchings <<
            " min matches (" << matching.minMatchingLength << "bp), " <<
            matching.nonOverlapping << " non-overlapping" << std::endl;
        
        if(matching.minMatchings == 0) {
            Log::critical() << "\tDoes not contain min matching " << 
                minMatchingsUsed + matching.minMatchings << ": " <<
                minMatchings[minMatchingsUsed].start << " - " <<
                minMatchings[minMatchingsUsed].start +
                minMatchings[minMatchingsUsed].length << " @ " <<
                minMatchings[minMatchingsUsed].location << std::endl;
                
            // There is clearly a problem if we find no contained minimal
            // matching.
            throw std::runtime_error(
                "Maximal matching contains no minimal matchings!");
        }
        
        // OK, now we know everything we need to know about this maximal unique
        // matching. We need to see which previous matchings it could possibly
        // be in a run with.
        
        for(size_t prevMatching = matchingNumber - 1; 
            prevMatching != (size_t) -1; prevMatching--) {
            
            // For each previous matching
            
            // Where do we end on the right?
            size_t ourEnd = matching.start + matching.length;
            // And where do they start on he left?
            size_t theirStart = maxMatchings[prevMatching].start;
            
            int64_t queryGapLength = (int64_t) theirStart - (int64_t) ourEnd;
            
            // TODO: These limits are asymetrical due to inclusive/exclusive
            // bounds.
            
            if(ourEnd > theirStart &&
                ourEnd - theirStart > maxAlignmentSize) {
             
                Log::debug() << matching << " can't come from " <<
                    maxMatchings[prevMatching] <<
                    " because it overlaps too much in query." << std::endl;
             
                // They overlap too far into us.
                continue;   
            }
            
            if(ourEnd <= theirStart &&
                theirStart - ourEnd > maxAlignmentSize) {
                
            
                // The distance between our right end and their left end is too
                // great.
                
                // If the distance in the query is too far, we can't connect at
                // all, and nothing further away can either.
                
                Log::debug() << matching << " can't come from " <<
                    maxMatchings[prevMatching] <<
                    " because it's too far upstream in query." << std::endl;
                
                break;
            }

            // We need to compare these two reference positions to see if these
            // could make any sense in the same run.
            TextPosition& curPos = matching.location;
            TextPosition& prevPos = maxMatchings[prevMatching].location;
            
            if(curPos.getText() == prevPos.getText() &&
                curPos.getOffset() < prevPos.getOffset()) {
            
                // They have to be on the same text and going in the correct
                // direction to chain, and they are.
                
                // Now we need to check to make sure the requested alignment
                // won't be too big on the referece side.
                TextPosition curEnd = curPos;
                curEnd.addLocalOffset(matching.length);
                
                // Pull out the reference offsets, which we can treat just lke
                // query ofsets for distance finding. Works because these
                // offsets are always on the stran (text) we are mapping to.
                ourEnd = curEnd.getOffset();
                theirStart = prevPos.getOffset();
                
                if(ourEnd > theirStart &&
                    ourEnd - theirStart > maxAlignmentSize) {
                 
                    Log::debug() << matching << " can't come from " <<
                        maxMatchings[prevMatching] <<
                        " because it overlaps too much in reference." <<
                        std::endl;
                 
                    // They overlap too far into us.
                    continue;   
                }
                
                if(ourEnd <= theirStart &&
                    theirStart - ourEnd > maxAlignmentSize) {
                    
                
                    // The distance between our right end and their left end is
                    // too great.
                    
                    Log::debug() << matching << " can't come from " <<
                        maxMatchings[prevMatching] <<
                        " because it's too far upstream in reference." <<
                        std::endl;
                    
                    continue;
                }
                
                // We'll have to look at the
                // bit between these two matchings.
                matchingGraph[matchingNumber].push_back(prevMatching);
                
                Log::debug() << matching << " can come from " << 
                    maxMatchings[prevMatching] << std::endl;
            } else {
                Log::debug() << matching << " can't come from " << 
                    maxMatchings[prevMatching] <<
                    " due to text/direction issue" << std::endl;
            }
            
        }
        
        // The next maximal matching has to use the minimal matchings that come
        // later.
        minMatchingsUsed += matching.minMatchings;
    }
    
    Log::info() << maxMatchings.size() << " maximal matchings exist." <<
        std::endl;
    
    // OK, now we have to do the DP, identify blocks that should create base to
    // base mappings, and make those.
    
    // Do the DP to find the good runs, and flag matchings in them.
    identifyGoodMatchingRuns(query, maxMatchings, matchingGraph);
    
    // Where will we keep unique matchings by base? We store them as sets of
    // TextPositions to which each base has been matched.
    std::vector<std::set<TextPosition>> matchings(query.size());
    
    // We also need a blacklist to note bases that were part of a maximal unique
    // match that was rejected. They can't be allowed to map inconsistently with
    // that maximal match, which translates into not being able to map at all
    // (since it was maximal) until credit comes along.
    std::vector<bool> blacklist(query.size());
    
    for(Matching& matching : maxMatchings) {
        // For each maximal matching, we may need to make some mappings.
    
        if(matching.canMatch) {
            // This mapping can produce matchings
    
            for(size_t i = matching.start;
                i < matching.start + matching.length; i++) {
                
                Log::trace() << "Matching " << i << " to " <<
                    matching.location << std::endl;
                
                // For each position it covers, record the matching on that
                // position.
                matchings[i].insert(matching.location);
                
                // Bump the matching's location over by 1 for the next base.
                // We're working on the copy that the range for loop made
                // rather than the original.
                matching.location.addLocalOffset(1);
                
                if(matching.blacklist) {
                    // We have to blacklist all these positions in the query
                    // so they can't map.
                    blacklist[i] = true;
                }
            }
        }
    }
    
    // Now we have to collate the per-base matchings and blacklists into
    // Mappings.
    
    // What mappings will we return? Use this vector. Make sure it is big
    // enough; it should default to unmapped.
    std::vector<Mapping> mappings(query.size());
    
    if(minMatchings.size() == 0) {
        // There were not any unique matchings actually. Don't bother finding
        // the places where we can get out. Just say everything is unmapped.
        return mappings;
    }
    
    // Otherwise we know there is at least one unique match, so we have to make
    // sure we aren't missing unique matches waiting to happen at the ends, by
    // not mapping the ends outside the inner edges of the outermost minimum
    // unique matchings.
    
    // What is the leftmost min unique matching? They are in order from right to
    // left.
    Matching leftmostMinMatching = minMatchings[minMatchings.size() - 1];
    
    // What is the rightmost min unique matching? They are in order from right
    // to left.
    Matching rightmostMinMatching = minMatchings[0];
    
    // What is the leftmost position that can't get a non-unique context over
    // the left edge?
    size_t leftmostMappable = leftmostMinMatching.start +
        leftmostMinMatching.length - 1;
        
    // What is the rightmost position that can't get a non-unique context over
    // the right edge?
    size_t rightmostMappable = rightmostMinMatching.start;
    
    if(unstable) {
        // Cheat out of the above and let us map right up to the edge if we have
        // a unique matching.
        leftmostMappable = 0;
        rightmostMappable = query.size() - 1;
    }
    
    Log::info() << "Able to map in range " << leftmostMappable << " - " <<
        rightmostMappable << std::endl;
    
    for(size_t i = leftmostMappable; i <= rightmostMappable; i++) {
        // For every mappable position between those two, inclusive.
        if(!blacklist[i] && matchings[i].size() == 1) {
            // If no matching for the position failed the filters, and the
            // position was matched to exactly one reference position...
            
            // TODO: Account for ranges here. Should be exactly one merged
            // position.
            
            // Grab that one TextPosition and make a Mapping to it for this
            // query base.
            mappings[i] = Mapping(*(matchings[i].begin()));
        } else if(matchings[i].size() > 1) {
            Log::debug() << "Conflict: " << matchings[i].size() <<
                " matchings for base " << i << std::endl;
            unmappedByConflict++;
        }
    }
    
    // Record the bases unmapped by conflict.
    stats.add("conflicted", unmappedByConflict);
    
    // Now we've made all the mappings we need.
    return mappings;
}

void NaturalMappingScheme::identifyGoodMatchingRuns(
    const std::string& query, std::vector<Matching>& maxMatchings,
    const std::vector<std::vector<size_t>>& matchingGraph) const {

    // This holds all the runs ending at a certain matching without getting too
    // many mismatches/edits.
    std::vector<std::vector<Run>> runs(maxMatchings.size());

    for(size_t i = 0; i < maxMatchings.size(); i++) {
        // For each matching from right to left
        
        // Get a reference to this matcing we're on.
        Matching& matching = maxMatchings[i];
        
        for(size_t previous : matchingGraph[i]) {
            // Consider arriving at matching i from matching previous.
            
            // Where does the previous block start (inclusive). This is also the
            // exclusive end of the gap.
            size_t prevStart = maxMatchings[previous].start;
            
            // How long is the gap between blocks i the query? This may be
            // negative if the blocks overlap due to a deletion.
            int64_t queryGapLength = prevStart - (int64_t) matching.start -
                (int64_t) matching.length;      
            
            // And now we need to find the same distance in the reference?
            TextPosition prevStartPos = maxMatchings[previous].location;
            TextPosition endPos = matching.location;
            endPos.addLocalOffset(matching.length);
            
            // Get the number of intervening (+) or shared (-) bases in the
            // reference. Works because these offsets are always on the stran
            // (text) we are mapping to.
            int64_t referenceGapLength = (int64_t) prevStartPos.getOffset() -
                (int64_t) endPos.getOffset();
            
            // How much would it cost to connect these blocks?
            size_t cost;
            
            if(queryGapLength == 0 && referenceGapLength == 0) {
                throw std::runtime_error(
                    "Grave error: Can't have query and reference matchings "
                    "both up against each other");
            } else if(queryGapLength <= 0 && referenceGapLength >= 0) {
                // They overlap in query and are separated in reference; charge
                // the bases that were deleted.
                cost = -queryGapLength + referenceGapLength;
            } else if(queryGapLength >= 0 && referenceGapLength <= 0) {
                // They overlap in the reference and are separated in the query;
                // charge the bases that were inserted.
                cost = queryGapLength - referenceGapLength;
            } else if(queryGapLength < 0 && referenceGapLength < 0) {
                // They overlap in both, but must be inconsistent. Charge the
                // difference.
                cost = std::abs(queryGapLength - referenceGapLength);
            } else {
                // They are separated by some distance in both.
                
                // Get a TextPosition to the first base in the gap
                TextPosition afterMatching = matching.location;
                afterMatching.addLocalOffset(matching.length);
            
                // The cost is the number of edits from there to the start of
                // the previous Matching. There can't be maxHammingDistance + 1
                // or more.
                cost = countEdits(query,
                    matching.start + matching.length, queryGapLength, afterMatching,
                    referenceGapLength, maxHammingDistance + 1);
            } 
            
            Log::debug() << "Arriving at " << maxMatchings[i] << " from " <<
                maxMatchings[previous]  << " (" << queryGapLength << "|" <<
                referenceGapLength << ") costs " << cost << std::endl;
            
            for(const Run& prevRun : runs[previous]) {
                // For each run ending with this other matching
                
                if(cost + prevRun.totalCost > maxHammingDistance) {
                    // Skip it if it costs to much to extend
                    continue;
                }
                
                // Otherwise, extend it with this new matching.
                Run newRun = prevRun;
                
                newRun.matchings.push_back(i);
                newRun.totalCost += cost;
                newRun.totalClearance += matching.nonOverlapping;
                
                // And add the extension as a run for this matching
                runs[i].push_back(newRun);
            }
        }
        
        // We also need a Run that is just this match
        Run startHere;
        startHere.matchings.push_back(i);
        startHere.totalCost = 0;
        startHere.totalClearance = matching.nonOverlapping;
        
        runs[i].push_back(startHere);
    }
    
    // Now we do a very inefficient thing: for each sufficiently good run, we
    // flag all the involved matchings as able to match up bases. This means
    // we'll possibly end up flagging the same matching repeatedly. TODO: make
    // this more efficient by somehow linking matchings up with there eventual
    // completed runs.
    for(size_t i = 0; i < maxMatchings.size(); i++) {
        // For each matching from right to left
        
        // Get a reference to this matcing we're on.
        Matching& matching = maxMatchings[i];
     
        for(Run& run: runs[i]) {
            // For each run that ends at this matching
            
            // Get the first and last matchings in the run
            Matching& firstMatching = maxMatchings[run.matchings[0]];
            Matching& lastMatching = maxMatchings[*(run.matchings.rbegin())];
                
            // Work out the length of the run
            size_t runLength = lastMatching.start + lastMatching.length -
                firstMatching.start;
            
            if(run.totalCost <= maxHammingDistance &&
                run.totalClearance >= minHammingBound &&
                runLength >= minContext) {
                
                Log::debug() << "Found acceptable " << runLength << 
                    "bp run of +" << run.totalClearance << ", -" << 
                    run.totalCost << ":" << std::endl;
                
                // The run has enough Hamming clearance, and doesn't have too
                // many mismatches/edits.
                
                for(size_t included : run.matchings) {
                    // For each matching in the run, flag it
                    maxMatchings[included].canMatch = true;
                    
                    Log::debug() << "\tContains matching " <<
                        maxMatchings[included] << std::endl;
                }
            
            } else {
                // Complain about insufficiently good runs.
                Log::debug() << "Unacceptable " << runLength << "bp run of +" <<
                     run.totalClearance << ", -" << run.totalCost << ":" <<
                     std::endl;
                for(size_t included : run.matchings) {
                    Log::debug() << "\tContains matching " <<
                        maxMatchings[included] << std::endl;
                }
            }
        }
    }
    
    // Now we have found all the possible good runs, and flagged all of their
    // matchings as good. 
    
    if(conflictBelowThreshold) {
    
        // But there are some maximal matchings that we did not flag as good,
        // but which we still want to cause conlict. So we blacklist those.
        
        for(Matching& matching : maxMatchings) {
            // For each matching
            if(!matching.canMatch && matching.length >= ignoreMatchesBelow) {
                // We can't use it, but we can't just ignore it, so make it
                // blacklisted.
                matching.blacklist = true;
                
                // Which means we need to mark it as important for making base-
                // to-base matchings.
                matching.canMatch = true;
            }
        }
    }
    
    // Now that we've flagged all the maximal unique matchings that can actually
    // produce useful base-to-base matchings, we need to look for the ones that
    // aren't flagged but might become flagged on extension, and which therefore
    // need to have their bases blacklisted so they don't map anywhere else
    // until that extension happens.
    
    // We're just going to blacklist any matchings that aren't flagged and that
    // are within maxAlignmentSize of the ends
        
    for(size_t i = 0; i < maxMatchings.size(); i++) {
        // For each matching coming in from the right
        
        // TODO: Only look at matchings near the ends somehow.
        
        if((maxMatchings[i].start < maxAlignmentSize || query.size() - 
            (maxMatchings[i].start + maxMatchings[i].length) <= 
            maxAlignmentSize) && !maxMatchings[i].canMatch && !unstable) {
        
            // We're too close to one end or the other, and aren't flagged as
            // part of a good enough run. And we care about stability.
            
            // This one could be in a run extending off the end, and it didn't
            // get flagged already.
            
            Log::debug() << "\t--- Blacklisting matching!" << std::endl;
            
            // This one needs to be blacklisted.
            maxMatchings[i].blacklist = true;
            // Which means we need to mark it as important for making base-to-
            // base matchings.
            maxMatchings[i].canMatch = true;
        
        }
    }
}

void NaturalMappingScheme::scan(SyntenyBlock& block,
    const std::string& query) const {
    
    // Scan through the SyntenyBlock, deciding what Mappings should make
    // matches, and what ones should be blacklisted.
    
    // We go through from right to left and do a DP algorithm. For each Mapping,
    // we look to see which previous Mappings it could have come from. We then
    // look at the runs for each, and tack on this new Mapping, retracting old
    // Mappings if the new Mapping won't fit. We also need to make sure to
    // deduplicate the runs, and to consider starting a new run here.
    
    // If this Mapping ends up in any run that passes our threshold, we flag it
    // as able to match up query bases.
    
    // We go along with an inchworm algorithm, sort of like we use to find the
    // matchngs in the first place: add as many as you can without going over
    // the mismatch threshold, and then, if you can't add one, throw off old
    // ones until you can.
    
    // Pull out a reference to the matching vector.
    std::vector<Matching>& matchings = block.maximalMatchings;
    
    // What's the index of the block that will be the next to remove? If it's
    // equal to the loop index, there aren't any left to remove.
    size_t nextToRemove = 0;
    
    // What's the index of the next matching that we have not yet flagged as
    // able to provide matchings? We use this so we don't have to fully scan
    // each run that makes matchings able to map.
    size_t nextToFlag = 0;
    
    // How many total non-overlapping minimal matches do we currently have
    // selected?
    size_t nonOverlapping = 0;
    
    // How many mismatches do we currently have selected?
    size_t mismatches = 0;
    
    for(size_t i = 0; i < matchings.size(); i++) {
        // We're going to extend with each block in turn.
        
        if(nextToRemove < i && 
            mismatches + matchings[i].mismatchesBefore > maxHammingDistance) {
            
            // We can retract (i.e. our range isn't empty), so we have to count
            // the mismatches before the next block, and extending with this
            // next block would put us over the max Hamming distance between the
            // query and the reference.
            
            Log::debug() << mismatches + matchings[i].mismatchesBefore << 
                " mismatches (+" << matchings[i].mismatchesBefore <<
                ") is too many (>" << maxHammingDistance << ")" << std::endl;
            
            
            Log::debug() << "Retracting matching " << nextToRemove <<
                " (" << matchings[nextToRemove].start << " - " <<
                matchings[nextToRemove].start +
                matchings[nextToRemove].length << "): -" <<
                matchings[nextToRemove + 1].mismatchesBefore <<
                " mismatches, -" << matchings[nextToRemove].nonOverlapping <<
                " minimal matches." << std::endl;
            
            // We can retract a block. Try that.
            nonOverlapping -= matchings[nextToRemove].nonOverlapping;
            
            if(nextToRemove != i - 1) {
                // We can drop the gap after the thing we are removing, too.
                mismatches -= matchings[nextToRemove + 1].mismatchesBefore;
            }
            
            // Next time we will remove the next thing, if we have to.
            nextToRemove++;
            
            // Try adding this matching again.
            i--;
            continue;                
        }
        
        Log::debug() << "Adding matching " << i << " (" << matchings[i].start << 
            " - " << matchings[i].start + matchings[i].length << "): +" <<
            ((nextToRemove < i) ? matchings[i].mismatchesBefore : 0) <<
            " mismatches, +" << matchings[i].nonOverlapping <<
            " minimal matches." << std::endl;
        
        // If we get here, we can add in this new matching.
        nonOverlapping += matchings[i].nonOverlapping;
        if(nextToRemove < i) {
            // We have the matching before this one, so we need to count the
            // gap.
            mismatches += matchings[i].mismatchesBefore;
        }
        
        // OK, we now have the matchings from nextToRemove to i, inclusive, and
        // we need to figure out if we pass the test.
        
        // Work out the length of the run, in bases. It's from the start of what
        // we just added on the left to the end of the thing we next have to
        // remove on the right.
        size_t runBases = matchings[nextToRemove].start +
            matchings[nextToRemove].length - matchings[i].start;
        
        
        Log::debug() << "Run (" << runBases << " bp): " << nextToRemove <<
            " - " << i << ": " << mismatches << " mismatches, " <<
            nonOverlapping << " minimal matches." << std::endl;
        
        if(nonOverlapping >= minHammingBound &&
            mismatches <= maxHammingDistance &&
            runBases > minContext) {
            
            // We're close enough to the reference and far enough from
            // everything else, and above the minimum context length. Flag
            // everything we have selected as able to provide matchings to
            // bases.
            
            for(size_t j = std::max(nextToRemove, nextToFlag); j <= i; j++) {
                // For each selected matching after the point where we have
                // flagged up to the newly added matching, flag it as able to
                // match up bases.
                matchings[j].canMatch = true;
                
                Log::debug() << "\t+++ Flagged " << matchings[j].start <<
                    " - " << matchings[j].start + matchings[j].length <<
                    " (+" << matchings[j].length <<
                    ") as useful for mapping!" << std::endl;
                
                // Note that we have flagged through here. Nothing un-flagged
                // and to the left of here needs to be flagged.
                nextToFlag = j + 1;
            }
        } else {
            // This run isn't good enough to let us actually use its component
            // matchings. But the matchings may get picked up by a logner run.
            Log::debug() << "\t--- Rejected!" << std::endl;
        }
    }
    
    if(conflictBelowThreshold) {
        // We're supposed to cause conflict for any maximal unique match
        // that doesn't make the cut.
        for(Matching& matching : matchings) {
            // For each maximal unique match
            if(!matching.canMatch) {
                // It can't match up bases, it needs to be blacklisted so
                // nothing else can take its bases.
                matching.canMatch = true;
                matching.blacklist = true;
            }    
        }
    }
    
    // Now that we've flagged all the maximal unique matchings that can actually
    // produce useful base-to-base matchings, we need to look for the ones that
    // aren't flagged but might become flagged on extension, and which therefore
    // need to have their bases blacklisted so they don't map anywhere else
    // until that extension happens.
    
    // We're just going to blacklist any matchings that aren't flagged and that
    // are within maxAlignmentSize of the ends
        
    for(size_t i = 0; i < matchings.size(); i++) {
        // For each matching coming in from the right
        
        // TODO: Only look at matchings near the ends somehow.
        
        if((matchings[i].start < maxAlignmentSize || query.size() - 
            (matchings[i].start + matchings[i].length) <= maxAlignmentSize) &&
            !matchings[i].canMatch && !unstable) {
        
            // We're too close to one end or the other, and aren't flagged as
            // part of a good enough run. And we cae about stability.
            
            // This one could be in a run extending off the end, and it didn't
            // get flagged already.
            
            Log::debug() << "\t--- Blacklisting matching!" << std::endl;
            
            // This one needs to be blacklisted.
            matchings[i].blacklist = true;
            // Which means we need to mark it as important for making base-to-
            // base matchings.
            matchings[i].canMatch = true;
        
        }
    }
}

size_t NaturalMappingScheme::countMismatches(const std::string& query,
    size_t queryStart, TextPosition referenceStart, size_t length, 
    int64_t threshold, char direction) const {

    // Count mismatches between both sequences starting at the inclusive start
    // positions and wandering off in the appropriate direction, until one or
    // the other ends, we see the threshold number of mismatches, or we hit the
    // length requested.
    
    Log::trace() << "Counting mismatches between " << referenceStart <<
        " and query index " << queryStart << " in " << length << " bases" <<
        std::endl;
    
    // Get the length of the reference text we are on.
    int64_t textLength = index.getContigLength(
        referenceStart.getContigNumber());
    
    size_t mismatchesFound = 0;
    
    // Count total bases looked at.
    size_t basesChecked = 0;
    
    // We're going to work on our arguments: advance referenceStart and
    // queryStart and knock down length.
    
    while(true) {
        if(referenceStart.getOffset() == (size_t) -1) {
            Log::trace() << "Ran off left of reference" << std::endl;
            break;
        }
        
        if(referenceStart.getOffset() >= textLength) {
            Log::trace() << "Ran off right of reference" << std::endl;
            break; 
        }
        
        if(queryStart == (size_t) -1) {
            Log::trace() << "Ran off left of query" << std::endl;
            break;
        }
        
        if(queryStart >= query.size()) {
            Log::trace() << "Ran off right of query" << std::endl;
            break;
        }
        
        if(length == (size_t) -1) {
            Log::trace() << "Ran out of length" << std::endl;
            break;
        }
        
        if(mismatchesFound >= threshold && threshold != -1) {
            Log::trace() << mismatchesFound << " hit threshold of " <<
                threshold << std::endl;
            break;
        }
        
        
        // Count the number of mismatches at this character.
        mismatchesFound +=
            (query[queryStart] != index.displayCached(referenceStart));
        
        // Update our arguments and tail-recurse
        length--;
        queryStart += direction;
        referenceStart.addLocalOffset(direction);
        
        basesChecked++;
    }
    
    Log::trace() << "Found " << mismatchesFound << " in " << basesChecked <<
        " bases checked." << std::endl;
    
    if((referenceStart.getOffset() == (size_t) -1 ||
        referenceStart.getOffset() >= textLength) && threshold != -1) {
        // We ran out of reference, and we have a threshold set. Running out of
        // reference is sufficiently bad that we will never be able to extend
        // anything over it with fewer than any number of mismatches. So we're
        // just going to return that threshold. TODO: this is kind of a hack,
        // just make our callers do their own bounds checking, or make a
        // canExtendOverEnd function.
        return threshold;
    }
    
    // We didn't fall off the reference and return our threshold, so return the
    // mismatches we actually found so far.
    return mismatchesFound;
}

size_t NaturalMappingScheme::countEdits(const std::string& query,
    size_t queryStart, size_t queryLength, TextPosition referenceStart,
    size_t referenceLength, int64_t threshold) const {
    
    // We can get more accurate costs for these for free, and keep 0s out of our
    // DP problem
    if(queryLength == 0) {
        // We have to delete every reference base.
        return referenceLength;
    }
    if(referenceLength == 0) {
        // We have to insert every query base.
        return queryLength;
    }
    
    if(((queryLength > referenceLength) &&
        (queryLength - referenceLength >= threshold)) ||
        ((referenceLength > queryLength) &&
        (referenceLength - queryLength >= threshold))) {
    
        // We know the difference in lengths is at least the threshold. In that
        // case, there must be at least threshold edits.
        return threshold;
        
    }
    
    if(queryLength > maxAlignmentSize ||
        referenceLength > maxAlignmentSize) {
    
        // This alignment would be bigger than the biggest one we want to do.
        Log::warning() << "Skipping " << queryLength << " x " <<
            referenceLength << " alignment" << std::endl;
            
        return threshold;
    }

    // Make the DP matrix for aligning these two strings. Indexed by query, then
    // by reference. Start it full of 0s.
    std::vector<std::vector<int>> matrix(queryLength + 1);
    for(std::vector<int>& row : matrix) {
        for(size_t i = 0; i < referenceLength + 1; i++) {
            row.push_back(0);
        }
    } 
    
    for(size_t i = 1; i < queryLength + 1; i++) {
        for(size_t j = 1; j < referenceLength + 1; j++) {
        
            // For each interior position, see if you should have a gap in one,
            // a gap in the other, or a match/mismatch
            
            // Don't advance in the query, and pay 1 for a gap.
            int queryGapCost = matrix[i - 1][j] + 1;
            // Don't advance in the reference, and pay 1 for a gap.
            int referenceGapCost = matrix[i][j - 1] + 1;
            
            // Advance in both, and pay 1 if this isn't a match.
            TextPosition referencePosition = referenceStart;
            referencePosition.addLocalOffset(j - 1);
            int matchMismatchCost = matrix[i - 1][j - 1] + 
                (query[queryStart + i - 1] !=
                index.displayCached(referencePosition));
                
            // Adopt the cost of whatever the best choice was.
            matrix[i][j] = std::min(std::min(queryGapCost, referenceGapCost),
                matchMismatchCost);
        }
    }
    
    // Read off the cost in the lower right
    return matrix[queryLength][referenceLength];
}

void NaturalMappingScheme::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // Map using the natural context scheme: get matchings from all the
    // unique-in-the-reference strings that overlap you.
    
    // Map the query naturally.
    std::vector<Mapping> naturalMappings = naturalMap(query);
        
    // How many bases have we mapped or not mapped (credit or not).
    size_t mappedBases = 0;
    
    // How many bases are mapped on credit?
    size_t creditBases = 0;
    // How many bases have conflicted credit?
    size_t conflictedCredit = 0;
    
    for(size_t i = 0; i < naturalMappings.size(); i++) {
        // For each query base
        
        if(naturalMappings[i].isMapped()) {
            // If it actually mapped...
            
            // Work out where to
            TextPosition candidate = naturalMappings[i].getLocation();
            
            // Dispatch the callback with this query index and this
            // TextPosition.
            callback(i, candidate);
                
            mappedBases++;
        }
    }
    
    if(credit) {
        // Each base zips matching bases inwards until it hits a mismatch,
        // or another mapped base. If zippings disagree, the base is not to
        // be mapped. Only bases between pairs of mapped bases can be
        // mapped.
        
        Log::info() << "Applying credit tolerating " << z_max <<
            " mismatches." << std::endl;
        
        // Make a set of all the zippings we will find for each position.
        std::vector<std::set<TextPosition>> zippings(
            naturalMappings.size());
    
        // Find the first mapped base (past the end if none)
        size_t firstMapped = 0;
        for(; firstMapped < naturalMappings.size() && 
            !naturalMappings[firstMapped].isMapped(); firstMapped++);
        
        // Find the last mapped base (size_t version of -1 if none)
        size_t lastMapped = naturalMappings.size() - 1;
        for(; lastMapped != (size_t) -1 && 
            !naturalMappings[lastMapped].isMapped(); lastMapped--);
        
        // This holds the index of the base currently providing credit.
        size_t provider;
        
        // How many mismatches have we seen since the credit provider?
        size_t mismatchesSeen;
        
        for(size_t i = firstMapped; i < naturalMappings.size() &&
            i <= lastMapped; i++) {
            
            // From the first to the last, do all the zipping off the right
            // sides of mapped bases.
            
            if(naturalMappings[i].isMapped()) {
                // This base is mapped and is now the credit provider
                provider = i;
                mismatchesSeen = 0;
            } else {
                // This base is not mapped. We know we saw a mapped base
                // already and thus have a credit provider.
                
                // What position would we be implied to be at?
                TextPosition implied = 
                    naturalMappings[provider].getLocation();
                implied.addLocalOffset((int64_t) i - (int64_t) provider);
                
                if(implied.getOffset() >= index.getContigLength(
                    implied.getContigNumber())) {
                    
                    // This position is implied to zip off the end of the
                    // reference text. Skip to the next one.
                    continue;
                }
                    
                if(index.displayCached(implied) != query[i]) {
                    // This is a mismatch
                    mismatchesSeen++;
                    
                    Log::trace() << index.displayCached(implied) << 
                        " at " << implied << " mismatches " << query[i] <<
                        std::endl;
                }
                
                if(mismatchesSeen <= z_max) {
                    // Not too many mismatches since the last credit
                    // provider. We can do credit.
                    
                    // TODO: not checking base identity here lets credit
                    // conflict even if it's wrong about base identity.
                    
                    // Zip this base to the position it is implied as being
                    // at by the credit provider.
                    zippings[i].insert(implied);
                    
                    Log::debug() << "Right credit zips " << i << " to " <<
                        implied << " from " << provider << std::endl;
                    
                }    
            }
            
        }
        
        for(size_t i = lastMapped; i != (size_t) -1 && i >= firstMapped;
            i--) {
            
            // From the last to the first, do all the zipping off the left
            // sides of mapped bases.
            
            // TODO: Unify this stateful logic with the above; it's
            // direction-independent.
            
            if(naturalMappings[i].isMapped()) {
                // This base is mapped and is now the credit provider
                provider = i;
                mismatchesSeen = 0;
            } else {
                // This base is not mapped. We know we saw a mapped base
                // already and thus have a credit provider.
                
                // What position would we be implied to be at?
                TextPosition implied = 
                    naturalMappings[provider].getLocation();
                implied.addLocalOffset((int64_t) i - (int64_t) provider);
                
                if(implied.getOffset() >= index.getContigLength(
                    implied.getContigNumber())) {
                    
                    // This position is implied to zip off the end of the
                    // reference text. Skip to the next one.
                    continue;
                }
                    
                if(index.displayCached(implied) != query[i]) {
                    // This is a mismatch
                    mismatchesSeen++;
                    
                    Log::trace() << index.displayCached(implied) <<
                        " at " << implied << " mismatches " << query[i] <<
                        std::endl;
                }
                
                if(mismatchesSeen <= z_max) {
                    // Not too many mismatches since the last credit
                    // provider. We can do credit.
                    
                    // TODO: not checking base identity here lets credit
                    // conflict even if it's wrong about base identity.
                    
                    // Zip this base to the position it is implied as being
                    // at by the credit provider.
                    zippings[i].insert(implied);
                    
                    Log::debug() << "Left credit zips " << i << " to " <<
                        implied << " from " << provider << std::endl;
                }    
            }
            
            // We don't need to do another pass, we can integrate here.
            
            if(zippings[i].size() == 1) {
                // This base got zipped to exactly one place.
                
                // Work out where to (the only element in the set)
                TextPosition candidate = *(zippings[i].begin());
                
                if(index.displayCached(candidate) == query[i]) {
                    // Dispatch the callback with this query index and this
                    // TextPosition.
                    callback(i, candidate);
                        
                    Log::debug() << "Credit agrees on " << i << std::endl;
                    
                    mappedBases++;
                    creditBases++;     
                } else {
                    // Merging to that position would merge mismatching
                    // bases.
                    Log::debug() << "Credit agrees on " << i <<
                        " but would merge a mismatch." << std::endl;
                }
                    
                
                    
            } else if(zippings[i].size() > 1) {
                // This base has conflicted credit.
                
                Log::debug() << "Credit disagrees on " << i << std::endl;
                
                conflictedCredit++;
            }
            
        }
    }
    
    // Track stats
    stats.add("mapped", mappedBases);
    stats.add("unmapped", query.size() - mappedBases);
    stats.add("credit", creditBases); // Mapped on credit
    stats.add("conflictedCredit", conflictedCredit); // Conflicted on credit
    
    Log::output() << "Mapped " << mappedBases << " bases, " << 
        creditBases << " on credit, " << conflictedCredit << 
        " bases with conflicting credit." << std::endl;
}
