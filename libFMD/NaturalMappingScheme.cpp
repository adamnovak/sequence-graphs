#include "NaturalMappingScheme.hpp"
#include "Log.hpp"

#include <vector>

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
    
    // Then we can look at all the maximal unique matches and apply them to
    // actual query bases.
    
    // Then we can apply credit in another function.
    
    // Where will we keep unique matchings by base? We store them as sets of
    // TextPositions to which each base has been matched.
    std::vector<std::set<TextPosition>> matchings(query.size());
    
    // We also need a blacklist to note bases that were part of a maximal unique
    // match that was rejected. They can't be allowed to map inconsistently with
    // that maximal match, which translates into not being able to map at all
    // (since it was maximal) until credit comes along.
    std::vector<bool> blacklist(query.size());
    
    // Get the min and max matchings.
    std::vector<Matching> minMatchings = findMinMatchings(query);
    std::vector<Matching> maxMatchings = findMaxMatchings(query);
    
    for(Matching m : minMatchings) { 
        Log::trace() << "Min matching: " << m.start << " - " <<
            m.start + m.length << std::endl;
    }
    
    // We need to work out which min matchings belong to each max matching. This
    // holds how many min matchings have been used up and were contained in
    // previous (i.e. further right) max matchings. Every min matching is
    // contained within exactly one max matching (but may overlap others).
    size_t minMatchingsUsed = 0;
    
    for(Matching matching : maxMatchings) {
        Log::debug() << "Max matching " << matching.start << " - " <<
            matching.start + matching.length << " @ " << matching.location <<
            std::endl;
    
        // For each maximal unique match (always nonoverlapping) from right to
        // left...
        
        // See if this max matching passes the criteria for inclusion.
        
        while(minMatchings[minMatchingsUsed].start + 
            minMatchings[minMatchingsUsed].length > matching.start +
            matching.length) {
            
            Log::debug() << "\tDiscard min matching " << minMatchingsUsed <<
                ": " << minMatchings[minMatchingsUsed].start << " - " <<
                minMatchings[minMatchingsUsed].start +
                minMatchings[minMatchingsUsed].length << " @ " <<
                minMatchings[minMatchingsUsed].location << std::endl;
            
            // Throw away min matchings while they end further right than our
            // end.
            minMatchingsUsed++;
        }
        
        // How many min matchings have we grabbed that are contained in this max
        // matching? There will be at least one.
        size_t minMatchingsTaken = 0;
        
        // How many of those are non-overlapping?
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
        
        while(minMatchingsUsed + minMatchingsTaken < minMatchings.size() && 
            minMatchings[minMatchingsUsed + minMatchingsTaken].start >=
            matching.start) {
            
            // What's the next min matching?
            const Matching& minMatching = minMatchings[minMatchingsUsed +
                minMatchingsTaken];
            
            Log::debug() << "\tContains min matching " << 
                minMatchingsUsed + minMatchingsTaken << ": " <<
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
                nonOverlapping++;
                previousLeftEndpoint = minMatching.start;
                
                Log::debug() <<
                    "\t\tMatching is on maximal non-overlapping path" <<
                    std::endl;
            }
            
            // Take matchings while they start further right than our start,
            // until we run out.
            minMatchingsTaken++;
            
            
        }
        
        Log::debug() << "\tTook " << minMatchingsTaken << " min matches, " <<
            nonOverlapping << " non-overlapping" << std::endl;
        
        if(minMatchingsTaken == 0) {
            Log::critical() << "\tDoes not contain min matching " << 
                minMatchingsUsed + minMatchingsTaken << ": " <<
                minMatchings[minMatchingsUsed + minMatchingsTaken].start <<
                " - " <<
                minMatchings[minMatchingsUsed + minMatchingsTaken].start +
                minMatchings[minMatchingsUsed + minMatchingsTaken].length <<
                " @ " <<
                minMatchings[minMatchingsUsed + minMatchingsTaken].location <<
                std::endl;
                
            // There is clearly a problem if we find no contained minimal
            // matching.
            throw std::runtime_error(
                "Maximal matching contains no minimal matchings!");
        }
        
        // We contain minMatchings[minMatchingsUsed] to
        // minMatchings[minMatchings + minMatchingsTaken], inclusive on the low
        // end.
        
        // What's the total length of all included minimal exact matches?
        size_t totalMinLength = 0;
        
        for(size_t i = minMatchingsUsed;
            i < minMatchingsUsed + minMatchingsTaken; i++) {
            
            // For each min matching that we contain, add in its length.
            totalMinLength += minMatchings[i].length;
        }
        
        // Calculate the average length of involved minimal unique matches.
        double averageMinLength = (double) totalMinLength / 
            (double) minMatchingsTaken;
            
        Log::debug() << "\tAverage minimal unique match length: " <<
            averageMinLength << std::endl;
    
        // This flag keeps track of whether we passed all the filters on
        // admissible max contexts to map on.
        bool passedFilters = true;
    
        if(passedFilters && matching.length < minContext) {
            // Min context filter is on and we have failed it.
            passedFilters = false;
            
            Log::debug() << "\tFailed min context filter" << std::endl;
        }
        
        if(passedFilters && multContext * averageMinLength >= matching.length) {
            // Mult context filter is on and we have failed it; we aren't at
            // least multContext times as long as the average length of the min
            // unique matches we contain.
            passedFilters = false;
            
            Log::debug() << "\tFailed mult context filter: " <<
                multContext * averageMinLength << " > " << matching.length <<
                std::endl;
        }
        
        if(passedFilters) {
            // We passed all the filters. Add in matchings for the bases this
            // maximal unique match covers.
            
            for(size_t i = matching.start; i < matching.start + matching.length;
                i++) {
                
                Log::trace() << "Matching " << i << " to " <<
                    matching.location << std::endl;
                
                // For each position it covers, record the matching on that
                // position.
                matchings[i].insert(matching.location);
                
                // Bump the matching's location over by 1 for the next base.
                // We're working on the copy that the range for loop made rather
                // than the original.
                matching.location.addLocalOffset(1);
            }
        } else {
            // We failed at least one filter.
            
            // We neet to blacklist all the bases in this match, so they can't
            // map (without credit). Any mapping from some other mutually
            // overlapping maximal unique match would necessarily be
            // inconsistent with this maximal unique match, and we should
            // prevent that matching even if we can't justify our own.
            
            for(size_t i = matching.start; i < matching.start + matching.length;
                i++) {
                // Blacklist each involved query position.
                blacklist[i] = true;
            }
        }   
        
        // The next maximal matching has to use the minimal matchings that come
        // later.
        minMatchingsUsed += minMatchingsTaken;
    }
    
    // What mappings will we return? Make sure it is big enough; it should
    // default to unmapped.
    std::vector<Mapping> mappings(query.size());
    
    if(minMatchings.size() == 0) {
        // There were not any unique matchings actually. Just say everything is
        // unmapped.
        return mappings;
    }
    
    // Otherwise we know there is at least one unique match.
    
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
        }
    }
    
    // Now we've made all the mappings we need.
    return mappings;
}

void NaturalMappingScheme::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // Map using the natural context scheme: get matchings from all the
    // unique-in-the-reference strings that overlap you.
    
    // Map the query naturally.
    std::vector<Mapping> naturalMappings = naturalMap(query);
        
    // How many bases have we mapped or not mapped (credit or not).
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    
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
                    
                    Log::trace() << "Right credit zips " << i << " to " <<
                        implied << std::endl;
                    
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
                    
                    Log::trace() << "Left credit zips " << i << " to " <<
                        implied << std::endl;
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
    
    Log::output() << "Mapped " << mappedBases << " bases, " << 
        creditBases << " on credit, " << conflictedCredit << 
        " bases with conflicting credit." << std::endl;

    
}
