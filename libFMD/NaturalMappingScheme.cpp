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
    
    // Then we can look at all the maximal unique matches (divided into synteny
    // blocks and filtered) and apply them to actual query bases.
    
    // Then we can apply credit in another function.
    
    // Get the min and max matchings.
    std::vector<Matching> minMatchings = findMinMatchings(query);
    std::vector<Matching> maxMatchings = findMaxMatchings(query);
    
    for(Matching m : minMatchings) { 
        Log::trace() << "Min matching: " << m.start << " - " <<
            m.start + m.length << " (+" << m.length << ")" << std::endl;
    }
    
    // We keep a map of SyntenyBlocks by (text, query offset) pairs.
    std::map<std::pair<size_t, int64_t>, SyntenyBlock> syntenyBlocks;
    
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
    
        // For each maximal unique match (always nonoverlapping) from right to
        // left...
        
        // See if this max matching is good enough to attach to the previous
        // one, or if we need to make a new SyntenyBlock.
        
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
        // matching. We need to throw it into the appropriate synteny block
        // based on the text (i.e. contig and orientation) to which it matches,
        // and the offset it implies between that text and the query.
        
        // What text are we mapped to?
        size_t text = matching.location.getText();
        
        // Where is the first base of the query string in that text? It may be
        // negative.
        int64_t queryOffset = (int64_t) matching.location.getOffset() -
            (int64_t) matching.start;

        Log::debug() << "\tGoing into synteny block " << text << ", " <<
            queryOffset << std::endl;
            
        // Go get the SyntenyBlock we should be adding to, or insert a new
        // default-constructed (i.e. empty) one if there isn't one already.
        SyntenyBlock& block = syntenyBlocks[std::make_pair(text, queryOffset)];
        
        if(block.maximalMatchings.size() > 0) {
            // We need to work out the mismatches in the gap between the end of
            // this new block and the start of the previous one.
            
            // Where does the last block start (inclusive). This is also the
            // exclusive end of the gap.
            size_t prevStart = block.maximalMatchings.rbegin()->start;
            // How long is the gap between blocks?
            size_t gapLength = prevStart - matching.start - matching.length;            
            
            // Get a TextPosition to the first base in the gap
            TextPosition afterMatching = matching.location;
            afterMatching.addLocalOffset(matching.length);
            
            // Go count up the mismatches in that gap.
            matching.mismatchesBefore = countMismatches(query,
                matching.start + matching.length, afterMatching, gapLength);
                
            Log::debug() << "\t" << matching.mismatchesBefore <<
                " associated mismatches" << std::endl;
        }
        
        // Add this matching and its statistics to the appropriate synteny
        // block.
        block.extendLeft(matching);
        
        // The next maximal matching has to use the minimal matchings that come
        // later.
        minMatchingsUsed += matching.minMatchings;
    }
    
    // OK, now we have to process those synteny blocks.
    
    // Where will we keep unique matchings by base? We store them as sets of
    // TextPositions to which each base has been matched.
    std::vector<std::set<TextPosition>> matchings(query.size());
    
    // We also need a blacklist to note bases that were part of a maximal unique
    // match that was rejected. They can't be allowed to map inconsistently with
    // that maximal match, which translates into not being able to map at all
    // (since it was maximal) until credit comes along.
    std::vector<bool> blacklist(query.size());
    
    for(auto& entry : syntenyBlocks) { 
        // Do the scan of each SyntenyBlock to figure out whether it matchings
        // should count, and whether they should be blacklisted.
        
        Log::debug() << "Scanning frame " << entry.first.first << ", " <<
            entry.first.second << std::endl;
        
        // Get a reference to the block, so we can update it in place.
        SyntenyBlock& block = entry.second;
        
        // Scan the block and set the canMatch and blacklist flags on its
        // matchings.
        scan(block, query);
        
        for(Matching matching : block.maximalMatchings) {
            // For each maximal matching in the block, we may need to make some
            // mappings.
        
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

void NaturalMappingScheme::scan(SyntenyBlock& block,
    const std::string& query) const {
    
    // Scan through the SyntenyBlock, deciding what Mappings should make
    // matches, and what ones should be blacklisted.
    
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
                matchings[nextToRemove].length << "): " <<
                matchings[nextToRemove].mismatchesBefore <<
                " mismatches, " << matchings[nextToRemove].nonOverlapping <<
                " minimal matches." << std::endl;
            
            // We can retract a block. Try that.
            nonOverlapping -= matchings[nextToRemove].nonOverlapping;
            
            if(nextToRemove != i - 1) {
                // We can drop the gap after the thing we are removing, too.
                mismatches -= matchings[nextToRemove].mismatchesBefore;
            }
            
            // Next time we will remove the next thing, if we have to.
            nextToRemove++;
            
            // Try adding this matching again.
            i--;
            continue;                
        }
        
        Log::debug() << "Adding matching " << i << " (" << matchings[i].start << 
            " - " << matchings[i].start + matchings[i].length << "): " <<
            matchings[i].mismatchesBefore << " mismatches, " <<
            matchings[i].nonOverlapping << " minimal matches." << std::endl;
        
        // If we get here, we can add in this new matching.
        nonOverlapping += matchings[i].nonOverlapping;
        if(nextToRemove < i) {
            // We have the matching before this one, so we need to count the
            // gap.
            mismatches += matchings[i].mismatchesBefore;
        }
        
        Log::debug() << "Run: " <<i << " (query pos " << matchings[i].start <<
            ") - " << nextToRemove << " (query pos " << 
            matchings[nextToRemove].start + matchings[nextToRemove].length <<
            ")" << std::endl;
            
        Log::debug() << "\t" << mismatches << " mismatches, " <<
            nonOverlapping << " non-overlapping minimal matches" << std::endl;
        
        // OK, we now have the matchings from nextToRemove to i, inclusive, and
        // we need to figure out if we pass the test.
        
        if(nonOverlapping >= minHammingBound &&
            mismatches <= maxHammingDistance) {
            
            // We're close enough to the reference and far enough from
            // everything else. Flag everything we have selected as able to
            // provide matchings to bases.
            
            for(size_t j = std::max(nextToRemove, nextToFlag); j <= i; j++) {
                // For each selected matching after the point where we have
                // flagged up to the newly added matching, flag it as able to
                // match up bases.
                matchings[j].canMatch = true;
                
                Log::debug() << "\tFlagged as useful for mapping" << std::endl;
                
                // Note that we have flagged through here. Nothing un-flagged
                // and to the left of here needs to be flagged.
                nextToFlag = j + 1;
            }
        }
    }
    
    // Now that we've flagged all the maximal unique matchings that can actually
    // produce useful base-to-base matchings, we need to look for the ones that
    // aren't flagged but might become flagged on extension, and which therefore
    // need to have their bases blacklisted so they don't map anywhere else
    // until that extension happens.
    
    // Grab the text and offset of the base before after first matching.
    TextPosition rightUnmatchedStart = matchings[0].location;
    rightUnmatchedStart.addLocalOffset(matchings[0].length);
    
    // Count the mismatches between the query starting to the right of
    // the first matching, and the reference starting to the right of
    // where the first matching maps. Count until we finish or get more
    // mismatches than are allowed with the ones we have selected. Count
    // left to right.
    size_t rightMismatches = countMismatches(query,
        matchings[0].start + matchings[0].length, rightUnmatchedStart,
        query.size() - matchings[0].length,
        maxHammingDistance + 1, 1);
        
    for(size_t i = 0; i < matchings.size(); i++) {
        // For each matching coming in from the right
        
        // Add in their mismatches (the rightmost will have 0).
        rightMismatches += matchings[i].mismatchesBefore;
        
        Log::debug() << rightMismatches << "/" << maxHammingDistance <<
            " mismatches at or right of query position " <<
            matchings[i].start + matchings[i].length << std::endl; 
        
        if(rightMismatches <= maxHammingDistance) {
            // This one needs to be blacklisted.
            matchings[i].blacklist = true;
            // Which means we need to mark it as important for making base-to-
            // base matchings.
            matchings[i].canMatch = true;
        } else {
            // We finally got too many mismatches before hitting this matching.
            // We can stop now.
            break;
        }
    }
    
    // And on the left? TODO: unify with above, somehow abstracting out
    // direction.
    TextPosition leftUnmatchedStart = matchings.rbegin()->location;
    leftUnmatchedStart.addLocalOffset(-1);
    
    // Count left from the left end of the last maximal unique matching, and see
    // how many mismatches there are (or if there are enough to block any new
    // matchings from connecting).
    size_t leftMismatches = countMismatches(query,
        matchings.rbegin()->start - 1, leftUnmatchedStart,
        matchings.rbegin()->start, maxHammingDistance + 1, -1);
    
    for(size_t i = matchings.size() - 1; i != (size_t) -1; i--) {
        // For each matching coming in from the left
        
        if(i != matchings.size() - 1) {
            // Add in their mismatches (the leftmost maximal unique matching
            // will not even have a corresponding field).
            leftMismatches += matchings[i + 1].mismatchesBefore;
        }
        
        Log::debug() << leftMismatches << "/" << maxHammingDistance <<
            " mismatches left of query position " <<
            matchings[i].start << std::endl; 
        
        if(leftMismatches <= maxHammingDistance) {
            // This one needs to be blacklisted.
            matchings[i].blacklist = true;
            // Which means we need to mark it as important for making base-to-
            // base matchings.
            matchings[i].canMatch = true;
        } else {
            // We finally got too many mismatches before hitting this matching.
            // We can stop now.
            break;
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
    
    // Get the length of the reference text we are on.
    int64_t textLength = index.getContigLength(
        referenceStart.getContigNumber());
    
    size_t mismatchesFound = 0;
    
    // We're going to work on our arguments: advance referenceStart and
    // queryStart and knock down length.
    
    while(
        // We're in-bounds on the reference
        referenceStart.getOffset() != (size_t) -1 &&
        referenceStart.getOffset() < textLength &&
        // We're in-bounds on the query
        queryStart != (size_t) -1 &&
        queryStart < query.size() &&
        // We haven't gotten to the end of length yet.
        length > 0 &&
        // We haven't found enough mismatches yet.
        mismatchesFound < threshold) {
        
        
        // Counjt the number of mismatches at this character.
        mismatchesFound +=
            (query[queryStart] != index.displayCached(referenceStart));
        
        // Update our arguments and tail-recurse
        length--;
        queryStart += direction;
        referenceStart.addLocalOffset(direction);
    }
    
    return mismatchesFound;
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
