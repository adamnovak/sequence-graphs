#include "NaturalMappingScheme.hpp"
#include "Log.hpp"

#include <vector>
#include <algorithm>
#include <cstdlib>

#include <boost/range/adaptor/reversed.hpp>

std::vector<Matching> NaturalMappingScheme::findMaxMatchings(
    const std::string& query) const {
 
    Log::info() << "Looking for max matchings" << std::endl;
 
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
        
        if(i < 100) {
            // See if we are somehow getting more results on extension.
            Log::info() << results.getLength(mask, ranges) <<
                " results before extending, " <<
                extended.getLength(mask, ranges) << " after" << std::endl;
        }
        
        if(results.isUnique(mask, ranges) && extended.isAmbiguous(mask, ranges)) {
            // We searched for something longer and got un-unique
            Log::critical() << "Got additional results by extending!" <<
                std::endl;
                
            Log::critical() << results.getLength(mask, ranges) <<
                " results before extending, " <<
                extended.getLength(mask, ranges) << " after" << std::endl;
                
            Log::critical() << "Query part: " <<
                query.substr(i, patternLength) << std::endl;
                
            Log::critical() << "Query range: " << i << " - " <<
                i + patternLength << std::endl;
                
            Log::critical() << "Before extending: " << results << std::endl;
            for(int64_t j = results.getForwardStart();
                j < results.getForwardStart() + results.getEndOffset() + 1;
                j++) {
                
                // Log each selected position before
                
                Log::critical() << index.display(j) << "\t" <<
                    index.locate(j) << "\t" << 
                    (mask == NULL ? -1 : mask->isSet(j)) << "\t" << 
                    (ranges == NULL ? -1 : ranges->isSet(j)) << "\t" <<
                    std::endl;
                
            }
            
            Log::critical() << "After extending: " << extended << std::endl;
            for(int64_t j = extended.getForwardStart();
                j < extended.getForwardStart() + extended.getEndOffset() + 1;
                j++) {
                
                // Log each selected position after
                
                Log::critical() << index.display(j) << "\t" <<
                    index.locate(j) << "\t" << 
                    (mask == NULL ? -1 : mask->isSet(j)) << "\t" << 
                    (ranges == NULL ? -1 : ranges->isSet(j)) << "\t" <<
                    std::endl;
                
            }
            
            Log::critical() << "Aborting run" << std::endl << std::flush;
            
            throw std::runtime_error("Got additional results by extending");
                
        }
        
        if(results.isUnique(mask, ranges) && extended.isEmpty(mask)) {
            // We are already a maximal unique match and can't extend any more.
            // Report ourselves.
            // Note that we can only ever do this once per left endpoint.
            
            // What BWT indices are selected?
            auto bwtIndices = results.getResults(mask);
            
            // We're going to find the lowest-text, lowest-offset correxponding
            // text position. Start out with the biggest a TextPosition can be.
            TextPosition smallest(std::numeric_limits<size_t>::max(),
                std::numeric_limits<size_t>::max());
            
            for(auto bwtIndex : bwtIndices) {
                // For each BWT index masked in (there must be at least 1), find
                // its corresponding TextPosition.
                TextPosition candidate = index.locate(bwtIndex);
                
                if(candidate < smallest) {
                    // If it has a lower text or offset, use it.
                    // Guaranteed to overwrite the original fake value.
                    smallest = candidate;
                }
            }
            
            // Make a Matching to this chosen place. Make sure we fix the left
            // endpoint in the query as i + 1, since we moved i left already and
            // we want to talk about where it was last loop.
            auto maxMatching = Matching(i + 1,  smallest, patternLength);
            
            // Send it off
            if(i < 100) {
                Log::info() << "Found max matching " << maxMatching << std::endl;
            }
            toReturn.push_back(maxMatching);
        } else if(results.isUnique(mask, ranges)) {
            // We're unique but not maximal on the left (extended is not empty
            // under mask)
            if(i < 100) {
                Log::info() << "Extending max matching left to " << i << " + " <<
                    patternLength << std::endl;
            }
        } else {
            // We aren't unique yet
            if(i < 100) {
                Log::info() << "Not yet unique at " << i << " + " << patternLength <<
                    std::endl;
            }
        }
        
        
        while(extended.isEmpty(mask)) {
            // If you can't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPosition retracted = results;
            // Make sure to drop characters from the total pattern length.
            index.retractRightOnly(retracted, --patternLength);
            
            if(i < 100) {
                Log::info() << "Retracting to " << patternLength << std::endl;
            }
            
            // Try extending again
            extended = retracted;
            index.extendLeftOnly(extended, query[i]);
            
            if(i < 100) {
                Log::info() << retracted.getLength(mask, ranges) <<
                    " results after retracting, " <<
                    extended.getLength(mask, ranges) << " after extending" <<
                    std::endl;
            }
            
            // Say that last step we came from retracted.
            results = retracted;
        }
    
        // We successfully extended.
        results = extended;
        // Increment the pattern length since we did actually extedn by 1. 
        patternLength++;
    }
    
    if(results.isUnique(mask, ranges)) {
        Log::info() << "Unique at left edge" << std::endl;
    
        // We are a maximal unique match butted up against the left edge. Report
        // it.
        // TODO: don't duplicate this code.
        // What BWT indices are selected?
        auto bwtIndices = results.getResults(mask);
        
        // We're going to find the lowest-text, lowest-offset correxponding
        // text position. Start out with the biggest a TextPosition can be.
        TextPosition smallest(std::numeric_limits<size_t>::max(),
            std::numeric_limits<size_t>::max());
        
        for(auto bwtIndex : bwtIndices) {
            // For each BWT index masked in (there must be at least 1), find
            // its corresponding TextPosition.
            TextPosition candidate = index.locate(bwtIndex);
            
            if(candidate < smallest) {
                // If it has a lower text or offset, use it.
                // Guaranteed to overwrite the original fake value.
                smallest = candidate;
            }
        }
        
        // Make a Matching to this chosen place.
        auto maxMatching = Matching(0,  smallest, patternLength);
        
        // Send it off
        Log::info() << "Found max matching " << maxMatching << std::endl;
        toReturn.push_back(maxMatching);
    } else {
        Log::info() << "Not unique at left edge" << std::endl;
    }
    
    // This works: if there were a longer unique match on the right, we would
    // not have retracted. And we explicitly check to see if there is a longer
    // unique match on the left.
    
    return toReturn;
}

std::vector<Matching>  NaturalMappingScheme::findMinMatchings(
    const std::string& query) const {

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
        
        while(retracted.isUnique(mask, ranges)) {
            // Retract until we would no longer be unique. Make sure to drop
            // characters from the total pattern length.
            results = retracted;
            index.retractRightOnly(retracted, (--patternLength) - 1);
            mustRetract = false;
        }
        
        if(results.isUnique(mask, ranges) &&
            retracted.isAmbiguous(mask, ranges) && !mustRetract) {
            
            // We found a minimally unique match starting at this position and
            // ending patternLength right from here. Just match it up against
            // any text that's selected. It would be a lot of work to go through
            // all the texts available to minimal matches, and we don't really
            // use their locations much anyway.
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

std::unordered_map<Matching, std::vector<std::pair<Matching, size_t>>>
    NaturalMappingScheme::generateMaxMatchingGraph(
    std::vector<Matching> maxMatchings, const std::string& query) const {
    
    // Bucket the matchings by text and diagonal. TODO: make this an
    // unordered_map and provide a way to hash pairs (which really seems like it
    // ought to be in the STL already...)
    std::map<std::pair<size_t, int64_t>, std::vector<Matching>> buckets;
    
    Log::info() << "Bucketing max matchings" << std::endl;
    
    for(const Matching& matching : maxMatchings) {
        // Claculate the diagonal: offset between the query and the reference
        int64_t diagonal = (int64_t) matching.start - 
            (int64_t) matching.location.getOffset();
            
        // Look up the vector for that text and that diagonal, and append this
        // Matching to it. Matchings will stay in query order.
        buckets[{matching.location.getText(), diagonal}].push_back(matching);
            
        // TODO: maybe make the buckets maps on start, end pairs to save some
        // copying around?
    }
    
    // Index each bucket for preceeding-range queries
    std::map<std::pair<size_t, int64_t>, IntervalIndex<Matching>> indices;
    
    Log::info() << "Indexing " << buckets.size() << " buckets" << std::endl;
    
    for(const auto& kv : buckets) {
        // For each bucket...
    
        // Make start, end keys for all the Matchings
        std::vector<std::pair<std::pair<size_t, size_t>, Matching>> withKeys;
        
        for(const Matching& matching : kv.second) {
            // For each matching, put it in the vector with its start, length
            // range.
            withKeys.push_back({{matching.start, matching.length}, matching});
        }
    
        // Make an IntervalIndex from the vector of mappings and keys. This
        // should be pre-sorted. Store it under the same text and diagonal as
        // the bucket.
        indices.emplace(kv.first, IntervalIndex<Matching>(withKeys));
    }
    
    // Now we have our index built. We need to turn this into a graph with
    // costs.
    std::unordered_map<Matching, std::vector<std::pair<Matching, size_t>>>
        graph;
        
    Log::info() << "Building graph" << std::endl;
    
    for(const Matching& matching : maxMatchings) {
        // For each matching
    
        // What are the text and diagonal?
        size_t text = matching.location.getText();
        int64_t diagonal = (int64_t) matching.start - 
            (int64_t) matching.location.getOffset();
    
        for(int64_t offset = - (int64_t) maxHammingDistance;
            offset <= (int64_t) maxHammingDistance; offset++) {
            // Scan the diagonals up and down by maxHammingDistance
            int64_t otherDiagonal = diagonal + offset;

            if(!indices.count({text, otherDiagonal})) {
                // Nothing else is on this diagonal.
                continue;
            }
            
            if(matching.start == 0 || 
                !indices[{text, otherDiagonal}].hasStartingBefore(
                matching.start - 1)) {
                
                // There is no matching on this diagonal starting strictly
                // before the matching we are interested in. We have to be
                // strict or we would find this matching before itself on its
                // own diagonal. TODO: Should we make this a strictly before
                // method and drop the -1?
                continue;
            }
            
            // Get the last matching starting before this one starts, in the
            // other diagonal.
            const Matching& previous = indices[{text, 
                otherDiagonal}].getStartingBefore(matching.start - 1).second;
            
            // How much are they separated in the query?            
            int64_t queryGapLength = (int64_t) matching.start -
                (int64_t) (previous.start + previous.length);
            
            // And in the reference (we know they are on the same text)? We can
            // just infer this from the diagonal offset and the query
            // separation.
            int64_t referenceGapLength = queryGapLength - offset;
            
            // How much should this gap cost?
            size_t gapCost;

            // We have 4 cases for the signs of these things.
            // +/+: do an alignment
            // +/-: Charge for indels to make up the length difference.
            // -/+: Charge for indels to make up the length difference.
            // -/-: Charge for indels to make up the length difference.
            
            // TODO: check overlaps against max alignment size.
            
            if(queryGapLength > 0 && referenceGapLength > 0) {
                // They are separated by some distance in both. They are on the
                // same text, so they must be syntenic.
                
                // Get a TextPosition to the first base in the gap
                TextPosition afterPrevious = previous.location;
                afterPrevious.addLocalOffset(previous.length);
            
                // The cost is the number of edits from there to the start of
                // the current Matching. There can't be maxHammingDistance + 1
                // or more, so fail fast if that looks like the case.
                gapCost = countEdits(query, previous.start + previous.length, 
                    queryGapLength, afterPrevious, referenceGapLength,
                    maxHammingDistance + 1);
                
            } else {
                // There is some overlap somewhere. They are still on the same
                // text and the previous one is before this one, so they are
                // still syntenic.
            
                // Charge for indels to make up the difference.
                gapCost = std::abs(queryGapLength - referenceGapLength);
                
                if(gapCost == 0) {
                    // This should never be 0, because there is overlap or a
                    // length difference.
                    throw std::runtime_error(
                        "Gap cost of 0 with no actual gap");
                }
            }
            
            if(gapCost <= maxHammingDistance) {
                // If the gap cost is low enough...
            
                // Record an edge from this matching to the previous one with
                // its cost.
                graph[matching].push_back({previous, gapCost});
                
            }
        }
        
        // Make sure to have a self edge at cost 0
        graph[matching].push_back({matching, 0});
    }
    
    // Give back the graph we have made.
    return graph;
}

std::unordered_map<Matching, std::vector<std::pair<Matching, size_t>>>
    NaturalMappingScheme::invertGraph(const std::unordered_map<Matching,
    std::vector<std::pair<Matching, size_t>>>& graph) const {
    
    // We need to invert the graph, so we need a new graph.
    std::unordered_map<Matching, std::vector<std::pair<Matching, size_t>>> inverted;
    
    for(const auto& kv : graph) {
        // For each source node and its edges
        for(const auto& edge : kv.second) {
            // For each {dest, cost} pair, add in the reverse edge
            inverted[edge.first].push_back({kv.first, edge.second});
        }
    }
    
    return inverted;
}

const Matching& NaturalMappingScheme::getMaxMatching(
    const IntervalIndex<Matching>& maxMatchings,
    const Matching& minMatching) const {
    
    // What max matching contains this min matching? It should be the last one
    // starting before it.
    
    if(!maxMatchings.hasStartingBefore(minMatching.start)) {
        // Somehow we have an orphan min matching.
        Log::critical() << "Min matching " << minMatching << 
            " starts too early to be in a max matching!" << std::endl <<
            std::flush;
            
        Log::critical() << "Next starting max matching is " <<
            maxMatchings.getStartingAfter(minMatching.start).second <<
            std::endl << std::flush;
            
        throw std::runtime_error(
            "Min matching starts too early");
    }

    // Find the matching it has to be    
    const Matching& maxMatching = maxMatchings.getStartingBefore(
        minMatching.start).second;
        
    if(maxMatching.start + maxMatching.length < 
        minMatching.start + minMatching.length) {
        
        // The max matching stops too early.
        
        Log::critical() << "Min matching " << minMatching << 
            " ends too late to be in a max matching!" << std::endl <<
            std::flush;
            
        Log::critical() << "Next ending max matching is " <<
            maxMatchings.getEndingAfter(minMatching.start +
            minMatching.length - 1).second << std::endl << std::flush;
            
        throw std::runtime_error(
            "Min matching ends too late");
    }
    
    Log::debug() << "Min matching " << minMatching << " has max matching " <<
        maxMatching << std::endl;
    
    // Return the matching we found.
    return maxMatching;
    
}

std::unordered_map<Matching, IntervalIndex<Matching>> 
    NaturalMappingScheme::assignMinMatchings(
    const IntervalIndex<Matching>& maxMatchings,
    const std::vector<Matching>& minMatchings) const {

    // We need to assign the min matchings to max matchings.
    
    // We keep a map from max matching to vector of min matchings with their
    // ranges, in ascending order.
    std::unordered_map<Matching,
        std::vector<std::pair<std::pair<size_t, size_t>, Matching>>>
        vectorForMax;
        
    Log::info() << "Assigning " << minMatchings.size() << 
        " min matchings to " << maxMatchings.size() << " max matchings" <<
        std::endl;
    
    for(const auto& minMatching : minMatchings) {
        // For each min matching, find the max matching it belongs to and put it
        // in that list. They come in in ascending order, so they also go out in
        // ascending order.
        
        // Find the max matching
        const Matching& maxMatching = getMaxMatching(maxMatchings, minMatching);
        // Give this min matching to it. TODO: have a helper function to
        // annotate matchings with their ranges like this, sine we do it in
        // several places and keep getting it wrong.
        vectorForMax[maxMatching].push_back({{minMatching.start, 
            minMatching.length}, minMatching});
    }
    
    // But we really need this map from max matching to IntervalIndex of min
    // matchings, instead of a map from max matching to vector of min matchings.
    std::unordered_map<Matching, IntervalIndex<Matching>> minsForMax;
    
    Log::info() << "Indexing " << vectorForMax.size() <<
        " min matching vectors" << std::endl;
    
    // Now we index the mins in each max.
    for(const auto& kv : vectorForMax) {
        // Just make an IntervalIndex of each vector.
        minsForMax.emplace(kv.first, IntervalIndex<Matching>(kv.second));
    }
    
    // Give back the assignments.
    return minsForMax;
}

std::unordered_map<Matching, std::vector<size_t>>
    NaturalMappingScheme::getMinMatchingChains(const std::unordered_map<
    Matching, std::vector<std::pair<
    Matching, size_t>>>& maxMatchingGraph, const std::unordered_map<
    Matching, IntervalIndex<
    Matching>>& minsForMax, const std::vector<
    Matching>& maxMatchings, bool isForward) const {
 
    Log::info() << "Constructing min matching chains (forward: " <<
        isForward << ")" << std::endl;
        
    for(const auto& kv : minsForMax) {
        Log::debug() << "Max " << kv.first << " has mins:" << std::endl;
        
        for(const auto& kv2 : kv.second) {
            Log::debug() << "\tMin: " << kv2.second << std::endl;
        }
    }
 
    // We have the max matching connectivity graph, the min matching
    // assignments, and the list of max matchings in order.
    
    // If we are going backwards, turn the graph around so we can come from max
    // matchings on the right instead.
    const auto& graph = isForward ? maxMatchingGraph : 
        invertGraph(maxMatchingGraph);
        
    // Make a DP table, storing achievable numbers of non-overlapping min
    // matchings in a chain by min matching and then by mismatch cost (<=
    // maxHammingDistance).
    std::unordered_map<Matching, std::vector<size_t>> table;
    
    for(size_t i = isForward ? 0 : maxMatchings.size() - 1;
        isForward ? i < maxMatchings.size() : i != (size_t) -1;
        isForward ? i++ : i--) {
        // For each max matching in the order we decided to go in (which would
        // be a whole lot easier to write if it were possible to decide between
        // std::vector::begin() and std::vector::rbegin() at runtime, I'll have
        // you know...)

        // Grab the matching
        const auto& maxMatching = maxMatchings[i];

        // Get the vector of all the min matchings in the current max matching,
        // in ascending order.
        const IntervalIndex<Matching>& mins = minsForMax.at(maxMatching);

        for(size_t j = isForward ? 0 : mins.size() - 1;
            isForward ? j < mins.size() : j != (size_t) -1;
            isForward ? j++ : j--) {
            // For each min matching in it, in the appropriate direction (at
            // least 1 must exist)
            const Matching& minMatching = mins[j].second;
            
            
            // Fill its DP table with 1 non-overlapping min match at any cost.
            table[minMatching] = std::vector<size_t>(maxHammingDistance + 1, 1);
            
            Log::debug() << "Initialized DP table for " << minMatching <<
                std::endl;
            
            for(const auto& edge : graph.at(maxMatching)) {
                // For each max matching we can pull from, and the cost to get
                // there (at least the self edge must exist)
                const Matching& prevMax = edge.first;
                size_t edgeCost = edge.second;
            
                // Get the index of all the min matchings in the previous max
                // matching, in ascending order.
                const IntervalIndex<Matching>& prevMins = minsForMax.at(
                    prevMax);
                    
                // There is only one min matching we need to consider coming
                // from: the last non-overlapping min matching in that max
                // matching (if we are going forward), or the first one (if we
                // are going backwards).
                
                // TODO: make hasEndingBefore and friends use open intervals.
                if(isForward) {
                    // We are going forward, so we need to make sure there's
                    // something left of this min matching in that max matching.
                    
                    if(minMatching.start == 0) {
                        Log::debug() << "Nothing in " << prevMax <<
                            " can end before our start of 0." << std::endl;
                        
                        // Skip on to the next graph edge.
                        continue;
                    }
                    
                    if(!prevMins.hasEndingBefore(minMatching.start - 1)) {
                        Log::debug() << "Nothing in " << prevMax <<
                        " ends at or left of " << minMatching.start - 1 <<
                        std::endl;
                        
                        // Skip on to the next graph edge.
                        continue;
                    }
                    
                } else {
                    // We are going backward, so we need to make sure there's
                    // something right of this min matching in that max
                    // matching.
                    
                    if(!prevMins.hasStartingAfter(minMatching.start + 
                        minMatching.length)) {
                    
                        Log::debug() << "Nothing in " << prevMax <<
                            " starts at or right of " << minMatching.start +
                            minMatching.length << std::endl;
                        
                        // Skip on to the next graph edge.    
                        continue;
                    
                    }
                    
                }
                
                // Get the min matching that we need to consider coming from.
                const Matching& prevMin = (isForward ? 
                    prevMins.getEndingBefore(minMatching.start - 1) : 
                    prevMins.getStartingAfter(minMatching.start + 
                    minMatching.length)).second;
            
                Log::debug() << "Could come to " << minMatching <<
                    " from " << prevMin << " with cost " << edgeCost <<
                    std::endl;
                        
                for(size_t k = edgeCost; k <= maxHammingDistance; k++) {
                    // Scan the range of our table that we can potentially
                    // fill in from the other table.
                    
                    // How long a run would we get if we incured total cost
                    // i and came from the corresponding place in the
                    // previous min match's table?
                    size_t newRunLength = table[prevMin][k - edgeCost] + 1;
                    
                    // If the new run is longer, use it.
                    table[minMatching][k] = std::max(table[minMatching][k],
                        newRunLength);
                        
                    // We will still need to fill in higher-cost entries
                    // later if we don't find anything better at the higher
                    // costs.
                }
            }
            
            // What's the longest run encountered so far scanning up through the
            // cost table?
            size_t longestRun = 0;
            for(size_t k = 0; k <= maxHammingDistance; k++) {
                // For each loaction in the cost table, see if we could instead
                // get a better run at a lower cost.
                
                // The table should have what's in it here if the longest run at
                // a lower cost, whichever is greater. And the longest run will
                // the what it was before or what we find here, whichever is
                // greater.
                longestRun = table[minMatching][k] = std::max(
                    table[minMatching][k], longestRun);
            }
        }
    }
    
    // Give back the filled DP table.
    return table;
}

std::set<Matching>
    NaturalMappingScheme::maxMatchingsInSyntenyRuns(
    const std::vector<Matching>& maxMatchings,
    const std::vector<Matching>& minMatchings,
    const std::string& query) const {

    // We need to make the connectivity graph, call the two directions of DP,
    // scan through the min matches, integrate the two directions, and flag the
    // corresponding max matchings. Then we need to make the set of max
    // matchings that pass.
    
    // Make the max matching cost graph
    auto graph = generateMaxMatchingGraph(maxMatchings, query);
    
    Log::info() << "Indexing all " << maxMatchings.size() << " max matchings" <<
        std::endl;
    
    // Make a vector of max matchings with their occupied intervals.
    std::vector<std::pair<std::pair<size_t, size_t>, Matching>> toIndex;
    for(const Matching& maxMatching : maxMatchings) {
        // Populate it with the inclusive start, end interval for each max
        // matching. Make sure we do (start, length) pairs as required.
        toIndex.push_back({{maxMatching.start, maxMatching.length},
            maxMatching});
    }
    
    if(toIndex.size() > 0) {
        Log::info() << "First max matching: " << (*toIndex.begin()).second <<
            std::endl;
        
        Log::info() << "Last max matching: " << (*toIndex.rbegin()).second <<
            std::endl;
    }
    
    // Index the max matchings by interval
    IntervalIndex<Matching> index(toIndex);
    
    // Assign all the min matchings to their max matchings, using the index we
    // just made.
    auto minsForMax = assignMinMatchings(index, minMatchings);
    
    for(const auto& kv : minsForMax) {
        Log::debug() << kv.first << " is a max." << std::endl;
    }
    
    // Do the DP looking left
    auto forwardChains = getMinMatchingChains(graph, minsForMax, maxMatchings,
        true);
        
    // And looking right
    auto reverseChains = getMinMatchingChains(graph, minsForMax, maxMatchings,
        false);
        
    // Make the set to return.
    std::set<Matching> toReturn;
    
    Log::info() << "Checking for passing max matchings" << std::endl;
    
    for(const Matching& maxMatching : maxMatchings) {
    
        Log::debug() << "Evaluating max matching " << maxMatching << std::endl;
    
        for(const auto& minMatchingRecord : minsForMax[maxMatching]) {
            // For each min matching, which we pull out of its index record...
            const Matching& minMatching = minMatchingRecord.second;
            
            // How long of a synteny chain is it in?
            size_t bestChain = 0;
            
            for(size_t forwardCost = 0; forwardCost <= maxHammingDistance;
                forwardCost++) {
                
                // For each way to apportion the mismatch cost between forward
                // and reverse directions...
                size_t reverseCost = maxHammingDistance - forwardCost;
                
                // How long a run can we get on that budget in each direction?
                // Subtract 1 so we don't double-count this min match.
                size_t runLength = forwardChains[minMatching][forwardCost] +
                    reverseChains[minMatching][reverseCost] - 1;
                   
               Log::debug() << "Cost " << forwardCost << "|" << reverseCost <<
                    ": chain with " <<
                    forwardChains[minMatching][forwardCost] << " + " <<
                    reverseChains[minMatching][reverseCost] << " - 1 = " <<
                    runLength << " min matchings" << std::endl;
                    
                bestChain = std::max(bestChain, runLength);
            }
            
            if(bestChain >= minHammingBound &&
                maxMatching.length >= minContext) {
                // This max matching has a min matching involved in a good
                // enough chain, and is itself sufficiently long, so keep it.
                toReturn.insert(maxMatching);
                
                Log::debug() << "Taking " << maxMatching << " with chain of " <<
                    bestChain << std::endl;
                
                // Skip the rest of the min matchings.
                break;
            } else {
                // We couldn't find a good enough chain through this min
                // matching.
                Log::debug() << "Best synteny chain involving " <<
                    minMatching << " is only " << bestChain << " long" <<
                    std::endl;
            }
                
        }
    
    }
    
    for(const auto& goodMatch : toReturn) {
        Log::debug() << goodMatch << " has a good enough synteny chain" <<
            std::endl;
    }
    
    // Give back the set of max matchings that have a min matching in a good
    // enough chain.
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
    
    // Flip them around to ascending order
    std::reverse(minMatchings.begin(), minMatchings.end());
    std::reverse(maxMatchings.begin(), maxMatchings.end());
    
    for(Matching m : minMatchings) { 
        Log::trace() << "Min matching: " << m.start << " - " <<
            m.start + m.length << " (+" << m.length << ")" << std::endl;
    }
    
    Log::info() << maxMatchings.size() << " maximal matchings exist." <<
        std::endl;
    
    // OK, now we have to do the DP, identify runs that should create base to
    // base mappings, and make those.
    
    // Which maximal matchings contain minimal matchings involved in
    // sufficiently good synteny runs?
    std::set<Matching> goodMaxMatchings = maxMatchingsInSyntenyRuns(
        maxMatchings, minMatchings, query);
    
    // Where will we keep unique matchings by base? We store them as sets of
    // TextPositions to which each base has been matched.
    std::vector<std::set<TextPosition>> matchings(query.size());
    
    // We also need a blacklist to note bases that were part of a maximal unique
    // match that was rejected. They can't be allowed to map inconsistently with
    // that maximal match, which translates into not being able to map at all
    // (since it was maximal) until credit comes along.
    std::vector<bool> blacklist(query.size());
    
    for(const Matching& matching : goodMaxMatchings) {
        // For each maximal matching in a sufficiently good block, make some
        // mappings.
    
        // Get a non-const TextPosition we can slide through the matching.
        TextPosition location = matching.location;
    
        for(size_t i = matching.start;
            i < matching.start + matching.length; i++) {
            
            Log::debug() << "Matching " << i << " to " << location << std::endl;
            
            // For each position it covers, record the matching on that
            // position.
            matchings[i].insert(location);
            
            // Bump the location over by 1 for the next base.
            location.addLocalOffset(1);
        }
    }
    
    for(const Matching& matching : maxMatchings) {
        // For each matching...
        if(!goodMaxMatchings.count(matching) &&
            matching.length >= ignoreMatchesBelow) {
            
            // We can't use it, because it doesn't have a MUS in a good enough
            // synteny run, but we can't just ignore it.
            
            // We need to check if it can reach the ends of the query
            // without too many mismatches, in which case we would need to
            // blacklist it to preserve stability.
            
            // TODO: Split this out into a function?
            
            // How many bases on the left of this MUM in the query do we care
            // about? We'll look up to half the max alignment size, since the
            // other aligned sequence has to be double this length.
            size_t leftQueryLength = std::min(matching.start,
                maxAlignmentSize/2);
            
            // Where do I have to start from in the reference for the end of the
            // query?
            TextPosition leftReferenceStart = matching.location;
            
            // And how many do we care about in the reference? We need 2x the
            // number in the query, or however many remain on the reference
            // text, whichever is smaller, in order to be sure of getting the
            // lowest-possible-cost one-side-justified alignment. TODO: how will
            // this ever work for graphs???
            size_t leftReferenceLength = std::min(2 * leftQueryLength,
                leftReferenceStart.getOffset());
            
            // Budge the reference start position over that much.
            leftReferenceStart.addLocalOffset(-leftReferenceLength);
            
            Log::debug() << "Doing alignment 0-" << leftQueryLength <<
                " against " << leftReferenceStart << " + " <<
                leftReferenceLength << std::endl;
            
            // What's the cost to reach out to the left end of the query? We
            // right justify the alignment, and we don't care about the exact
            // value if it's more than our maxHammingDistance.
            size_t leftEndCost = countEdits(query, 0,  leftQueryLength,
                leftReferenceStart, leftReferenceLength, maxHammingDistance + 1, 
                false, true);
            
            // TODO: We repeat most of this stuff in the opposite polarity for
            // the right side. Unify somehow?
            
            // Where does the bit after the MUM in the query start?
            size_t rightQueryStart = matching.start + matching.length;
            
            // How many bases on the right of this MUM in the query do we care
            // about? Again, we need to use double this from the reference, so
            // we go up to half the max size.
            size_t rightQueryLength = std::min(query.size() - rightQueryStart,
                maxAlignmentSize/2);
                
            // Where do I have to start from in the reference for the start of the
            // query?
            TextPosition rightReferenceStart = matching.location;
            rightReferenceStart.addLocalOffset(matching.length);
                
            // On the right side, we need again either twice the query length,
            // or however much is available.
            size_t rightReferenceLength = std::min(2 * rightQueryLength, 
                index.getContigLength(rightReferenceStart.getContigNumber()) -
                rightReferenceStart.getOffset());
            
            Log::debug() << "Doing alignment " << rightQueryStart << "-" << 
                rightQueryStart + rightQueryLength << " against " <<
                rightReferenceStart << " + " << rightReferenceLength <<
                std::endl;
            
            // What's the cost to reach out to the right end of the query? We
            // left justify the alignment, and we don't care about the exact
            // value if it's more than our maxHammingDistance.
            size_t rightEndCost = countEdits(query, rightQueryStart, 
            rightQueryLength, rightReferenceStart, rightReferenceLength,
            maxHammingDistance + 1, true, false);
            
            Log::debug() << matching << " has end costs " << leftEndCost <<
                " in " << leftQueryLength << "|" << leftReferenceLength <<
                " bases on the left and " << rightEndCost << " in " <<
                rightQueryLength << "|" << rightReferenceLength <<
                " bases on the right" << std::endl;
            
            // TODO: what if there are plenty of differences, but not in range
            // given the sizer of the alignment we are willing to do? Should we
            // check over at the far ends instead or something?
            
            if((!unstable && (leftEndCost <= maxHammingDistance ||
                rightEndCost <= maxHammingDistance)) ||
                conflictBelowThreshold) {
                
                // There aren't enough mismatches between here and the ends of
                // the query, or we've opted to blacklist every little MUM that
                // doesn't pass.
                
                
                
                for(size_t i = matching.start;
                    i < matching.start + matching.length; i++) {
                    
                    // Blacklist every query position in the matching.
                    blacklist[i] = true;
                    
                    Log::debug() << "Blacklisting base " << i << 
                        " for participation in max matching " << matching <<
                        std::endl;
                    
                    // TODO: Make this count as conflict if it stops any bases
                    // from mapping.
                    
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
    
    // What is the leftmost min unique matching? They are in order from left to
    // right.
    Matching leftmostMinMatching = minMatchings[0];
    
    // What is the rightmost min unique matching? They are in order from left to
    // right.
    Matching rightmostMinMatching = minMatchings[minMatchings.size() - 1];
    
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

size_t NaturalMappingScheme::countEdits(const std::string& query,
    size_t queryStart, size_t queryLength, TextPosition referenceStart,
    size_t referenceLength, int64_t threshold, bool leftJustify,
    bool rightJustify) const {
    
    // We can get more accurate costs for these for free, and keep 0s out of our
    // DP problem
    if(queryLength == 0) {
        // We have to delete every reference base.
        Log::debug() << "Alignment is just to delete the reference" <<
            std::endl; 
        // This costs 1 per reference base, if we're justifying the query to
        // both ends. If we let the query not make it to one end of the
        // reference, then this is free.
        return (leftJustify && rightJustify) ? referenceLength : 0;
    }
    if(referenceLength == 0) {
        // We have to insert every query base.
        Log::debug() << "Alignment is just to insert the query" << std::endl; 
        return queryLength;
    }
    
    if((leftJustify && rightJustify) && (((queryLength > referenceLength) &&
        (queryLength - referenceLength >= threshold)) ||
        ((referenceLength > queryLength) &&
        (referenceLength - queryLength >= threshold)))) {
    
        // We are justifying on both ends, and we know the difference in lengths
        // is at least the threshold. In that case, there must be at least
        // threshold edits.
        Log::debug() << 
            "Alignment length difference so great we don't have to do it" << 
            std::endl;
        return threshold;
        
    }
    
    if(queryLength > maxAlignmentSize ||
        referenceLength > maxAlignmentSize) {
    
        // This alignment would be bigger than the biggest one we want to do.
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
    
    // Initialize the first row and column of the matrix with the costs for
    // reference and query gaps on the left edge.
    for(size_t i = 0; i < queryLength + 1; i++) {
        // How much is a gap in the reference at the left (unbused query bases)?
        // Always 1.
        matrix[i][0] = i;
    }
    for(size_t j = 0; j < referenceLength + 1; j++) {
        // How much is a gap in the query at the left (unused reference bases)?
        // Depends on if we're left-justifying or not.
        matrix[0][j] = leftJustify * j;
    }
    
    // Now we computew and dump our DP matrix.
    Log::debug() << "Actually doing alignment" << std::endl;
    
    auto& headerLine = Log::debug() << std::setw(3) << "#" << std::setw(3) <<
        "#";
    for(size_t j = 1; j < referenceLength + 1; j++) {
        // Dump all the reference characters as a line
        TextPosition referencePosition = referenceStart;
        referencePosition.addLocalOffset(j - 1);
        headerLine << std::setw(3) << index.displayCached(referencePosition);
    }
    headerLine << std::endl;
    
    auto& firstRow = Log::debug() << std::setw(3) << "#";
    for(size_t j = 0; j < referenceLength + 1; j++) {
        // Dump the first row of the matrix.
        firstRow << std::setw(3) << matrix[0][j];
    }
    firstRow << std::endl;
    
    for(size_t i = 1; i < queryLength + 1; i++) {
    
        auto& matrixLine = Log::debug() << std::setw(3) <<
            query[queryStart + i - 1] << std::setw(3) << matrix[i][0];  
        
        for(size_t j = 1; j < referenceLength + 1; j++) {
        
            // For each interior position, see if you should have a gap in one,
            // a gap in the other, or a match/mismatch
            
            // Don't advance in the reference, and pay 1 for a gap.
            int referenceGapCost = matrix[i - 1][j] + 1;
            
            // Don't advance in the query, and pay 1 for a gap.
            int queryGapCost = matrix[i][j - 1] + 1;
            
            if(!rightJustify && i == queryLength) {
                // Don't charge the +1 for gaps in the query (unused reference
                // bases) on the right end.
                queryGapCost = matrix[i][j - 1];
            }
            
            // Advance in both, and pay 1 if this isn't a match.
            TextPosition referencePosition = referenceStart;
            referencePosition.addLocalOffset(j - 1);
            int matchMismatchCost = matrix[i - 1][j - 1] + 
                (query[queryStart + i - 1] !=
                index.displayCached(referencePosition));
                
            // Adopt the cost of whatever the best choice was.
            matrix[i][j] = std::min(std::min(queryGapCost, referenceGapCost),
                matchMismatchCost);
                
            matrixLine << std::setw(3) << matrix[i][j];
        }
        matrixLine << std::endl;
    }
    
    Log::debug() << "Alignment done!" << std::endl;
    
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
