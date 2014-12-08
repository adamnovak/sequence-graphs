#include "NaturalMappingScheme.hpp"
#include "Log.hpp"

#include <vector>

std::vector<Mapping> NaturalMappingScheme::naturalMap(
    const std::string& query) const {
    
    // TODO: This could be dramatically simpler if we had a nice way to ban
    // things from mapping without making sure to never match.
    
    // What mappings will we return? Make sure it is big enough; it should
    // default to unmapped.
    std::vector<Mapping> mappings(query.size());
    
    // Where will we keep matchings? We store them as vectors of TextPositions
    // to which each base has been matched.
    std::vector<std::vector<TextPosition>> matchings(query.size());
    
    // When we retract bases, we'll call this to handle collating matchings into
    // mappings.
    std::function<void(size_t)> makeMapping = [&](size_t retracted) {
        Log::debug() << "Making mappings for " << retracted << std::endl;
        
        // We found all the strings that can overlap this retracted
        // character, so map it if possible.
        if(matchings[retracted].size() > 0) {
            // There were some matchings for this base
            
            // Do they all agree?
            bool matchingsAgree = true;
            
            for(size_t i = 1; i < matchings[retracted].size(); i++) {
                // For each of them after the first
                
                if(matchings[retracted][i] != matchings[retracted][0]) {
                    // Note if it disagrees
                    matchingsAgree = false;
                    break;
                }
            }
            
            if(matchingsAgree) {
                // Make a mapping to this place all the matchings agree on.
                mappings[retracted] = Mapping(matchings[retracted][0]);
            } else {
                // If they disagree we stay unmapped.
                mappings[retracted] = Mapping();
            }
        } else {
            // If there are no matchings we also stay unmapped.
            mappings[retracted] = Mapping();
        }
    };
    
    // Start with everything selected.
    FMDPosition results = index.getCoveringPosition();
    
    // What is the end (inclusive) of the shortest unique string up against the
    // left edge? Start past the end.
    int64_t leftmostMappable = query.size();
    
    // And the beginning of the shortest unique string up against the right
    // edge? Start past the beginning.
    int64_t rightmostMappable = -1;
    
    // How many characters are currently searched?
    size_t patternLength = 0;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--) {
        // For each position in the query from right to left
        
        // We're going to extend backward with this new base.
        FMDPosition extended = results;
        index.extendLeftOnly(extended, query[i]);
        
        
        while(extended.isEmpty(mask)) {
            // If you can't extend, retract until you can. TODO: Assumes we can
            // find at least one result for any character.
            
            // We're retracting this base, so nothing more can overlap it. Safe
            // to map.
            if((int64_t) i <= rightmostMappable) {
                // This base is lefter than the rightmost base that can't sneak
                // ambiguity out the right.
                
                // Make some mappings for it if it was consistently matched.
                makeMapping(i + patternLength);
            }
            
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
        
        if(extended.getLength(mask) == 1) {
            // We are currently unique.
            
            // Grab the BWT index of the only selected BWT position.
            // TODO: this is going to be super-hard to do with ranges.
            int64_t bwtIndex = extended.getResult(mask);
            
            // Work out the text position we have found. This corresponds to
            // the newly added character.
            TextPosition location = index.locate(bwtIndex);
            
            if(results.getLength(mask) == 1 && patternLength > minContext) {
                // We were previously unique, and we didn't just add the 1 base
                // to bring us up to minContext. We just need to add a matching
                // for this base.
                
                Log::debug() << "New unique position " << location << " len " << 
                    patternLength << " = " << i << std::endl;
                
                if(patternLength >= minContext) {
                    // This is long enough to care about.
                    
                    Log::debug() << "\tAdding matching " <<
                        matchings[i].size() << " for " << i << " = " <<
                        location << std::endl;
                    
                    // Say we match this query character to this text position.
                    matchings[i].push_back(location);
                }
                
            } else {
                // We are newly unique, or just added the 1 base to bring us up
                // to minContext.
                
                if((int64_t) i > rightmostMappable) {
                    // This is the rightmost (first) mappable base, since it
                    // can't get an ambiguous string to the right edge.
                    rightmostMappable = i;
                }
                
                Log::debug() << "Gained uniqueness at " << location << 
                    " len " << patternLength << " = " << i << std::endl;
                    
                Log::debug() << "Pattern length " << patternLength <<
                    " vs min context " << minContext <<
                    " vs rightmost mappable position " << rightmostMappable <<
                    std::endl;
                
                if(patternLength >= minContext) {
                    // This is long enough to care about.
                
                    for(size_t j = 0; 
                        j < patternLength && i + j <= rightmostMappable; j++) {
                        // Go through and add a matching for each selected base
                        // that isn't so far right it was involved in non-unique
                        // things touching the edge.
                            
                        // Make a new TextPosition to express where we are
                        // mapping to.
                        TextPosition mappedLocation = location;
                        mappedLocation.addLocalOffset(j);
                        
                        Log::debug() << "\tAdding matching " << 
                            matchings[i].size() << " for " << i + j << " = " <<
                            mappedLocation << std::endl;
                        
                        // Match to it.
                        matchings[i + j].push_back(mappedLocation);
                    }
                }
                
            }
            
        } else if(i + patternLength < query.size()) {
            // We are not unique, but we don't but up against the right edge.
            // This means that we can't get an ambiguous string from this new
            // base to the end of the query.
            
            if((int64_t) i > rightmostMappable) {
                //Note that we can match for this current base, if we ever find
                // a unique string for it.
                rightmostMappable = i;
            }
            
        }
        
        // Save the results after doing this position so we can do the next
        // query position.
        results = extended;
    }
    
    // Now we have made it to the left edge of the query, and seen all the
    // unique strings that matter.
    
    Log::debug() << "Doing final retractions" << std::endl;
    
    while(results.getLength(mask) == 1 && patternLength != 0) {
        // Retract until we aren't unique again. (We assume we are guaranteed to
        // at some point become non-unique.)
        
        if((int64_t) patternLength <= (int64_t) rightmostMappable) {
            // This base is lefter than the rightmost base that can't sneak
            // ambiguity out the right.
            
            // Make some mappings for it if it was consistently matched.
            makeMapping(patternLength);
        }
        
        // Retract that base and move on to the next one. TODO: Can we just get
        // the point at which we would get more results and do the whole run to
        // there at once somehow?
        FMDPosition retracted = results;
        index.retractRightOnly(retracted, --patternLength);
        
        // Update results with the retracted version.
        results = retracted;
    }
    
    Log::debug() << "Ended at " << patternLength << " with " <<
        results.getLength(mask) << " results." << std::endl;
        
    Log::debug() << "Non-unique substring: " <<
        query.substr(0, patternLength) << std::endl;
        
    Log::debug() << "Result count check: " <<
        LOG_LAZY(index.count(query.substr(0, patternLength)).getLength(mask)) <<
        std::endl;
    
    // Once we retract until the search is not unique, we're at a base that can
    // sneak ambiguity out the left, so we shouldn't map it.
    
    // Return the mappings.
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
