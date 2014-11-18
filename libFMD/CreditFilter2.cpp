#include "CreditFilter2.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>
#include "Log.hpp"
#include "MismatchResultSet.hpp"

/**
 * Given the forward strand of a reference text, a flag for whether we are
 * talking about the reverse strand, and a 0-based start position on that
 * strand, along with a query string and a 0-based start position in that, count
 * the number of mismatches in the next length characters.
 *
 * Note that no bounds checking is done.
 */
size_t countMismatches(const std::string& reference, bool referenceStrand, 
    size_t referenceStart, const std::string& query, size_t queryStart, 
    size_t length) {
    
    // We start with no mismatches.
    size_t mismatches = 0;
    
    for(size_t i = 0; i < length; i++) {
        // For each comparison we need to make
        
        // What character do we want from the reference?
        char refChar;
        if(referenceStrand) {
            // Pull from the reverse strand of the reference.
            refChar = complement(reference[reference.size() - 1 - 
                referenceStart - i]);
        } else {
            // Pull from the forward strand
            refChar = reference[referenceStart + i];
        }
        
        // Add in any mismatches.
        mismatches += (refChar != query[queryStart + i]);
    }

    // Return the total number of mismatches.    
    return mismatches;
}

CreditFilter2::CreditFilter2(const FMDIndex& index, 
    const GenericBitVector& ranges, size_t z_max, const GenericBitVector* mask):
    index(index), ranges(ranges), z_max(z_max), mask(mask),
    disambiguate(index) {
    
    // Nothing to do!
}

std::vector<Mapping> CreditFilter2::apply(
    const std::vector<Mapping>& leftMappings,
    const std::vector<Mapping>& rightMappings, const std::string& query) {

    // First disambiguate everything. This also makes combined left and right
    // contexts.
    std::vector<Mapping> disambiguated = disambiguate.apply(leftMappings,
        rightMappings);

    // Find the left and right sentinels (must be mapeped on the left/right
    // respectively, and not have the other side disagree). Maybe this can be
    // solved by disambiguation.
    
    int64_t leftSentinel = -1;
    int64_t rightSentinel = -1;
    
    // We need to track the lengths of the sentinel words for later
    size_t leftWordLength = 0;
    size_t rightWordLength = 0;
    
    for(size_t i = 0; i < disambiguated.size(); i++) {
        // Look in from the left for the left sentinel
        if(leftMappings[i].isMapped() && 
            disambiguated[i].isMapped()) {
            
            // This is a left-mapped thing on the left, and may be a sentinel.
            // To be sure, we have to find its word.
            
            // Get the length of its longest possible word, which has its right
            // end at i
            leftWordLength = disambiguated[i].getLeftMaxContext();
            
            // Clip out the word that the base mapped on on the left.
            std::string word = query.substr(i + 1 - leftWordLength,
                leftWordLength);
            
            Log::debug() << "Considering left sentinel word " << word <<
                std::endl;
            
            // Count how many times the word appears within the right number of
            // mismatches.
            // Count how many times the word appears within the right number of
            // mismatches.
            MismatchResultSet count = index.mismatchCount(ranges, word,
                z_max, mask);
            
            if(count.isMapped(ranges, mask)) {
                // This word appears exactly once within the specified number of
                // mismatches. So this is the leftmost left sentinel.
                leftSentinel = i;
                Log::info() << "Left sentinel found at " << leftSentinel << 
                    std::endl;
                break;
            } else {
                Log::debug() << "Left word candidate had " <<
                    count.getLength(mask) << " positions under mask " << mask <<
                    " with " << z_max << " mismatches." << std::endl;
            }
        }
    }
    
    for(size_t i = disambiguated.size() - 1; i != (size_t) -1; i--) {
        // Look in from the right for the right sentinel
        if(rightMappings[i].isMapped() && 
            disambiguated[i].isMapped()) {
            
            // This is a right-mapped thing on the right, and may be a sentinel.
            // To be sure, we have to find its word.
            
            // Get the length of its longest possible word, which has its left
            // end at i
            rightWordLength = disambiguated[i].getRightMaxContext();
            
            // Clip out the word that the base mapped on on the right.
            std::string word = query.substr(i, rightWordLength);
            
            Log::debug() << "Considering right sentinel word " << word <<
                std::endl;
            
            // Count how many times the word appears within the right number of
            // mismatches.
            MismatchResultSet count = index.mismatchCount(ranges, word,
                z_max, mask);
            
            if(count.isMapped(ranges, mask)) {
                // This word appears exactly once within the specified number of
                // mismatches. So this is the rightmost right sentinel.
                rightSentinel = i;
                Log::info() << "Right sentinel found at " << rightSentinel << 
                    std::endl;
                break;
            } else {
                Log::debug() << "Right word candidate had " <<
                    count.getLength(mask) << " positions under mask " << mask <<
                    " with " << z_max << " mismatches." << std::endl;
            }
        }
    }
    
    if(leftSentinel == -1 || rightSentinel == -1 || 
        rightSentinel <= leftSentinel) {
        
        // We didn't find one of the sentinels, or there is no space between
        // them. We can't do credit between them. Just disambiguate without
        // applying credit.
        Log::info() << "No sequence between sentinels. No credit applied." <<
            std::endl;
        return disambiguated;
    }
    
    // Make a function that tells us whether a base's maximal left context can
    // extend over the left sentinel's word.
    auto extendsOverLeft = [&](int i) -> bool {
        // The max context length counts the base being tested as 1. The base -
        // the left sentinel gets the number of characters between them (not
        // including one end), and then the right word length (which includes
        // the left sentinel) covers the left sentinel. So if the left max
        // context reaches back to exactly the left end of the left word, these
        // will be equal.
        return disambiguated[i].getLeftMaxContext() >= 
            i - leftSentinel + leftWordLength;
    };
    
    // And another to see if a base's maximal right context can extend over the
    // right sentinel's word.
    auto extendsOverRight = [&](int i) -> bool {
        return disambiguated[i].getRightMaxContext() >= 
            rightSentinel - i + rightWordLength;
    };
    
    // Get the positions for the two sentinels
    TextPosition leftPosition = disambiguated[leftSentinel].getLocation();
    TextPosition rightPosition = disambiguated[rightSentinel].getLocation();
    
    // The sentinels are a pair if neither can get a context over the other, or
    // if they are consistent. So they are not a pair if one can get a context
    // over the other and they are not consistent (i.e. map to places other than
    // the offset between them suggests).
    
    if((extendsOverRight(leftSentinel) || extendsOverLeft(rightSentinel)) && 
        !leftPosition.isConsistent(rightPosition, 
        rightSentinel - leftSentinel)) {
    
        // One can get over the other, and they aren't the right distance apart
        // on the same strand. This is not a sentinel pair.
        
        // If the outermost sentinels are not a pair, no sentinel pair exists on
        // this strand, so we can give up on credit.
        Log::warning() << "Sentinels don't pair! No credit applied." <<
            std::endl;
        return disambiguated;
    }
    
    // Now we need to figure out which bases are credit providers. Mapped bases
    // between sentinels that are either consistent with or unable to reach over
    // each sentinel are credit providers.
    
    // This holds which bases in the string are credit providers
    std::vector<bool> isCreditProvider(disambiguated.size(), false);
    
    // How many are there?
    size_t totalCreditProviders = 0;
    
    for(size_t i = leftSentinel; i <= rightSentinel; i++) {
        // For each position between the sentinels (including them)...
        
        if(!disambiguated[i].isMapped()) {
            // Not mapped, so doesn't count
            continue;
        }
        
        // Where did the base map?
        TextPosition basePosition = disambiguated[i].getLocation();
        
        if(extendsOverLeft(i) && 
            !leftPosition.isConsistent(basePosition, i - leftSentinel)) {
            
            // The maximal left context extends over the left sentinel word, but
            // our mapping is inconsistent with the left sentinel. Can't give
            // credit; we might unmap.
            continue;
        }
        
        if(extendsOverRight(i) && 
            !basePosition.isConsistent(rightPosition, rightSentinel - i)) {
            
            // The maximal right context extends over the right sentinel word,
            // but our mapping is inconsistent with the right sentinel. Can't
            // give credit; we might unmap.
            continue;
        }
        
        // If we get here, we can give credit.
        isCreditProvider[i] = true;
        totalCreditProviders++;
        
        Log::trace() << "Position " << i << " is a credit provider." <<
            std::endl;
        
    }
    
    // The sentinels always have to come out as credit providers.
    assert(isCreditProvider[leftSentinel]);
    assert(isCreditProvider[rightSentinel]);
    
    Log::info() << "Assigning credit from " << totalCreditProviders <<
        " providers." << std::endl;
    
    // Now we can actually try applying credit.
    
    // The obvious way is for each unmapped base to go looking for credit
    // providers. But that comes out as n^2. So we need a different approach.
    
    // What we're going to do is read left to right, keeping a list of credit
    // providers. New credit providers not consistent with anything on the list
    // are added; if they are consistent with something they replace that thing.
    // As we read along, we count up mismatches against the positions implied by
    // each credit provider, and drop those that get too many mismatches. When
    // we hit an unmapped base, we give it credit from the left if we have
    // exactly one credit provider for it.
    
    // Then we do the same thing from the right, and make sure they agree.
    
    // We keep a table of all the credits which are assigned to bases. This
    // holds the Mapping (which may be unmapped) that credit wants to give each
    // base, or nothing if credit has nothing to say yet about that base.
    std::map<size_t, Mapping> credits;
    
    for(int direction = 1; direction >= -1; direction -= 2) {
        // Looking left to right first, then looking right to left.
    
        Log::info() << "Applying credit in direction " << direction <<
            std::endl;
    
        // We keep a table of how many mismatches each active credit provider
        // has incurred going this direction, by provider position. The keys of
        // this are our active credit providers.
        std::map<size_t, size_t> activeMismatches;
        
        for(size_t i = (direction > 0 ? 0 : disambiguated.size() - 1); 
            (direction > 0 ? i < disambiguated.size() : i != (size_t) -1);
            i += direction) {
            // For each base, reading in the appropriate direction...
            
            for(auto iterator = activeMismatches.begin(); 
                iterator != activeMismatches.end(); 
                iterator = (iterator->second > z_max ? 
                    // This credit provider has too many mismatches. Delete it
                    // and move to the next one.
                    activeMismatches.erase(iterator)
                : 
                    // Else keep it and move to the next one.
                    ++iterator)) {
                
                // The body runs before the update logic, so we can set the
                // mismatch count here and have it acted upon there.
                
                // Where is the credit provider?
                size_t provider = iterator->first;
                
                // Where did it map?
                TextPosition providerPos = 
                    disambiguated[provider].getLocation();
                
                // Grab a reference to the reference string that the provider
                // wants to put us on. Make sure to convert from text to contig.
                const std::string& referenceText = index.displayContigCached(
                    providerPos.getText() / 2);
                    
                // Advance the provider mapping position to get our implied
                // mapping position.
                providerPos.addOffset(i - provider);
                
                // Work out how many mismatches exist at exactly this one
                // implied position. Should be 0 or 1. TODO: See if we still
                // need this helper function and if it can be refactored.
                size_t mismatches = countMismatches(referenceText, 
                    providerPos.getStrand(), providerPos.getOffset(), query, i,
                    1);
            
                // Add the possible mismatches into the total mismatch strikes
                // against this credit provider.
                iterator->second += mismatches;
                
                // If it gets too many mismatches, the loop update logic will
                // remove it when advancing to the next credit provider.
                
                if(iterator->second > z_max) {
                    // We're going to drop this, so talk about it.
                    Log::debug() << "Credit provider " << iterator->first <<
                        " drops out at position " << i << std::endl;
                }
            }
            
            // TODO: Apply credit to this base from the left, if applicable.
            
            if(isCreditProvider[i]) {
                // We found a credit provider
                
                // Pull out its TextPosition
                TextPosition mapped = disambiguated[i].getLocation();
                
                // We haven't found any credit providers consistent with this
                // one.
                bool found = false;
                
                for(auto& keyValue : activeMismatches) {
                    // Check against all the other current credit providers
                    
                    // Pull out the old TextPosition
                    TextPosition existing =
                        disambiguated[keyValue.first].getLocation();
                    
                    if(existing.isConsistent(mapped, i - keyValue.first)) {
                        // This new mapping is consistent with the old mapping.
                        found = true;
                        
                        if(keyValue.second > 0) {
                            // The old mapping has accumulated some mismatches
                            // that we can skip if we give credit from the new
                            // one instead.
                            
                            Log::debug() << "Credit provider at " << i << 
                                " supplants that at " << keyValue.first <<
                                std::endl;
                            
                            // Remove the old entry
                            activeMismatches.erase(keyValue.first);
                            
                            // Add the new entry with 0 mismatches.
                            activeMismatches[i] = 0;
                            
                            // Leave the loop since we modified the map.
                            break;
                        } else {
                            Log::trace() << "Credit provider at " << i << 
                                " subsumed by that at " << keyValue.first <<
                                std::endl;
                        }
                        // Else the old provider is no different than the new one in
                        // terms of how it gives credit, so we don't need to mess
                        // with the map.
                        
                    }
                }
                
                
                if(!found) {
                    // There's nothing in the map currently to take care of this
                    // credit. Add the new entry with 0 mismatches.
                    // TODO: Can I refactor to do this once only?
                    activeMismatches[i] = 0;
                    
                    Log::debug() << "New credit provider at position " << i <<
                        std::endl;
                }
            } else if(i > leftSentinel && i < rightSentinel &&
                !leftMappings[i].isMapped() && !rightMappings[i].isMapped()) {
                // This base is between sentinels, and unmapped due to having no
                // mapping on either side (i.e. is unmapped and not
                // conflictingly mapped). We can apply credit to it.
                
                Log::trace() << "Position " << i <<
                    " eligible for credit from " << activeMismatches.size() <<
                    " providers" << std::endl; 
                
                if(activeMismatches.size() > 0) {
                    // What mapping are we going to make?
                    Mapping mapping;
                        
                    if(activeMismatches.size() == 1) {
                        // There is exactly one outstanding source of credit.
                        
                        // Grab its position
                        size_t provider = activeMismatches.begin()->first;
                        
                        // Where did it map?
                        TextPosition implied =
                            disambiguated[provider].getLocation();
                        
                        // TODO: do the check against giving credit to things
                        // with the wrong bases here?
                            
                        // Advance the provider mapping position to get our
                        // implied mapping position.
                        implied.addOffset(i - provider);
                        
                        Log::debug() << "Credit in direction " << direction <<
                            " maps " << i << " to " << implied << std::endl;
                            
                        // Make a mapped Mapping
                        mapping = Mapping(implied);
                        
                        
                    } else {
                        // There are multiple conflicting credits for this base
                        // from this side. TODO: Check to see if some of them
                        // can't apply to us because they imply a merge of two
                        // different bases.
                        
                        // Leave the unmapped Mapping to be added.
                        
                        Log::debug() << "Credit in direction " << direction <<
                            " conflicts at " << i << std::endl;
                    }
                    
                    // Now mapping is filled in with whatever mapping from
                    // credit we give (mapped somewhere or conflicted).
                    
                    if(credits.count(i)) {
                        // This base was mapped on credit in the other
                        // direction.
                        
                        if(credits[i] != mapping) {
                            // Change the mapping already there to be unmapped
                            // if it doesn't match this one.
                            credits[i] = Mapping();
                        }
                    } else {
                        // This base needs a new mapping, since we had something
                        // to say about it.
                        credits[i] = mapping;
                    }
                    
                }
            }
        }
    }
    
    // This is going to hold our output
    std::vector<Mapping> toReturn;
    
    for(size_t i = 0; i <= leftSentinel; i++) {
        // For each base before or at the left sentinel, disambiguate normally.
        toReturn.push_back(disambiguated[i]);
    }
    
    // How many bases successfully map on credit?
    size_t mappedOnCredit = 0;
    // How many bases got conflicting credit?
    size_t conflictedCredit = 0;

    for(size_t i = leftSentinel + 1; i < rightSentinel; i++) {    
        // For each base between the sentinels
    
        if(disambiguated[i].isMapped()) {
            // If it is mapped on one or more sides, disambiguate normally.
            toReturn.push_back(disambiguated[i]);   
        } else if(credits.count(i) > 0) {
            // Some credit was offered to it. Take it.
            toReturn.push_back(credits[i]);
            
            if(credits[i].isMapped()) {
                // Credit was a real mapping
                mappedOnCredit++;
            } else {
                // Credit was conflicted.
                conflictedCredit++;
            }
            
        } else {
            // No way to map it, say it is unmapped.
            toReturn.push_back(Mapping());
        }
        
    }
    
    // Dump some info about credit mappings
    Log::info() << mappedOnCredit << 
        " bases mapped on credit (before identity check)" << std::endl;
    Log::info() << conflictedCredit << 
        " bases did not map due to conflicting credit" << std::endl;
    
    for(size_t i = rightSentinel; i < disambiguated.size(); i++) {
        // Then for each base at or after the right sentinel, disambiguate
        // normally.
        toReturn.push_back(disambiguated[i]);
    }
    
    // Return the credit-mapped stuff
    return std::move(toReturn);
    
}
