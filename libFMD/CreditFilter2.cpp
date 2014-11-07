#include "CreditFilter2.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>
#include "Log.hpp"

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
    const GenericBitVector& ranges, size_t z_max): index(index), ranges(ranges),
    z_max(z_max), disambiguate(index) {
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
            
            // Get the length of its word, which has its right end at i
            leftWordLength = disambiguated[i].getLeftMinContext();
            
            // Clip out the word that the base mapped on on the left.
            std::string word = query.substr(i + 1 - leftWordLength,
                leftWordLength);
            
            if(index.misMatchCount(ranges, word, z_max).is_mapped) {
                // This word appears exactly once within the specified number of
                // mismatches. So this is the leftmost left sentinel.
                leftSentinel = i;
                Log::info() << "Left sentinel found at " << leftSentinel << 
                    std::endl;
                break;
            }
        }
    }
    
    for(size_t i = disambiguated.size() - 1; i != (size_t) -1; i--) {
        // Look in from the right for the right sentinel
        if(rightMappings[i].isMapped() && 
            disambiguated[i].isMapped()) {
            
            // This is a right-mapped thing on the right, and may be a sentinel.
            // To be sure, we have to find its word.
            
            // Get the length of its word, which has its left end at i
            rightWordLength = disambiguated[i].getRightMinContext();
            
            // Clip out the word that the base mapped on on the right.
            std::string word = query.substr(i, rightWordLength);
            
            if(index.misMatchCount(ranges, word, z_max).is_mapped) {
                // This word appears exactly once within the specified number of
                // mismatches. So this is the rightmost right sentinel.
                rightSentinel = i;
                Log::info() << "Right sentinel found at " << rightSentinel << 
                    std::endl;
                break;
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
        // the left sentinel) covers the left sentinel. So if the right max
        // context reaches back to exactly the left end of the left word, these
        // will be equal.
        return disambiguated[i].getLeftMaxContext() >= 
            i - leftSentinel + leftWordLength;
    };
    
    // And another to see if a base;'s maximal right context can extend over the
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
        
    }
    
    // The sentinels always have to come out as credit providers.
    assert(isCreditProvider[leftSentinel]);
    assert(isCreditProvider[rightSentinel]);
    
    // Now we can actually try applying credit.
    
    // Find the max left and right contexts we need to worry about checking
    // consistency from.
    size_t maxLeftContext = 0;
    size_t maxRightContext = 0;
    
    for(size_t i = 0; i < disambiguated.size(); i++) {
        // Take the max left and right context when we find them.
        
        Log::debug() << disambiguated[i].getLeftMaxContext() << ", " << 
            disambiguated[i].getRightMaxContext() << std::endl;
        
        maxLeftContext = std::max(maxLeftContext, 
            disambiguated[i].getLeftMaxContext());
        maxRightContext = std::max(maxRightContext, 
            disambiguated[i].getRightMaxContext());
    }
    
    Log::debug() << "Max context sizes: " << maxLeftContext << "|" << 
        maxRightContext << std::endl;
    
    // This is going to hold our output
    std::vector<Mapping> toReturn;
    
    for(size_t i = 0; i <= leftSentinel; i++) {
        // For each base before or at the left sentinel, disambiguate normally.
        toReturn.push_back(disambiguated[i]);
    }

    for(size_t i = leftSentinel + 1; i < rightSentinel; i++) {    
        // For each base between the sentinels
    
        if(disambiguated[i].isMapped()) {
            // If it is mapped on one or more sides, disambiguate normally.
            toReturn.push_back(disambiguated[i]);   
        } else {
        
            Log::trace() << "Trying to credit map base " << i << std::endl;
        
            // Set this to true if you find a base that implies a position for
            // this one by its right context.
            bool rightFound = false;
            // This is where that base would place this one
            TextPosition rightCreditPosition(0, 0);
            // Set this if right contexts all place this base in one spot
            bool rightConsistent = true;
        
            for(size_t j = i - 1; j != (size_t) -1 && 
                (int64_t) j >= (int64_t) i - (int64_t) maxRightContext; j--) {
                // Look left from here until you find a base that maps and
                // places us.
                Log::trace() << "Checking base " << j << " on right" << 
                    std::endl;
                
                if(!isCreditProvider[j]) {
                    // This base either didn't map or possibly could unmap, so
                    // don't take credit from it.
                    continue;
                }
                
                // Go get the text we would map on. This is always strand 0.
                // TODO: Proper weakref cache thingy.
                const std::string& referenceText = index.displayContigCached(
                    disambiguated[j].getLocation().getText());
                    
                // Work out what offset we should have on that strand.
                int64_t offset = (int64_t) i - (int64_t) j + 
                    (int64_t) disambiguated[j].getLocation().getOffset(); 
                    
                if(offset < 0 || offset >= referenceText.size()) {
                    // This base can't actually give credit to us since it would
                    // put us off the end of its strand.
                    continue;
                }
                
                // Make the mapping we imply
                TextPosition implied(disambiguated[j].getLocation().getText(),
                    offset);
                
                // OK the base can give credit, but can it give credit on its
                // particular text at this distance, even if it happens to be
                // unmapped on this side? Remember, we need credit from unmapped
                // sides to make sure credit is consistent with sentinels.
                // Otherwise we would have to special-case that credit was
                // consistent with sentinels.
                
                // We need to count the number of mismatches between the query
                // string and the reference, from the credit provider to here.
                // So in the reference we start at the credit provider's mapping
                // location (on the appropriate strand), and in the query we
                // start at the credit provider, and we go forward until we get
                // to this base.
                size_t mismatches = countMismatches(referenceText, 
                    disambiguated[j].getLocation().getStrand(), 
                    disambiguated[j].getLocation().getOffset(), query, j, 
                    i - j + 1);
                
                if(mismatches > z_max) {
                    // No credit because too many mismatches.
                    continue;
                }
                
                Log::trace() << "Base " << j << " places base " << i << 
                    " at " << implied << " by right context with " <<
                    mismatches << " mismatches" << std::endl;
                
                if(!rightFound) {
                    // This is the first one. Make sure all the others match.
                    rightFound = true;
                    rightCreditPosition = implied;
                } else if (rightCreditPosition != implied) {
                    // There are two or more implied locations from right
                    // contexts.
                    rightConsistent = false;
                    break;
                }
                
            }
            
            // And the same thing for left contexts
            
            // Set this to true if you find a base that implies a position for
            // this one by its left context.
            bool leftFound = false;
            // This is where that base would place this one
            TextPosition leftCreditPosition(0, 0);
            // Set this if left contexts all place this base in one spot
            bool leftConsistent = true;
        
            for(size_t j = i + 1; j < leftMappings.size() && 
                j < i + maxLeftContext; j++) {
                
                // Look right from here until you find a base that maps and
                // places us.
                
                Log::trace() << "Checking base " << j << " on left" << 
                    std::endl;
                
                if(!isCreditProvider[j]) {
                    // This base either didn't map or possibly could unmap, so
                    // don't take credit from it.
                    continue;
                }
                
                // Go get the text we would map on. This is always strand 0.
                // TODO: Proper weakref cache thingy.
                const std::string& referenceText = index.displayContigCached(
                    disambiguated[j].getLocation().getText());
                    
                // Work out what offset we should have on that strand.
                int64_t offset = (int64_t) i - (int64_t) j + 
                    (int64_t) disambiguated[j].getLocation().getOffset(); 
                    
                if(offset < 0 || offset >= referenceText.size()) {
                    // This base can't actually give credit to us since it would
                    // put us off the end of its strand.
                    continue;
                }
                
                // Make the mapping we imply
                TextPosition implied(disambiguated[j].getLocation().getText(),
                    offset);
                
                // OK the base can give credit, but can it give credit on its
                // particular text at this distance, even if it happens to be
                // unmapped on this side? Remember, we need credit from unmapped
                // sides to make sure credit is consistent with sentinels.
                // Otherwise we would have to special-case that credit was
                // consistent with sentinels.
                
                // We need to count the number of mismatches between the query
                // string and the reference, from here to the credit provider.
                // So in the reference we start at the implied mapping location
                // (on the appropriate strand), and in the query we start at us,
                // and we go forward until we get to the credit provider.
                size_t mismatches = countMismatches(referenceText, 
                    implied.getStrand(), implied.getOffset(), query, i, 
                    j - i + 1);
                    
                if(mismatches > z_max) {
                    // No credit because too many mismatches.
                    continue;
                }
                
                Log::trace() << "Base " << j << " places base " << i << 
                    " at " << implied << " by left context with " <<
                    mismatches << " mismatches" << std::endl;
                
                if(!leftFound) {
                    // This is the first one. Make sure all the others match.
                    leftFound = true;
                    leftCreditPosition = implied;
                } else if (leftCreditPosition != implied) {
                    // There are two or more implied locations from left
                    // contexts.
                    leftConsistent = false;
                    break;
                }
                
            }
            
            // TODO: somehow call into the disambiguation filter here?
        
            if(leftFound && leftConsistent) {
                if(rightFound && rightConsistent) {
                    // We have credit from both the left and the right. Do they
                    // agree?
                    
                    if(leftCreditPosition == rightCreditPosition) {
                        // They agree! Record a successful mapping. TODO: do we
                        // want to send credit info along here as well?
                        toReturn.push_back(Mapping(leftCreditPosition));
                    } else {
                        // They disagree. Fail the mapping
                        toReturn.push_back(Mapping());
                    }
                } else {
                    // We have credit only from the left.
                    toReturn.push_back(Mapping(leftCreditPosition));
                }
            } else if(rightFound && rightConsistent) {
                // We have credit only from the right
                
                toReturn.push_back(Mapping(rightCreditPosition));
            } else {
                // We have credit from nowhere
                toReturn.push_back(Mapping());
            }
        }
        
    }
    
    for(size_t i = rightSentinel; i < disambiguated.size(); i++) {
        // Then for each base at or after the right sentinel, disambiguate
        // normally.
        toReturn.push_back(disambiguated[i]);
    }
    
    // Return the credit-mapped stuff
    return std::move(toReturn);
    
}
