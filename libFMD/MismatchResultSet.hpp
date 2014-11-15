#ifndef MISMATCHRESULTSET_HPP
#define MISMATCHRESULTSET_HPP

#include <set>
#include "GenericBitVector.hpp"
#include "MismatchResult.hpp"


/**
 * Represents a set of MismatchResult objects being used for some single search.
 */
class MismatchResultSet {
public:
    
    // Default copy constructor, assignment operator, and so on are OK.
    
    /**
     * Make a new MismatchResultSet that searches on the given index.
     */
    MismatchResultSet(const FMDIndex& index);
    
    /**
     * Extend all the results on the left with matches and mismatches for the
     * given base.
     */
    void extendLeft(char base);
    
    /**
     * Retract all the results on the right by 1 base.
     */
    void retractRight();
    

    /**
     * Remove all results which have more than z_max mismatches, or which have
     * no BWT positions selected under the given mask.
     */
    void filter(size_t z_max, GenericBitVector* mask = NULL);
    
    /**
     * Return the one range number in the range vector overlapped by results, or
     * -1 if no such unique range exists. Only count positions masked on in the
     * mask, if present.
     */
    int64_t range(const GenericBitVector& ranges, 
        GenericBitVector* mask = NULL);
    
protected:
    // We need to keep track of the index we use
    const FMDIndex& index;
    // We need to keep our MismatchResults
    std::set<MismatchResult> results;
};

#endif
