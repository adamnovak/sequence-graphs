#ifndef MISMATCHRESULTSET_HPP
#define MISMATCHRESULTSET_HPP

#include <set>
#include "GenericBitVector.hpp"
#include "MismatchResult.hpp"

// We need to know FMDIndex exists
class FMDIndex;


/**
 * Represents a set of MismatchResult objects being used for some single search.
 */
class MismatchResultSet {
public:
    
    // Default copy constructor is OK, but assignment won't work because we
    // contain a reference.
    
    /**
     * Make a new MismatchResultSet that searches on the given index.
     */
    MismatchResultSet(const FMDIndex& index);
    
    /**
     * And the default copy constructor.
     */
    MismatchResultSet(const MismatchResultSet& other) = default;
    
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
    void filter(size_t z_max, const GenericBitVector* mask = NULL);
    
    /**
     * Return the one range number in the range vector overlapped by results, or
     * -1 if no such unique range exists. Only count positions masked on in the
     * mask, if present.
     */
    int64_t range(const GenericBitVector& ranges, 
        const GenericBitVector* mask = NULL) const;
        
    /**
     * Return if the result set is empty of actual BWT positions under the given
     * mask.
     */
    bool isEmpty(const GenericBitVector* mask = NULL) const;
    
    /**
     * Return the number of positions actually found, under the given mask.
     */
    int64_t getLength(const GenericBitVector* mask = NULL) const;
    
    /**
     * Return true if we have found exactly one range, false otherwise.
     */
    bool isMapped(const GenericBitVector& ranges, 
        const GenericBitVector* mask = NULL) const;
    
protected:
    // We need to keep track of the index we use
    const FMDIndex& index;
    // We need to keep our MismatchResults
    std::set<MismatchResult> results;
};

#endif
