#ifndef NATURALMAPPINGSCHEME_HPP
#define NATURALMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"
#include "Mapping.hpp"
#include "Log.hpp"

/**
 * Mapping scheme implementing Benedict's "natural" mapping scheme. If all
 * the unique-in-the-reference strings overlapping a query base agree about
 * where the abse should be, and you can't get from the base off the end of
 * the query with multiple results, map the base.
 */
class NaturalMappingScheme: public MappingScheme {

public:
    /**
     * Inherit the constructor.
     */
    using MappingScheme::MappingScheme;
    
    /**
     * Map the given query string according to the natural mapping scheme with
     * optional inexact credit. When a mapping is found, the callback function
     * will be called with the query base index, and the TextPosition to which
     * it maps in the forward direction.
     *
     */
    virtual void map(const std::string& query,
        std::function<void(size_t, TextPosition)> callback) const override;
        
    
    // Now come the scheme parameters and their default values.
    
    /**
     * What's the minimum factor by which the maximum unique string length must
     * exceed the minimum?
     */
    double multContext = 0.0;
    
    /**
     * What's the minimum run length in bases to use?
     */
    size_t minContext = 0;
    
    /**
     * Should maximal unique matchings not in an admissible run still produce
     * conflict?
     */
    bool conflictBelowThreshold = false;
    
    /**
     * How many mismatches can be tolerated when applying credit?
     */
    size_t z_max = 0;
    
    /**
     * Should credit be used to create more mappings?
     */
    bool credit = false;
    
    /**
     * What is the minimum length that a maximal unique match has to be to even
     * contribute to conflicts? Maximal unique matches shorter than this will be
     * completely ignored, and won't be able to stop any bases from mapping on
     * other maximal unique matches.
     */
    size_t ignoreMatchesBelow = 0;
    
    /**
     * How many non-overlapping minimal unique matches does a maximal unique
     * match (or region of a synteny block) need in order to map?
     *
     * This is a lower bound on the maximal unique match's Hamming distance from
     * every other place in the reference except the place it mapped to.
     */
    size_t minHammingBound = 0;
    
    /**
     * How many mismatches are we allowed to have in gaps between maximal unique
     * matches that are working together to map?
     */
    size_t maxHammingDistance = 0;
    
    /**
     * How big of an alignment will we tolerate doing for edit distance
     * calculations?
     */
    size_t maxAlignmentSize = 10;
    
    /**
     * Should we enable potentially unstable (as in, mappings can change on
     * extension, not as in may crash the program) mapping hacks to increase
     * coverage when mapping reads?
     */
    bool unstable = false;
    
protected:
    /**
     * Map the given query string, producing a vector of Mappings. Does not
     * include credit yet.
     */
    std::vector<Mapping> naturalMap(const std::string& query) const;
    
    /**
     * Keep track of an occurrence of a unique string which matches a query
     * position to a reference position.
     */
    struct Matching {
        /**
         * Where does it start in the query?
         */
        size_t start;
        /**
         * Where is it in the reference?
         */
        TextPosition location;
        /**
         * How long is it?
         */
        size_t length;
        
        /**
         * Make a new Matching.
         */
        inline Matching(size_t start, TextPosition location, size_t length): 
            start(start), location(location), length(length) {
        }
        
        // We also have some metadata we use for maximal Matchings.
        // TODO: Make that its own type?
        
        /**
         * How many minimal matchings are in this maximal matching?
         */
        size_t minMatchings = 0;
        
        /**
         * What is their total length?
         */
        size_t minMatchingLength = 0;
        
        /**
         * How many of them are non-overlapping?
         */
        size_t nonOverlapping = 0;
        
        /**
         * How many mismatches do you have to collect before grabbing this
         * maximal unique matching in its synteny block? Assumes you are coming
         * in from the right.
         */
        size_t mismatchesBefore = 0;
        
        /**
         * Is this maximal match flagged as able to produce matchings?
         */
        bool canMatch = false;
        
        /**
         * Is this maximal match supposed to create conflicting matchings but
         * not map? If this is true, canMatch must also be true.
         */
        bool blacklist = false;
         
    };
    
    /**
     * A Run of Matchings (referenced by indices) with a certain total min
     * unique matchings and total Hamming/edit distance.
     */
    struct Run {
        /**
         * Holds indices of the Matchings in this run.
         */
        std::deque<size_t> matchings;
        
        /**
         * One shorter than the matchings vector. Holds the Hamming/edit distance
         * cost for the gap after the corresponding entry in the matchings
         * vector.
         */
        std::deque<size_t> costs;
        
        /**
         * Holds the total number of non-overlapping minimal unique matchings in
         * the run. TODO: Need to account for overlap in a run when maximal
         * unique matches in the same run overlap.
         */
        size_t totalClearance;
        
        /**
         * Total Hamming/edit distance between the query ant the reference for
         * this run.
         */
        size_t totalCost;
        
        /**
         * Equality comparison so we can put these in a set.
         */
        inline bool operator==(const Run& other) const {
            return(matchings == other.matchings && costs == other.costs &&
                totalClearance == other.tolalClearance &&
                totalCost == other.totalCost);
        }
        
        /**
         * Ordering defined so we can put these in a set.
         */
        inline bool operator<(const Run& other) const {
            return matchings < other.matchings ||
                matchings == other.matchings && 
                (costs < other.costs || 
                costs == other.costs && 
                (totalClearance < other.totalClearance || 
                totalClearance == other.totalClearance && 
                totalCost == other.totalCost));
        }
        
    };
    
    /**
     * Keep track of a run of maximal exact matchings which we have decided form
     * a synteny block, and ought to be presented to the filters as a group.
     * This keeps track of all of the state that is needed to implement the
     * various filters, as well as the information needed to actually make the
     * mappings.
     */
    struct SyntenyBlock {
        /**
         * What maximal matchings are in this block (in right to left order)?
         */
        std::vector<Matching> maximalMatchings;
    
        // Default constructor is OK, we always start empty.
        
        // Default copy constructor, assignment operator, and so on are fine.
        
        /**
         * Extend a SyntenyBlock with the next maximal unique Matching to the
         * left. Takes the new Matching, which must have its metadata filled in.
         * Modifies the block in place.
         */
        inline void extendLeft(Matching maximal) {
            // And add in the matching
            maximalMatchings.push_back(maximal);
        }
        
    };
    
    /**
     * Find all of the maximal unique matchings between query string characters
     * and the reference, in descending order by left endpoint.
     */
    std::vector<Matching> findMaxMatchings(const std::string& query) const;
        
    /**
     * Find all of the minimal unique matchings between query string characters
     * and the reference, in descending order by left endpoint.
     */
    std::vector<Matching> findMinMatchings(const std::string& query) const;
    
    /**
     * Scan through the given SyntenyBlock and do the inchworm algorithm to flag
     * all of the maximal matchings that pass the scheme's configured
     * thresholds. Also make sure to blacklist any maximal unique matches that
     * could eventually pass those thresholds on extension.
     */
    void scan(SyntenyBlock& block, const std::string& query) const;
    
    /**
     * Count the mismatches between a query and a text in the index, in a range.
     * If the range would go out of the query, stop counting there. If the range
     * would go out of the reference, return the threshold number of mismatches
     * (if set), or the total counted. If threshold is not -1, will count until
     * the end or until that many mismatches are found. Direction can be 1 or -1
     * and determines the direction to count in (forwards or backwards from the
     * inclusive start positions given).
     */
    size_t countMismatches(const std::string& query, size_t queryStart,
        TextPosition referenceStart, size_t length,
        int64_t threshold = -1, char direction = 1) const;
        
    /**
     * Compute the edit distance between the specified region of the query and
     * the specified region of the reference. If the distance would be greater
     * than the threshold, return th threshold value instead. Also returns the
     * threshold value if the alignment that would need to be done is larger
     * than maxAlignmentSize.
     */
    size_t countEdits(const std::string& query, size_t queryStart,
        size_t queryLength, TextPosition referenceStart, size_t referenceLength,
        int64_t threshold = -1) const;
    
    
};
   
#endif
