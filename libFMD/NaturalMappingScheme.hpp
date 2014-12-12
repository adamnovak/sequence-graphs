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
     * How many mismatches can be tolerated when applying credit?
     */
    size_t z_max = 0;
    
    /**
     * Should credit be used to create more mappings?
     */
    bool credit = false;
    
    /**
     * Should synteny blocking be used to bridge mismatches?
     */
    bool synteny = false;
    
    /**
     * Should costs for breaks between maximal matches in a synteny block be a
     * flat 1, as opposed to scaling with gap length?
     */
    bool flatCost = false;
    
    /**
     * What's the minimum context length to match on?
     */
    size_t minContext = 0;
    
    /**
     * What is the minimum length that a maximal unique match has to be to even
     * contribute to conflicts? Maximal unique matches shorter than this will be
     * completely ignored, and won't be able to stop any bases from mapping on
     * other maximal unique matches.
     */
    size_t ignoreMatchesBelow = 0;
    
    /**
     * How many non-overlapping minimal unique matches does a maximal unique
     * match need in order to map?
     *
     * This is a lower bound on the maximal unique match's Hamming distance from
     * every other place in the reference except the place it mapped to.
     */
    size_t minHammingBound = 0;
    
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
    
        // TODO: No forgery?
        /**
         * What is the first (leftmost) base in the query in this synteny block?
         */
        size_t start = 0;
        
        /**
         * How many many bases are in the block?
         */
        size_t length = 0;
        
        /**
         * How many minimal exact matchings are contained by the maximal exact
         * matchings of this synteny block?
         */
        size_t minimalExactMatchings = 0;
        
        /**
         * How many *non-overlapping* minimal exact matchings are contained by
         * the maximal exact matchings of this SyntenyBlock?
         */
        size_t nonOverlapping = 0;
        
        /**
         * How many bases of minimal exact matchings are contained in the
         * maximal exact matchings of this block? We need this for calculating
         * their average length.
         */
        size_t totalMinimalExactMatchingBases = 0;
        
        // How many non-overlapping minimal unique matchings were contained in
        // the last maximal unique matching that was added? This is used to
        // compute if it is worth connecting a new maximal unique matching to
        // the block.
        size_t lastNonOverlapping = 0;
        
        /**
         * How many mismatches are crossed by this SyntenyBlock?
         */
        size_t mismatches = 0;
        
        // Default constructor is OK, we always start empty.
        
        // Default copy constructor, assignment operator, and so on are fine.
        
        /**
         * Extend a SyntenyBlock with the next maximal unique Matching to the
         * left. Takes the new Matching, the number of minimal exact matchings
         * it contains, the number of total bases in those matchings, and the
         * numebr of those matchings that are non-overlapping. Also takes the
         * cost to incur, which only the caller knows.
         */
        inline SyntenyBlock extendLeft(Matching maximal,
            size_t newMinimalExatMatchings,
            size_t newMinimalExactMatchingBases,
            size_t newNonOverlapping, size_t cost) const {
            
            // Make a new SyntenyBlock to extend.
            SyntenyBlock toReturn = *this;
            
            // Incur the mismatch cost
            toReturn.mismatches += cost;
            
            if(toReturn.maximalMatchings.empty()) {
                // We were empty, so grab only the bases from this new Matching.
                toReturn.length = maximal.length;
            } else {
                // Grab all the bases through to the start of this new leftmost
                // Matching.
                toReturn.length += toReturn.start - maximal.start;
            }
            
            // Set the start position
            toReturn.start = maximal.start;
            
            // Update all the statistics
            toReturn.minimalExactMatchings += newMinimalExatMatchings;
            toReturn.totalMinimalExactMatchingBases +=
                newMinimalExactMatchingBases;
            toReturn.nonOverlapping += newNonOverlapping;
            
            // Remember the nonOverlapping score of the last added matching.
            toReturn.lastNonOverlapping = newNonOverlapping;
            
            // And add in the matching
            toReturn.maximalMatchings.push_back(maximal);
            
            // Return it.
            return toReturn;
        }
        
        /**
         * Is this next maximal Matching to the left consistent with this
         * synteny block, and therefore elligible to be a member of it? It is
         * not consistent if it overlaps into our block, or if when we get the
         * locations the offsets are wrong.
         */
        inline bool isConsistent(Matching maximal) const {
            if(maximalMatchings.empty()) {
                // Anything is consistent with an empty SyntenyBlock.
                return true;
            }
        
            // How many bases from the start of this new thing to the leftmost
            // thing we are holding?
            size_t observedOffset = maximalMatchings.rbegin()->start -
                maximal.start;
            
            // Get the position at which the new matching maps    
            TextPosition position = maximal.location;
            // And budge it forward by that same offset            
            position.addLocalOffset(observedOffset);
            
            // Make sure it ends up being the corresponding location.
            return position == maximalMatchings.rbegin()->location;
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
     * Given a vector of SyntenyBlocks, decide if each one passes the filters or
     * not.
     */
    std::vector<bool> filter(const std::vector<SyntenyBlock>& blocks) const;
    
    /**
     * What should the cost be for attacking this Matching to this SyntenyBlock?
     * Part of the MappingScheme since it needs scheme parameters.
     */
    inline size_t cost(Matching maximal, const SyntenyBlock& block) const {
        if(block.maximalMatchings.empty()) {
            // It costs nothing to connect if we're empty.
            return 0;
        }
    
        if(flatCost) {
            // Make each mismatch block count for 1.
            return 1;
        } else {
            // Use the length as the cost.
            return block.start - (maximal.start + maximal.length);
        }
        // TODO: Try scanning and counting up actual mismatches.
    }
    
    /**
     * Given that it's possible to connect this new maximal unique matching
     * containing the given number of non-overlapping minimal unique matchings
     * to this synteny block, is it worth it? Part of the MappingScheme since it
     * needs scheme parameters.
     */
    inline bool worthConnecting(Matching maximal,
        size_t newNonOverlapping, const SyntenyBlock& block) const {
        
        // If we are empty, it's always worth it.
        if(block.maximalMatchings.empty()) {
            return true;
        }
        
        // A connection to the most recent maximal matching we have is worth
        // it (in both directions) if the cost of the connection (in
        // mismatches crossed) is less than or equal to the number of
        // minimal unique matchings contained in each maximal unique
        // matching on its own.
        size_t connectionCost = cost(maximal, block);
        return connectionCost <= newNonOverlapping &&
            connectionCost <= block.lastNonOverlapping;
    }
    
};
   
#endif
