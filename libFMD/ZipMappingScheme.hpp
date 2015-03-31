#ifndef ZIPMAPPINGSCHEME_HPP
#define ZIPMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"
#include "Mapping.hpp"
#include "IntervalIndex.hpp"
#include "Log.hpp"
#include "Matching.hpp"

#include <iomanip>
#include <queue>
#include <unordered_map>

/**
 * Mapping scheme supporting mapping to graphs, where you zip together maximal
 * matches looking out in both directions from each position. Finds a subset of
 * the unique substrings for each base that the natural mapping scheme finds,
 * but will never map a base that the actual natural to-graph mapping scheme
 * would call as conflicted.
 */
class ZipMappingScheme: public MappingScheme {

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
        
        
    // Mapping scheme parameters
    
    /**
     * Determines whether the retraction dynamic programming code will be used
     * to try and find shorter left and right contexts that are unique taken
     * together.
     */
    bool useRetraction = true;
    
    /**
     * What is the minimum total context length to accept?
     */
    size_t minContextLength = 20;
    
    /**
     * What is the maximum number of BWT merged ranges we will check for a place
     * where left and right contexts will agree on, at each retraction step.
     */
    size_t maxRangeCount = 10;
    
    /**
     * What is the maximum number of bases that we are willing to extend through
     * when trying to confirm if a context unique on one side is consistent with
     * the other?
     * TODO: Make this work for more than just the very top level.
     */
    size_t maxExtendThrough = 20;
    
protected:

    /**
     * Represents a retraction that needs to be checked for results.
     */
    struct DPTask {
        /**
         * Left search to look at.
         */
        FMDPosition left;
        /** 
         * Right search to look at.
         */
        FMDPosition right;
        /**
         * Left search interval from the last step. If not equal to left, then
         * this task represents a retraction on the left.
         */
        FMDPosition lastLeft;
        /**
         * Right search interval from the last step. If not equal to right, then
         * this task represents a retraction on the right.
         */
        FMDPosition lastRight;
        /**
         * How many bases are searched on the left?
         */
        size_t leftContext;
        /**
         * How many bases are searched on the right?
         */
        size_t rightContext;
        /**
         * What TextPositions are known to be selected by lastLeft? Ignored for
         * the root. These are going to have their strands backwards to the
         * actual positions that correspond to the searched *left* context.
         */         
        std::set<TextPosition> lastLeftPositions;
        /**
         * What TextPositions are known to be selected by lastRight? Ignored for
         * the root.
         */  
        std::set<TextPosition> lastRightPositions;
        /**
         * Has this DPTask only ever been retracted right since the root?
         */
        bool isRightEdge;
        
        /**
         * Make a new DPTask for the root of the DP tree (with nothing yet
         * retracted).
         */
        inline DPTask(const FMDPosition& initialLeft, size_t initialLeftContext,
            const FMDPosition& initialRight, size_t initialRightContext): 
            left(initialLeft), right(initialRight), lastLeft(left), 
            lastRight(right), leftContext(initialLeftContext), 
            rightContext(initialRightContext), isRightEdge(true) {
                
            // lastLeftPositions and lastRightPositions will be ignored since
            // this is the root.
        }
        
        /**
         * Return a copy retracted on the left (if false) or the right (if
         * true). It's the caller's responsibility to fix up lastLeftPositions
         * and/or lastRightPositions as appropriate, depending on which side was
         * retracted on to get the DPTask now being retracted.
         */
        inline DPTask retract(const FMDIndex& index, bool isRight) const {
            // Make the copy
            DPTask retracted = *this;
            
            if(isRight) {
                // Save history
                retracted.lastRight = retracted.right;
                // Retract the right context.
                retracted.rightContext = index.retractRightOnly(
                    retracted.right);
            } else {
                // Save history
                retracted.lastLeft = retracted.left;
                // Retract the left context (still on its local right).
                retracted.leftContext = index.retractRightOnly(
                    retracted.left);
                // This retracted copy has now been retracted on the left, and
                // no longer needs to throw off right retractions.
                retracted.isRightEdge = false;
            }
            
            // Give back the modified copy
            return retracted;
        }
    }; 

    /**
     * Use the inchworm algorithm to find the longest right context present in
     * the reference for each base in the query. Results are in the same order
     * as the characters in the string, and consist of an FMDPosition of search
     * results and a context length.
     */
    std::vector<std::pair<FMDPosition, size_t>> findRightContexts(
        const std::string& query) const;
        
    /**
     * Returns true if the given unique search result can be extended through
     * the given context string on the other side. That string must end with
     * the base shared with the search result set (which is not used to extend
     * again), and is extended through right to left.
     */
    bool canExtendThrough(FMDPosition context,
        const std::string& opposingQuery) const;
    
    /**
     * Evaluate the retraction represented by the given DPTask (which may be
     * either a root DPTask or a retraction). If it is possible to ascertain
     * whether the DPTask contains 0, 1, or many shared results, it will return
     * true and a set of 0, 1, or many results. If it is not possible to make
     * that determination (for example, because the number of newly selected
     * positions on the retracted side is too large) it will return false and an
     * unspecified set (in which case the position being examined needs to be
     * banned from mapping.
     *
     * Is not responsible for deciding whether a finding of results by one
     * DPTask means that other DPTasks do not need to be executed.
     */
    std::pair<bool, std::set<TextPosition>> exploreRetraction(
        const DPTask& task, std::queue<DPTask>& taskQueue,
        const std::string& query, size_t queryBase) const;
    
    /**
     * Explore all retractions of the two FMDPositions. Return either an empty
     * set (if no explored retraction finds overlapping positions between the
     * two sides), a set with one element (if we find exactly one such
     * overlap), or a set with two (or more) elements (if we find multiple
     * overlaps).
     */
    std::set<TextPosition> exploreRetractions(const FMDPosition& left,
        size_t patternLengthLeft, const FMDPosition& right,
        size_t patternLengthRight, const std::string& query,
        size_t queryBase) const;    
    
};

#endif
