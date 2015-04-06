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
     * If true, abort a mapping if you come to a retraction where it is too hard
     * to tell if the left and right have any overlaps. If false, just stop
     * searching that DP branch but continue all the other branches.
     */
    bool giveUpIfHard = false;
    
    /**
     * What is the minimum total context length to accept?
     */
    size_t minContextLength = 0;
    
    /**
     * What is the maximum number of BWT merged ranges we will check for a place
     * where left and right contexts will agree on, at each retraction step.
     */
    size_t maxRangeCount = (size_t) 100;
    
    /**
     * What is the maximum number of bases that we are willing to extend through
     * when trying to confirm if a context unique on one side is consistent with
     * the other?
     */
    size_t maxExtendThrough = (size_t) 100;
    
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
         * How many bases were searched on the left before the last retraction?
         */
        size_t lastLeftContext;
        /**
         * How many bases were searched on the right before the last retraction?
         */
        size_t lastRightContext;
        /**
         * Each DPTask keeps a reference back to the shared map linking left
         * context length to the set of positions selected on the left.
         */
        std::map<size_t, std::set<TextPosition>>& leftContextToPositions;
        /**
         * Each DPTask keeps a reference back to the shared map linking right
         * context length to the set of positions selected on the right.
         */
        std::map<size_t, std::set<TextPosition>>& rightContextToPositions;
        /**
         * Has this DPTask only ever been retracted right since the root?
         */
        bool isRightEdge;
        /**
         * Is the DP task still usable for set comparison, or are we only using
         * the extend-through heuristic because one of the sets on one side
         * would have gotten too big?
         */
        bool setsValid;
        
        /**
         * This is the correct way to read the set of TextPositions that this
         * DPTask's parent had selected on the left.
         * May never be called on a root task.
         */
        inline const std::set<TextPosition>& getLastLeftPositions() const {
            return leftContextToPositions.at(lastLeftContext);
        }
        
        /**
         * This is the correct way to read the set of TextPositions that this
         * DPTask's parent had selected on the right.
         * May never be called on a root task.
         */
        inline const std::set<TextPosition>& getLastRightPositions() const {
            return rightContextToPositions.at(lastRightContext);
        }
        
        /**
         * Set the set of TextPositions which are selected by this DPTask on the
         * left.
         */
        inline void setLeftPositions(std::set<TextPosition>&& leftPositions) {
            
            leftContextToPositions[leftContext] = std::move(leftPositions);
        }
        
        /**
         * Set the set of TextPositions which are selected by this DPTask on the
         * right.
         */
        inline void setRightPositions(std::set<TextPosition>&& rightPositions) {
            
            rightContextToPositions[rightContext] = std::move(rightPositions);
        }
        
        /**
         * Make a new DPTask for the root of the DP tree (with nothing yet
         * retracted).
         */
        inline DPTask(const FMDPosition& initialLeft, size_t initialLeftContext,
            std::map<size_t, std::set<TextPosition>>& leftContextToPositions, 
            const FMDPosition& initialRight, size_t initialRightContext,
            std::map<size_t, std::set<TextPosition>>& rightContextToPositions): 
            left(initialLeft), right(initialRight), lastLeft(left), 
            lastRight(right), leftContextToPositions(leftContextToPositions), 
            rightContextToPositions(rightContextToPositions),
            leftContext(initialLeftContext), rightContext(initialRightContext),
            lastLeftContext((size_t) -1), lastRightContext((size_t) -1),
            isRightEdge(true), setsValid(true) {
                
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
                retracted.lastRightContext = retracted.rightContext;
                // Retract the right context.
                retracted.rightContext = index.retractRightOnly(
                    retracted.right);
            } else {
                // Save history
                retracted.lastLeft = retracted.left;
                retracted.lastLeftContext = retracted.leftContext;
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
