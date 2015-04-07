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
    
    /**
     * How many unique non-overlapping one-sided MUSes need to be in the one-
     * sided MUMs of a base for that base to map? TODO: these sould be non-
     * overlapping, and restricted to the range of bases used to map the base in
     * question.
     */
    size_t minUniqueStrings = 0;
    
protected:
    
    /**
     * Represents an entire DP job for evaluating all the retractions of a left
     * and a right context.
     *
     * Operates in retraction space, where retracting to the next point at which
     * you would get more BWT positions selected on a side is a retraction.
     */
    struct DPTable {
        // Holds all the stuff we need to know about a retraction on the left,
        // or on the right.
        struct SideRetractionEntry {
            /**
             * What is actually selected still?
             */
            FMDPosition position;
            /**
             * How many bases are still searched.
             */
            size_t contextLength;
            /**
             * Are the sets below actually filled in? Or was the number of
             * ranges we would have had to have looked at to fill them too
             * great?
             */
            bool setsValid = true;
            /**
             * What positions are selected in total by the current FMDPosition?
             */
            std::set<TextPosition> selected;
            /**
             * Of the above, which are from newly selected ranges?
             */
            std::set<TextPosition> newlySelected;
            
            /**
             * Default constructor that leaves all the sets empty.
             */
            inline SideRetractionEntry() {
                // Nothing to do!
            }
            
            /**
             * Start out with an entry reflecting no retraction at all.
             * Fills in the sets if possible.
             */
            inline SideRetractionEntry(const FMDPosition& unretracted,
                size_t contextLength, const FMDIndexView& view,
                size_t maxRangeCount): position(unretracted), 
                contextLength(contextLength) {
                
                if(view.getApproximateNumberOfRanges(position) <=
                    maxRangeCount) {
                
                    // We can visit everything we have selected
                    selected = view.getTextPositions(position);
                    // Copy it all to the newly selected set too.
                    newlySelected = selected;
                    // And say our sets we just made are valid.
                    setsValid = true;
                } else {
                    // We can't afford to fill in our sets.
                    setsValid = false;
                }
            }
            
            /**
             * Create an entry from this one after retracting to the next place
             * where more results are to be had under the given view. TODO:
             * assumes such a place actually exists.
             *
             * Gives up on sets if it would need to visit more than
             * maxRangeCount ranges.
             */
            inline SideRetractionEntry retract(const FMDIndexView& view,
                size_t maxRangeCount) {
                
                Log::debug() << "Retracting a SideRetractionEntry" << std::endl;
                
                // Copy and retract our current FMDPosition
                FMDPosition retracted = position;
                size_t retractedContext = view.getIndex().retractRightOnly(
                    retracted);
                
                // Make a new entry
                SideRetractionEntry toReturn;
                toReturn.position = retracted;
                toReturn.contextLength = retractedContext;
                
                toReturn.setsValid = setsValid;
                if(setsValid) {
                    // We have sets we can build on, but we need to see if we've
                    // selected too much more stuff.
                    
                    size_t newRanges = view.getApproximateNumberOfNewRanges(
                        position, retracted);
                        
                    Log::debug() << "Will have " << newRanges <<
                        " new ranges" << std::endl;
                        
                    if(newRanges <= maxRangeCount) {
                        // We can safely look at all the newly selected stuff.
                        // TODO: replace this with some sort of multi-level set,
                        // or just find a way not to use the old retraction ever
                        // again, so we don't have to do a copy here.
                        toReturn.selected = selected;
                        
                        // Go find what is newly selected and save it
                        toReturn.newlySelected = view.getNewTextPositions(
                            position, retracted);
                            
                        Log::debug() << toReturn.newlySelected.size() <<
                            " new positions found" << std::endl;
                            
                        for(const auto& item : toReturn.newlySelected) {
                            // Copy it all to the selected set too.
                            toReturn.selected.insert(item);
                        }
                        
                    } else {
                        // We can't fill in the sets on the new retraction.
                        toReturn.setsValid = false;
                    }
                }
                
                // We've populated the retraction. Return it.
                return toReturn;
            }
                
        };
        
        /**
         * Represents a retraction that needs to be checked for results.
         */
        struct DPTask {
            /**
             * What number retraction on the left is being considered?
             */
            size_t leftIndex = 0;
            /**
             * What number retraction on the right is being considered?
             */
            size_t rightIndex = 0;
            /**
             * Has this DPTask only ever been retracted right since the root?
             */
            bool isRightEdge = true;
            /**
             * Are we a retraction on the right? If not we might have retracted
             * on the left, or we might be a root.
             */
            bool retractedRight = false;
            /**
             * Are we a retraction on the left? If not we might have retracted
             * on the right, or we might be a root.
             */
            bool retractedLeft = false;
            
            /**
             * Return a copy retracted on the left (if false) or the right (if
             * true).
             */
            inline DPTask retract(bool isRight) const {
                // Make the copy
                DPTask retracted = *this;
                
                if(isRight) {
                    // Retract to the next available spot on the right
                    retracted.rightIndex++;
                    retracted.retractedRight = true;
                    retracted.retractedLeft = false;
                } else {
                    // Retract to the next available spot on the left
                    retracted.leftIndex++;
                    // Which means we have left the right edge
                    retracted.isRightEdge = false;
                    retracted.retractedRight = false;
                    retracted.retractedLeft = true;
                }
                
                // Give back the modified copy
                return retracted;
            }
        }; 
        
        /**
         * Holds the retraction data for the left side.
         */
        std::vector<SideRetractionEntry> leftRetractions;
        /**
         * Holds the retraction data for the right side.
         */
        std::vector<SideRetractionEntry> rightRetractions;
        /**
         * Holds the DPTasks that need to be done.
         */
        std::queue<DPTask> taskQueue;
        
        /**
         * Initialize the DP table with the first left and right retractions.
         */
        inline DPTable(const FMDPosition& left, size_t patternLengthLeft,
            const FMDPosition& right, size_t patternLengthRight,
            const FMDIndexView& view, size_t maxRangeCount) {
            
            leftRetractions.push_back(SideRetractionEntry(left,
                patternLengthLeft, view, maxRangeCount));
                
            rightRetractions.push_back(SideRetractionEntry(right,
                patternLengthRight, view, maxRangeCount));
                
            // Make sure to start at the root
            taskQueue.push(DPTask());
        }
        
        /**
         * Get the left or right retraction entry (as appropriate) for the given
         * retraction number on that side.
         *
         * TODO: find a better way to pass the parameters you actually need to
         * do a retraction.
         */
        inline const SideRetractionEntry& getRetraction(bool isRight,
            size_t retractionNumber, const FMDIndexView& view,
            size_t maxRangeCount) {
            
            // Figure out what side to look at
            auto& retractions = isRight ? rightRetractions : leftRetractions;
            
            while(retractions.size() <= retractionNumber) {
                // Until we are out to that retraction, do all the required
                // retractions.
                retractions.push_back(
                    retractions[retractions.size() - 1].retract(view, 
                    maxRangeCount));
            }
            
            // Now we know it's computed, so return it.
            return retractions[retractionNumber];
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
     * Do the "Activity Selection Problem". Given a vector of [start, end]
     * ranges, find the maximum number of them that you can select without any
     * of the selected ranges overlapping.
     *
     * The input ranges may or may not be sorted. The minimum start value is 0.
     */
    size_t selectActivities(
        std::vector<std::pair<size_t, size_t>> ranges) const;
    
    /**
     * Evaluate the retraction represented by the given DPTask (which may be
     * either a root DPTask or a retraction). If it is possible to ascertain
     * whether the DPTask contains 0, 1, or many shared results, it will return
     * true and a set of 0, 1, or many results. If it is not possible to make
     * that determination (for example, because the number of newly selected
     * positions on the retracted side is too large) it will return false and an
     * unspecified set.
     *
     * Is not responsible for deciding whether a finding of results by one
     * DPTask means that other DPTasks do not need to be executed.
     */
    std::pair<bool, std::set<TextPosition>> exploreRetraction(
        const DPTable::DPTask& task, DPTable& table, const std::string& query,
        size_t queryBase) const;
    
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
