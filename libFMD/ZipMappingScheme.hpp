#ifndef ZIPMAPPINGSCHEME_HPP
#define ZIPMAPPINGSCHEME_HPP

#include "MappingScheme.hpp"
#include "Mapping.hpp"
#include "IntervalIndex.hpp"
#include "Log.hpp"
#include "Matching.hpp"
#include "CreditStrategy.hpp"

#include <iomanip>
#include <queue>
#include <unordered_map>

/**
 * Mapping scheme supporting mapping to graphs, where you zip together maximal
 * matches looking out in both directions from each position. Finds a subset of
 * the unique substrings for each base that the natural mapping scheme finds,
 * but will never map a base that the actual natural to-graph mapping scheme
 * would call as conflicted.
 *
 * Works with either an FMDPosition or an FMDPositionGroup as a SearchType.
 */
template<typename SearchType>
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
     * sided MUMs of a base for that base to map?
     */
    size_t minUniqueStrings = 0;
    
    /**
     * How many total mismatches may be in both sides together of a base's
     * unique context? Only used when templated on FMDPositionGroup.
     */
    size_t mismatchTolerance = 0;
    
    // We need to define this with a new; otherwise our inhereted constructor
    // gets deleted, since it can't default-construct the CreditStrategy without
    // an FMDIndexView. This is to work around a compiler bug
    // <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=62310> which as of 5/14/15
    // is still open.
    // TODO: when that is closed, change this all to:
    // CreditStrategy* credit{view, false};
private:
    CreditStrategy* creditPointer = new CreditStrategy(view, false);
    
public:
    /**
     * Controls credit application. One might want to set its maxMismatches and
     * enabled parameters. Starts disabled by default.
     */
    CreditStrategy& credit = *creditPointer;
    
    /**
     * Destructor for working around bug in GCC.
     */
    inline virtual ~ZipMappingScheme() {
        delete creditPointer;
    }
    
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
            SearchType selection;
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
             * What positions are selected in total by the current search?
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
            inline SideRetractionEntry(const SearchType& unretracted,
                size_t contextLength, const FMDIndexView& view,
                size_t maxRangeCount): selection(unretracted), 
                contextLength(contextLength) {
                
                if(selection.getApproximateNumberOfRanges(view) <=
                    maxRangeCount) {
                
                    // We can visit everything we have selected
                    selected = selection.getTextPositions(view);
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
                
                // Copy and retract our current search
                SearchType retracted = selection;
                size_t retractedContext = retracted.retractRightOnly(view);
                
                // Make a new entry
                SideRetractionEntry toReturn;
                toReturn.selection = retracted;
                toReturn.contextLength = retractedContext;
                
                toReturn.setsValid = setsValid;
                if(setsValid) {
                    // We have sets we can build on, but we need to see if we've
                    // selected too much more stuff.
                    
                    size_t newRanges =
                        retracted.getApproximateNumberOfNewRanges(view,
                        selection);
                        
                    Log::debug() << "Will have " << newRanges <<
                        " new ranges" << std::endl;
                        
                    if(newRanges <= maxRangeCount) {
                        // We can safely look at all the newly selected stuff.
                        // TODO: replace this with some sort of multi-level set,
                        // or just find a way not to use the old retraction ever
                        // again, so we don't have to do a copy here.
                        toReturn.selected = selected;
                        
                        // Go find what is newly selected and save it
                        toReturn.newlySelected = retracted.getNewTextPositions(
                            view, selection);
                            
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
        inline DPTable(const SearchType& left, size_t patternLengthLeft,
            const SearchType& right, size_t patternLengthRight,
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
     * as the characters in the string, and consist of a search and a context
     * length.
     *
     * If reverse is true, outputs results in its vector in reverse order.
     */
    std::vector<std::pair<SearchType, size_t>> findRightContexts(
        const std::string& query, bool reverse = false) const;
        
    /**
     * Returns true if the given unique search result can be extended through
     * the given context string on the other side. That string must end with
     * the base shared with the search result set (which is not used to extend
     * again), and is extended through right to left.
     */
    bool canExtendThrough(SearchType context,
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
     * Given the left and right contexts, finds the shortest unique left and
     * right contexts, and produces an index that can be used to easily find the
     * number of non-overlapping such intervals contained in any range. The
     * index stores, for every query base, the end position of the soonest-
     * ending unique context that starts at or after that base. By looking up
     * your range start, going to 1 later than that position, doing the lookup
     * again, and so on, you can do activity selection fairly easily.
     */
    std::vector<size_t> createUniqueContextIndex(
        const std::vector<std::pair<SearchType, size_t>>& leftContexts,
        const std::vector<std::pair<SearchType, size_t>>& rightContexts) const;
    
    /**
     * Do activity selection using an index from createUniqueContextIndex. Given
     * the start and end positions of a range (inclusive), and an index, finds
     * the maximum number of non-overlapping activities (i.e. unique contexts)
     * that can fit in the range. If threshold is specified, finds that number
     * or the threshold, whichever is smaller.
     */
    size_t selectActivitiesIndexed(size_t start, size_t end,
        const std::vector<size_t>& index, size_t threshold = (size_t) -1) const;
    
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
        const typename DPTable::DPTask& task, DPTable& table,
        const std::string& query, size_t queryBase) const;
    
    /**
     * Explore all retractions of the two searches. Return either an empty
     * Mapping (if no explored retraction finds overlapping positions between
     * the two sides), a Mapping with its TextPosition and context lengths set
     * (if we find exactly one such overlap), or an empty Mapping (if we find
     * multiple overlaps).
     */
    Mapping exploreRetractions(const SearchType& left,
        size_t patternLengthLeft, const SearchType& right,
        size_t patternLengthRight, const std::string& query,
        size_t queryBase) const;    
    
};

// Now we have to do all the template definitions, since they need to be in the
// headers unless they're specializations.

template<typename SearchType>
size_t ZipMappingScheme<SearchType>::selectActivities(
    std::vector<std::pair<size_t, size_t>> ranges) const {
    
    stats.add("activitySelectionRuns", 1);
    
    // First we have to sort the ranges by end, ascending.
    std::sort(ranges.begin(), ranges.end(), [](std::pair<size_t, size_t> a,
        std::pair<size_t, size_t> b) -> bool {
        
        // We'll give true if the first range ends before the second.
        // We also return true if there's a tie but a starts first.
        return (a.second < b.second) || (a.second == b.second &&
            a.first < b.first);
        
    });
    
    // How many non-overlapping things have we found?
    size_t found = 0;
    // What's the minimum start time we can accept?
    size_t nextFree = 0;
    
    for(const auto& range : ranges) {
        Log::trace() << "Have " << range.first << " - " << range.second <<
            std::endl;
        // For each range (starting with those that end soonest)...
        if(range.first >= nextFree) {
            // Take it if it starts late enough
            
            Log::trace() << "Taking " << range.first << " - " << range.second <<
                std::endl;
            
            nextFree = range.second + 1;
            found++;
            
            stats.add("activitiesSelected", 1);
        }
    }

    return found;
}

template<typename SearchType>
std::vector<size_t> ZipMappingScheme<SearchType>::createUniqueContextIndex(
    const std::vector<std::pair<SearchType, size_t>>& leftContexts,
    const std::vector<std::pair<SearchType, size_t>>& rightContexts) const {
    
    // Create vectors of the minimum unique left and right context lengths, or
    // (size_t) -1 if no such lenght exists.
    std::vector<size_t> minUniqueLeftContexts;
    std::vector<size_t> minUniqueRightContexts;
    
    for(const auto& rangeAndLength : leftContexts) {
        // Find the min unique lengths for the left contexts.
        
        SearchType position = rangeAndLength.first;
        size_t length = rangeAndLength.second;
        
        // What's the min unique length for this base? We haven't found any yet.
        size_t minUniqueLength = (size_t) -1;
        
        while(position.isUnique(view)) {
            // Retract to a point where we may not be unique
            length = position.retractRightOnly(view);
            // We know we are unique at 1 more character than that
            minUniqueLength = length + 1;
        }
        
        if(minUniqueLength == 0) {
            throw std::runtime_error("Too much retracting!");
        }
        
        minUniqueLeftContexts.push_back(minUniqueLength);
        
        Log::trace() << "Min left context " <<
            minUniqueLeftContexts.size() - 1 << " is " << minUniqueLength <<
            " vs " << length << " selecting " <<
            position.getTextPositions(view).size() << std::endl;
        
    }
    
    for(const auto& rangeAndLength : rightContexts) {
        // And again for the right contexts. TODO: make a function to not repeat
        // all this code.
        
        SearchType position = rangeAndLength.first;
        size_t length = rangeAndLength.second;
        
        // What's the min unique length for this base? We haven't found any yet.
        size_t minUniqueLength = (size_t) -1;
        
        while(position.isUnique(view)) {
            // Retract to a point where we may not be unique
            length = position.retractRightOnly(view);
            // We know we are unique at 1 more character than that
            minUniqueLength = length + 1;
        }
        
        minUniqueRightContexts.push_back(minUniqueLength);
    }
    
    // Now we're going to make a sorted list of (start, end) pairs for these
    // one-sided MUSes
    std::vector<std::pair<size_t, size_t>> uniqueContexts;
    
    for(size_t i = 0; i < minUniqueLeftContexts.size(); i++) {
        // For every base
        
        if(minUniqueLeftContexts[i] != (size_t) -1) {
            // If it has a left context,  Save this one-sided MUS
            uniqueContexts.emplace_back(
                i - minUniqueLeftContexts[i] + 1, i);
      
            if(minUniqueLeftContexts[i] == 0) {
                throw std::runtime_error("Left context runs backwards!");
            }
        }
        
        if(minUniqueRightContexts[i] != (size_t) -1) {
            // If it has a right context, save this one-sided MUS
            uniqueContexts.emplace_back(i,
                i + minUniqueRightContexts[i] - 1);
                
            if(minUniqueRightContexts[i] == 0) {
                throw std::runtime_error("Right context runs backwards!");
            }
        }
    }
    
    // The "first" range is the one that ends earliest.
    auto order = ([](std::pair<size_t, size_t> a, 
        std::pair<size_t, size_t> b) -> bool {
        
        // We'll give true if the first range ends before the second.
        // We also return true if there's a tie but a starts first.
        return (a.second < b.second) || (a.second == b.second &&
            a.first < b.first);
        
    });
    
    // Sort the unique contexts by start, descending, so we can pop stuff off
    // the front when we reach its start position.
    std::sort(uniqueContexts.begin(), uniqueContexts.end(),
        [](std::pair<size_t, size_t> a, std::pair<size_t, size_t> b) -> bool {
            // Define a total ordering
            return (a.first > b.first) || (a.first == b.first &&
                a.second > b.second);
    });
    
    // Go through the unique contexts from right to left, and save for each base
    // the end position of the earliest-ending unique context that starts at or
    // after it.
    std::vector<size_t> endOfBestActivity(leftContexts.size());
    
    // We're going to iterate through the unique contexts in order of decreasing
    // start index.
    auto nextActivity = uniqueContexts.begin();
    
    // We'll keep an iterator to the one that starts here or to the right and
    // ends soonest, but for now that's nothing.
    auto bestActivity = uniqueContexts.end();
    
    // We go from right to left. We have a current best (i.e. earliest-ending)
    // interval. When we reach the start point of an interval, we see if that
    // interval ends before our current best one. If so we replace the best
    // interval. If not we keep it and ignore the new interval.
    
    for(size_t i = leftContexts.size() - 1; i != (size_t) -1; i--) {
        // For each base from right to left
        
        while(nextActivity != uniqueContexts.end() && 
            (*nextActivity).first >= i) {
            // For every activity starting at or after here
            
            Log::trace() << "Considering activity " << (*nextActivity).first <<
                "-" << (*nextActivity).second << " which starts at or after " <<
                i << std::endl;
            
            // This activity starts here or later
            if(bestActivity == uniqueContexts.end() ||
                (*nextActivity).second < (*bestActivity).second) {
                
                // This is a sooner-ending activity than any we have so far.
                // Take it.
                bestActivity = nextActivity;
                
                Log::trace() << "It is the soonest ending" << std::endl;
                
            }
            // We accepted or rejected this activity, so go on to the next
            // lefter one to consider.
            ++nextActivity;
        }
        
        // Now bestActivity points to the soonest-ending activity that starts
        // here or later.
        if(bestActivity == uniqueContexts.end()) {
            // There is no such activity. Use (size_t) -1 as a no-such-value
            // value.
            endOfBestActivity[i] = (size_t) -1;
            
            Log::trace() << "No soonest ending activity starting at or after " <<
                i << std::endl;
        } else {
            // We found such an activity, so record its endpoint. We can go one
            // right from there and look up the end of the non-overlapping
            // activity that we would chain with.
            endOfBestActivity[i] = (*bestActivity).second;
            
            Log::trace() << "Soonest ending activity starting at or after " <<
                i << " ends at " << endOfBestActivity[i] << std::endl;
        }
    }
    
    // Giev back the index vector we made, which simplifies activity selection
    // problems immensely.
    return endOfBestActivity;
}

template<typename SearchType>
size_t ZipMappingScheme<SearchType>::selectActivitiesIndexed(size_t start,
    size_t end, const std::vector<size_t>& index, size_t threshold) const {
    
    // What base are we at?
    size_t position = start;
    // How many non-overlapping things have we found?
    size_t found = 0;

    while(found < threshold && position < index.size()) {
        // Look for activities until we find enough. A threshold of (size_t) -1
        // wil have us look forever. Also make sure we don't sneak out of the
        // index with that +1 later.
    
        // Jump to the end of the next activity we want to do.
        position = index[position];
        
        if(position > end) {
            // This activity runs too long and we can't do it. This also handles
            // (size_t) -1, the no-more-activities value, which is larger than
            // any other size_t.
            break;
        }
        
        // Otherwise we found an activity that we can do.
        found++;
        
        // Look for activities that start after this one ends.
        position++;
    }
    
    // Say how many activities were found.
    return found;
}

template<typename SearchType>
Mapping ZipMappingScheme<SearchType>::exploreRetractions(
    const SearchType& left, size_t patternLengthLeft, const SearchType& right,
    size_t patternLengthRight, const std::string& query,
    size_t queryBase) const {

    stats.add("basesAttempted", 1);

    // List the unique TextPosition we found, if we found one. Holds more than
    // one TextPosition (though not necessarily all of them) if we're ambiguous,
    // and none if we have no results. Holds TextPositions for right contexts.
    std::set<TextPosition> found;
    
    // What are the max left and right contexts that were used to map the base
    // anywhere? Note that they may not be part of the same string that maps the
    // base there, but they would have to be if we accounted for crossover.
    size_t maxLeftContext = 0;
    size_t maxRightContext = 0;

    // Make a DP table for this base. TODO: make DPTable remember the extra
    // parameters.
    DPTable table(left, patternLengthLeft, right, patternLengthRight,
        view, maxRangeCount);
        
    // I can do my DP by always considering retracting on the left from a state,
    // and only considering retracting on the right on the very edge of the
    // space (i.e. if the left is still full-length).
        
    // Make a map of the minimum right context (r) we will ever have to explore
    // for any left context (l) this size or smaller. We can use lower_bound for
    // the lookups so we only have to insert l,r pairs where we have overlap.
    std::map<size_t, size_t> minRightContext;
        
    // A function to see if we need to text a certain combination of left and
    // right context lengths, or if it's a more general context of one we
    // already had overlapping results for
    auto needToTest = [&](const typename DPTable::DPTask& task) {
        // Work out the context lengths you would have on the left and right for
        // this task (and probably actually compute the retracted sets too).
        const typename DPTable::SideRetractionEntry& left = table.getRetraction(
            false, task.leftIndex, view, maxRangeCount);
        const typename DPTable::SideRetractionEntry& right =
            table.getRetraction(true, task.rightIndex, view, maxRangeCount);
        
        if(left.contextLength == 0 || right.contextLength == 0) {
            // Way too short
            return false;
        }
        
        // Find the min right context for things with left contexts greater than
        // or equal to this one.
        auto minRightIterator = minRightContext.lower_bound(left.contextLength);
        
        if(minRightIterator != minRightContext.end() &&
            right.contextLength <= (*minRightIterator).second) {
            // We've retracted far enough that we don't need to process this
            // retraction, since it's a less general context than one we already
            // found results for (on the left it's a lower bound, and we just
            // saw it's a lower bound on the right).
            
            Log::debug() << "Skipping " << left.contextLength << ", " <<
                right.contextLength << " as it is covered by " <<
                (*minRightIterator).first << ", " <<
                (*minRightIterator).second << std::endl;
            
            // Note that there was a covered retraction.
            stats.add("retractionCovered", 1);
            
            return false;
        }
        
        // We do need to process this retraction.
        return true;
    };
        
    while(table.taskQueue.size() > 0) {
    
        // Grab the task in the table
        const typename DPTable::DPTask& task = table.taskQueue.front();
    
        if(!needToTest(task)) {
            // Skip queued tasks that have becomne redundant.
            table.taskQueue.pop();
            continue;
        }
    
        // Grab the left and right retraction entries.
        // TODO: Stop repeating this code somehow.
        const typename DPTable::SideRetractionEntry& left = table.getRetraction(
            false, task.leftIndex, view, maxRangeCount);
        const typename DPTable::SideRetractionEntry& right =
            table.getRetraction(true, task.rightIndex, view, maxRangeCount);
            
        // Which we needed to pull out the context lengths on which we may be
        // mapping.
        size_t leftContext = left.contextLength;
        size_t rightContext = right.contextLength;
        
        // Which we in turn need to calculate the max left and right context
        // lengths observed for the base.
        maxLeftContext = std::max(maxLeftContext, leftContext);
        maxRightContext = std::max(maxRightContext, rightContext);
    
        // Run the first task in the queue, possibly adding more, and getting
        // some results.
        auto flagAndSet = exploreRetraction(task, table, 
            query, queryBase);
            
        // Drop that task we just did. We can't use task anymore now, or left or
        // right.
        table.taskQueue.pop();
        
        if(!flagAndSet.first) {
            // We encountered something too hard to do.
            stats.add("tooHardRetraction", 1);
            
            if(giveUpIfHard) {
                // We have to abort mapping.
                Log::debug() << "Aborting mapping base " << queryBase <<
                    " because it is too hard." << std::endl;
                // Return an empty mapping.
                return Mapping();
            } else {
                Log::debug() << "Ignoring hard task" << std::endl;
            }
        }
        
        for(auto position : flagAndSet.second) {
            // Put in all the places we found matchings to to.
            found.insert(position);
        }
        
        if(found.size() > 1) {
            // We're already ambiguous. Short circuit.
            Log::debug() << "Already ambiguous, not retracting any more" <<
                std::endl;
            stats.add("ambiguous", 1);
            return Mapping();
        }
        
        if(!useRetraction) {
            // Only explore the very first task, which is the root and not
            // actually retracting on either side.
            break;
        }
        
    }
    
    Log::debug() << "Found " << found.size() << " locations" << std::endl;
    Log::debug() << "Used " << maxLeftContext << ", " << maxRightContext <<
        " context" << std::endl;
        
    if(found.size() == 1) {
        // We mapped to one place, on these contexts.
        stats.add("unambiguous", 1);
        return Mapping(*(found.begin()), maxLeftContext, maxRightContext);
    } else {
        // We mapped to nowhere, because we're ambiguous. TODO: log
        // ambiguousness.
        stats.add("ambiguous", 1);
        return Mapping();
    }
} 

template<typename SearchType>
void ZipMappingScheme<SearchType>::map(const std::string& query,
    std::function<void(size_t, TextPosition)> callback) const {
    
    // Get the right contexts
    Log::info() << "Looking for right contexts..." << std::endl << std::flush;
    auto rightContexts = findRightContexts(query, false);
    
    // And the left contexts (which is the reverse of the contexts for the
    // reverse complement, which we can produce in backwards order already).
    Log::info() << "Flipping query..." << std::endl << std::flush;
    std::string complement = reverseComplement(query);
    
    Log::info() << "Looking for left contexts..." << std::endl << std::flush;
    auto leftContexts = findRightContexts(complement, true);
    
    // This is going to hold, for each query base, the endpoint of the soonest-
    // ending minimally unique context that starts at or after that base.
    Log::debug() << "Creating unique context index" << std::endl << std::flush;
    auto uniqueContextIndex = createUniqueContextIndex(leftContexts,
        rightContexts);
    
    
    Log::info() << "Exploring retractions..." << std::endl;
    
    // We're going to store up all the potential mappings, and then filter them
    // down. This stores true and the mapped-to TextPosition if a base would map
    // before filtering, and false and an undefined TextPosition otherwise.
    std::vector<Mapping> mappings;
    
    // For each pair, figure out if the forward and reverse searches select
    // one single consistent TextPosition.
    for(size_t i = 0; i < query.size(); i++) {
    
        Log::debug() << "Base " << i << " = " << query[i] << " (+" << 
            leftContexts[i].second << "|+" << rightContexts[i].second << 
            ") selects " << leftContexts[i].first << " and " << 
            rightContexts[i].first << std::endl;
            
        // Go look at these two contexts in opposite directions, and all their
        // (reasonable to think about) retractions, and see whether this base
        // belongs to 0, 1, or multiple TextPositions.
        Mapping mapping = exploreRetractions(
            leftContexts[i].first, leftContexts[i].second,
            rightContexts[i].first, rightContexts[i].second, query, i);
        
        // TODO: If we can't find anything, try retracting a few bases on one or
        // both sides until we get a shared result.
        
        if(mapping.isMapped()) {
            // We map!
            Log::debug() << "Index " << i << " maps to " << mapping << std::endl;
        } else {
            // Too few results until we retracted back to too many
            Log::debug() << "Index " << i << " is not mapped." << std::endl;
       }
       // Save the mapping
       mappings.push_back(mapping);
        
    }
    
    Log::info() << "Applying filter..." << std::endl << std::flush;
    
    // We're going to filter everything first and then do the callbacks.
    std::vector<Mapping> filtered;
    
    for(size_t i = 0; i < mappings.size(); i++) {
        // Now we have to filter the mappings
        
        // Grab the mapping for this query index.
        const Mapping& mapping = mappings[i];
        
        if(!mapping.isMapped()) {
            // Skip unmapped query bases
            filtered.push_back(mapping);
            continue;
        }
        
        // What range of bases are within our left and right MUMs?
        size_t rangeStart = i - mapping.getLeftMaxContext() + 1;
        size_t rangeEnd = i + mapping.getRightMaxContext() - 1; // Inclusive
        
        // Do activity selection
        size_t nonOverlapping = selectActivitiesIndexed(rangeStart, rangeEnd,
            uniqueContextIndex, minUniqueStrings);
            
        if(nonOverlapping < minUniqueStrings) {
            Log::debug() <<
                "Dropping mapping " << i <<
                " due to having too few unique strings." << std::endl;
            // Report a non-mapping Mapping
            stats.add("filterFail", 1);
            filtered.push_back(Mapping());
        } else {
            // Report all the mappings that pass.
            stats.add("filterPass", 1);
            filtered.push_back(mapping);
        }
        
    }
    
    if(credit.enabled) {
        // If credit isn't disabled, run it. The cannonical check for enabled-
        // ness is in the CreditStrategy itself, but no point running it if it's
        // going to do nothing.
        Log::info() << "Applying credit..." << std::endl;
        credit.applyCredit(query, filtered);
    }
    
    Log::info() << "Sending callbacks..." << std::endl << std::flush;
    for(size_t i = 0; i < filtered.size(); i++) {
        if(filtered[i].isMapped()) {
            // Send a callback for everything that passed the filter.
            callback(i, filtered[i].getLocation());
        }
    }
    
}

#endif
