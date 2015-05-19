#ifndef FMDPOSITIONGROUP_HPP
#define FMDPOSITIONGROUP_HPP

#include "FMDPosition.hpp"
#include "FMDIndex.hpp"
#include "FMDIndexView.hpp"
#include <vector>
#include <set>
#include <deque>

/**
 * Represents a collection of FMDPositions participating in some kind of
 * mismatch-tolerant breadth-first search. Developed for graph credit (which can
 * be greedy and not always have to consider mismatches at every point they are
 * possible) but could also work for an all-mismatches search.
 *
 * Has the ability to be initialized with a bunch of FMDPositions (which could
 * correspond to certain un-located graph nodes).
 */
class FMDPositionGroup {

public:

    /**
     * Start out an empty FMDPositionGroup.
     */
    FMDPositionGroup();

    /**
     * Start out a new FMDPositionGroup containing the given FMDPositions.
     * Retract must not be used, because since we started with unknown numbers
     * of characters searched, we won't know how to pick up mismatches
     * correctly.
     *
     * The intervals must not overlap, and must all correspond to search strings
     * of the same length if retract is ever going to be used.
     *
     * TODO: Look into checking LCP.
     */
    FMDPositionGroup(const std::vector<FMDPosition>& positions);
    
    /**
     * Create an FMDPositionGroup selecting an entire view's index.
     */
    FMDPositionGroup(const FMDIndexView& view);
    
    // Default copy constructor/move constructor/destructor/assignment operator
    // are all OK.
        
    /**
     * Extend all the FMDPositions left, using the given view. If there are any
     * results extending with the correct character, take those. Otherwise,
     * extend with all of the mismatching characters. If any FMDPosition
     * accumulates too many mismatches, drop it.
     *
     * The extension is "greedy" because if it's possible to have a match, it
     * takes the match, even if that would lead it down a path with, overall,
     * more mismatches.
     *
     * Returns true if an exact match exists anywhere, and false if a mismatch
     * had to be used.
     */
    bool extendGreedy(const FMDIndexView& view, char correctCharacter,
        size_t maxMismatches);
        
    /**
     * Extend all FMDPositions left, using the given view. Explores all possible
     * mismatches, if they yield results and do not end up costing too much
     * total.
     */
    bool extendFull(const FMDIndexView& view, char correctCharacter,
        size_t maxMismatches);
        
    /**
     * Extend left with no mismatches.
     */
    void extendLeftOnly(const FMDIndexView& view, char character);
    
    /**
     * Retract to the given search string length, on the right.
     */
    void retractRightOnly(const FMDIndexView& view, size_t newLength);
    
    /**
     * Retract until more BWT positions are selected. They may not actually be
     * different merged ranges or masked in.
     */
    size_t retractRightOnly(const FMDIndexView& view);
    
    /**
     * How many mismatches are used max by any FMDPosition in the group?
     */
    size_t mismatchesUsed() const;
    
    /**
     * Is nothing selected under the given view?
     */
    bool isEmpty(const FMDIndexView& view) const;
    
    /**
     * Is exactly one merged position selected under the given view?
     */
    bool isUnique(const FMDIndexView& view) const;
    
    /**
     * Is more than one merged position selected under the given view?
     */
    bool isAmbiguous(const FMDIndexView& view) const;
    
    /**
     * Get the unique TextPosition selected under the viven view. isUnique()
     * must be true.
     */
    TextPosition getTextPosition(const FMDIndexView& view) const;
    
    /**
     * Provides an overestimate of the number of ranges selected under the given
     * view.
     */
    size_t getApproximateNumberOfRanges(const FMDIndexView& view) const;
    
    /**
     * Provides an overestimate of the number of new ranges selected under the
     * given view, relative to the given old FMDPosition.
     *
     * Neither FMDPositionGroup may contain overlapping intervals.
     */
    size_t getApproximateNumberOfNewRanges(const FMDIndexView& view,
        const FMDPositionGroup& old);
        
    /**
     * Get the TextPositions selected under the given view.
     */
    std::set<TextPosition> getTextPositions(const FMDIndexView& view) const;
    
    /**
     * Get the new TextPositions selected under the given view, relative to the
     * given old FMDPositionGroup. Some of these may already have been selected
     * in the old group.
     *
     * Neither FMDPositionGroup may contain overlapping intervals.
     */
    std::set<TextPosition> getNewTextPositions(const FMDIndexView& view,
        const FMDPositionGroup& old) const;
    
protected:

    /**
     * Keep track of an FMDPosition with its mismatch information.
     */
    struct AnnotatedFMDPosition {
        /**
         * What is the actual BWT interval?
         */
        FMDPosition position;
        
        /**
         * How many mismatches were used searching it?
         * TODO: Redundant with mismatch offset queue length.
         */
        size_t mismatches;
        
        /**
         * How many characters are searched (extended and not retracted)?
         */
        size_t searchedCharacters;
        
        /**
         * Stores first the distance from the retractable (right) end of the
         * searched string to the first mismatch, followed by the distance from
         * that mismatch to the next, and so on.
         */
        std::deque<size_t> mismatchOffsets;
        
        /**
         * Constructor to wrap up an FMDPosition.
         */
        inline AnnotatedFMDPosition(const FMDPosition& position):
            position(position), mismatches(0), searchedCharacters(0),
            mismatchOffsets() {
            
            // Nothing to do!
        }
        
        /**
         * Constructor for extension of a parent.
         */
        inline AnnotatedFMDPosition(const FMDPosition& position,
            const AnnotatedFMDPosition& parent,
            bool isMismatch): position(position), 
            mismatches(parent.mismatches + isMismatch), 
            searchedCharacters(parent.searchedCharacters + 1), 
            mismatchOffsets(parent.mismatchOffsets) {
            
            if(isMismatch) {
                // We added a mismatch at the left.
                if(mismatchOffsets.size() > 0) { 
                    // We come after a mismatch that's already there. After
                    // removing it, how many more characters would we need to
                    // retract to remove this one?
                    mismatchOffsets.push_back(searchedCharacters -
                        mismatchOffsets.back());
                } else {
                    // There are no mismatches before this new one. How many
                    // characters would we have to retract to remove it?
                    mismatchOffsets.push_back(searchedCharacters);
                }
            }
            
        }
        
        /**
         * Constructor for retraction from a parent.
         */
        inline AnnotatedFMDPosition(const FMDPosition& position,
            const AnnotatedFMDPosition& parent, size_t newLength): 
            position(position), mismatches(parent.mismatches), 
            searchedCharacters(newLength), 
            mismatchOffsets(parent.mismatchOffsets) {
            
            // How many characters were retracted?
            size_t droppedLength = parent.searchedCharacters -
                searchedCharacters;
            
            while(mismatches > 0 && droppedLength > 0) {
            
                if(mismatchOffsets.front() > droppedLength) {
                    // If there were more than this many characters left to the
                    // next mismatch, just charge for them.
                    mismatchOffsets.front() -= droppedLength;
                } else {
                    // If there were fewer, delete that mismatch and then see if
                    // we get the next one too.
                    droppedLength -= mismatchOffsets.front();
                    mismatchOffsets.pop_front();
                    mismatches--;
                }
            }
        }
        
        /**
         * Equality comparison so this can go in a set.
         */
        inline bool operator==(const AnnotatedFMDPosition& other) const {
            return position == other.position &&
                mismatches == other.mismatches &&
                searchedCharacters == other.searchedCharacters;
            // We don't have to check the offset lists because they are
            // determined by our other parameters, assuming we have a consistent
            // history for the group in terms of what was searched.
        }
        
        /**
         * Order comparison so this can go in a set.
         */
        inline bool operator<(const AnnotatedFMDPosition& other) const {
            return position < other.position || 
                (position == other.position && (mismatches < other.mismatches ||
                (mismatches == other.mismatches && searchedCharacters <
                other.searchedCharacters)));
            // We don't have to check the offset lists because they are
            // determined by our other parameters, assuming we have a consistent
            // history for the group in terms of what was searched.
        }
    };

    /**
     * Holds all the FMDPositions and their mismatch counts.
     */
    std::set<AnnotatedFMDPosition> positions;
};

#endif
