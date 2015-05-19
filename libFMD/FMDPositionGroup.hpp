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
     * Start out an empty FMDPositionGroup on the given view of an index.
     */
    FMDPositionGroup(const FMDIndexView& view);

    /**
     * Start out a new FMDPositionGroup containing the given FMDPositions.
     * Retract must not be used, because since we started with unknown numbers
     * of characters searched, we won't know how to pick up mismatches
     * correctly.
     *
     * TODO: Look into checking LCP.
     */
    FMDPositionGroup(const FMDIndexView& view,
        const std::vector<FMDPosition>& positions);
        
    /**
     * Extend all the FMDPositions. If there are any results extending with the
     * correct character, take those. Otherwise, extend with all of the
     * mismatching characters. If any FMDPosition accumulates too many
     * mismatches, drop it.
     *
     * The extension is "greedy" because if it's possible to have a match, it
     * takes the match, even if that would lead it down a path with, overall,
     * more mismatches.
     *
     * Returns true if an exact match exists anywhere, and false if a mismatch
     * had to be used.
     */
    bool extendGreedy(char correctCharacter, size_t maxMismatches);
    
    /**
     * Retract by a single base.
     *
     * TODO: work out how to retract the minimum necessary for one FMDPosition
     * to show more results.
     */
    void retractOne();
    
    /**
     * Are there no non-empty FMDPositions in the group?
     */
    bool isEmpty() const;
    
    /**
     * Are all non-empty FMDPositions selecting the exact same graph position?
     */
    bool isUnique() const;
    
    /**
     * Get the unique TextPosition selected. isUnique() must be true.
     */
    TextPosition getTextPosition() const;
    
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
    
    /**
     * We hold a reference to the FMDIndexView (and thus the FMDIndex) that we
     * work with, so that we can encapsulate the search logic ourselves.
     */
    const FMDIndexView& view;
    
    

};

#endif
