#ifndef FMDPOSITIONGROUP_HPP
#define FMDPOSITIONGROUP_HPP

#include "FMDPosition.hpp"
#include "FMDIndex.hpp"
#include "FMDIndexView.hpp"
#include <vector>

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
     */
    void extendGreedy(char correctCharacter, size_t maxMismatches);
    
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
     * Holds all the FMDPositions and their mismatch counts.
     */
    std::vector<std::pair<FMDPosition, size_t>> positions;
    
    /**
     * We hold a reference to the FMDIndexView (and thus the FMDIndex) that we
     * work with, so that we can encapsulate the search logic ourselves.
     */
    const FMDIndexView& view;
    
    

};

#endif
