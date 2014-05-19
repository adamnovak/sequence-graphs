#ifndef FMDINDEXITERATOR_HPP
#define FMDINDEXITERATOR_HPP

#include <utility>
#include <deque>

#include "FMDPosition.hpp"

// Forward declaration for circular dependency
class FMDIndex;

/**
 * Represents an iterator over the suffix tree defined by an FMDIndex, yielding
 * the suffix tree suffixes and associated FMDPositions for all suffixes of a
 * given length.
 *
 * The FMDPositions obtained from this iterator are in SA space, but internally
 * the iterator works in BWT space.
 */
class FMDIndexIterator
{
  // TODO: Dramaticallt re-organize this class to not pass data up and down
  // the call stack by putting it in member variables.

public:
    /**
     * Make a new FMDIndexIterator that iterates over the suffix tree of the given
     * FMDIndex at the given depth. If beEnd is set, it skips right to the end
     * to be a 1-past-the-end sentinel. If reportDeadEnds is set, will also
     * include suffixes at depths shorter than the given depth that end in an
     * end of text. For those suffixes, the reverse ranges of the bi-intervals
     * iterated over will not be valid!
     *
     * Depth may not be 0.
     */
    FMDIndexIterator(const FMDIndex& parent, size_t depth, bool beEnd=false,
        bool reportDeadEnds=false);
    
    /**
     * Copy the given FMDIterator.
     */
    FMDIndexIterator(const FMDIndexIterator& toCopy);
    
    /**
     * Pre-increment. Advance and return a reference to ourselves.
     */
    FMDIndexIterator& operator++();
    
    /**
     * Post-increment (which is distinguished by the imaginary int argument).
     * Increment ourselves, and return the previous version by value.
     */
    FMDIndexIterator operator++(int);
    
    /**
     * Dereference operator: get the string suffix and the SA-space FMDPosition
     * corresponding to it.
     */
    std::pair<std::string, FMDPosition> operator*() const;
    
    /**
     * Equality check operator.
     */
    bool operator==(const FMDIndexIterator& other) const;
    
    /**
     * Inquality check operator.
     */
    bool operator!=(const FMDIndexIterator& other) const;
    
    // Default assignment operator is fine.
  
private:
    /**
     * Keep a reference to the parent.
     */
    const FMDIndex& parent;
    
    /**
     * How deep should this iterator go?
     */
    size_t depth;
    
    /**
     * Should this iterator only yield things of the appropriate depth? Or
     * should it also yield dead ends of a shallower depth caused by contexts
     * that run into the end of the text?
     */
    bool reportDeadEnds;
    
    /**
     * This holds a stack for the depth-first search, containing the FMDPosition
     * reached at each level, and the number of the base in BASES that was
     * recursed on to produce that level (so that when we go back to that stack
     * frame we can recurse on the next character). It isn't a real stack
     * because we do need to traverse it to do things like iterator equality.
     */
    std::deque<std::pair<FMDPosition, size_t> > stack;
    
    /**
     * Holds the string corresponding to the current FMDPosition on top of the
     * stack.
     */
    std::string pattern;
    
    /**
     * Keep track of the next value the iterator should return.
     */
    std::pair<std::string, FMDPosition> toYield;
    
    /**
     * Set the next value to be returned when the iterator is dereferenced.
     */
    void yield(std::pair<std::string, FMDPosition> value);
    
    /**
     * Run the depth-first search until we either find a non-empty FMDPosition
     * at the correct depth, or run out of tree.
     */
    void search();
    
    /**
     * Recurse on the base with the given number. Returns true if that was
     * actually done, or false if that would have resulted in an empty range.
     */
    bool recurse(size_t baseNumber);
    
    /**
     * Try recursing with all base numbers, starting from the given one. Take
     * the first one that succeeds and return true, or if none succeed return
     * false.
     */
    bool tryRecurse(size_t baseNumber);
    
    /**
     * Try recursing to a non-empty interval at the iterator's specified depth,
     * starting by exploring the given base number and continuing, if necessary,
     * with base numbers increasing from there. Returns true if a max-depth non-
     * empty interval was found, and false otherwise.
     */
    bool tryRecurseToDepth(size_t baseNumber);
    
    
    /**
     * Pop the top stack frame.
     */
    std::pair<FMDPosition, size_t> pop();
    
    
    
};

#endif
