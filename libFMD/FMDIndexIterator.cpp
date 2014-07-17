#include <string>
#include <utility>
#include <stdexcept>
#include "FMDIndexIterator.hpp"
#include "FMDIndex.hpp"
#include "util.hpp"
#include "debug.hpp"

// Stuff for FMDIndexIterators that traverse the suffix tree.

FMDIndexIterator::FMDIndexIterator(const FMDIndex& parent, size_t depth,
    bool beEnd, bool reportDeadEnds): 
    parent(parent), depth(depth), reportDeadEnds(reportDeadEnds), stack(),
    pattern() {
    
    // By default we start out with empty everything, which is what we should
    // have at the end.

    if(!beEnd)  {
        // Start at the beginning.
        search();
    }
}

FMDIndexIterator::FMDIndexIterator(const FMDIndexIterator& toCopy): 
    parent(toCopy.parent), depth(toCopy.depth), 
    reportDeadEnds(toCopy.reportDeadEnds), stack(toCopy.stack),
    pattern(toCopy.pattern) {
    
    // Already made a duplicate stack. Nothing to do.
}

FMDIndexIterator& FMDIndexIterator::operator++() {
    // Advance
    search();

    // Return ourselves.
    return *this;
}

FMDIndexIterator FMDIndexIterator::operator++(int) {
    // Copy ourselves
    FMDIndexIterator copy(*this);

    // Advance
    search();

    // Return the copy
    return copy;
}

std::pair<std::string, FMDPosition> FMDIndexIterator::operator*() const  { 
    // Just grab the value we stored.
    return toYield;
}

void FMDIndexIterator::yield(std::pair<std::string, FMDPosition> value) {
    // Save the value for returning when dereferenced.
    toYield = value;
}

bool FMDIndexIterator::operator==(const FMDIndexIterator& other) const {
    return 
        // We have the same parent addresses
        &parent == &(other.parent) && 
        // And go to the same depth
        depth == other.depth && 
        // And both report dead ends or not
        reportDeadEnds == other.reportDeadEnds &&
        // And are at the same depth
        stack.size() == other.stack.size() && 
        // And followed the same path to get there
        std::equal(stack.begin(), stack.end(), other.stack.begin()) && 
        // And we have the same string pattern (which the above should imply)
        pattern == other.pattern;
}

bool FMDIndexIterator::operator!=(const FMDIndexIterator& other) const {
    // Just use the equality check.
    return !(*this == other);
}

void FMDIndexIterator::search() {
    // TODO: Unify with tryRecurseToDepth by making its topDepth a parameter.

    if(stack.size() == 0)  {
        // Try to recurse down to depth, starting with base 0 at the root
        // branch.
        tryRecurseToDepth(0);

        // Now we are either at the correct depth at the leftmost suffix tree
        // node, or in the state we started in (which is equal to end). So we
        // are done.
    } else if(stack.size() == depth) {
        // We were at the right depth, meaning we don't need to recurse further.

        do {
            // Pop the bottom frame (already explored).
            std::pair<FMDPosition, size_t> lastFrame = pop();

            // Try recursing to the next thing at that same level and then down
            // to the required depth.
            if(tryRecurseToDepth(lastFrame.second + 1)) {
                // We got something at the required depth, or ran into a shorter
                // end of text suffix to report.
                return;
            }

            // Otherwise, we didn't get something to yield in this subtree. That
            // means this subtree has been exhausted and we should pop and try
            // the next one over, which we will do on the next loop iteration.
        } while(stack.size() > 0);

        // If we get here, we have gone all the way back up and popped and tried
        // replacing everything with no more results. That means we have
        // finished the search, and should be equal to end.

    } else if(reportDeadEnds) {
        // We weren't at the right depth; we were at a node that happened to
        // show up followed by a text end somewhere. Continue recursing from
        // here.

        if(tryRecurseToDepth(0)) {
            // We can just go down from here.
            return;
        }

        // Otherwise we need to pop up and continue the main recursion loop from
        // above. TODO: Unify.

        do {
            // Pop the bottom frame (already explored).
            std::pair<FMDPosition, size_t> lastFrame = pop();

            // Try recursing to the next thing at that same level and then down
            // to the required depth.
            if(tryRecurseToDepth(lastFrame.second + 1)) {
                // We got something at the required depth, or another shorter
                // suffix followed by end of text.
                return;
            }

            // Otherwise, we didn't get something at the required depth in this
            // subtree. That means this subtree has been exhausted and we should
            // pop and try the next one over, which we will do on the next loop
            // iteration.
        } while(stack.size() > 0);    

    } else {
        // We broke something.
        throw std::runtime_error("Iterator was at wrong depth");
    }
}

bool FMDIndexIterator::recurse(size_t baseNumber) {
    if(baseNumber >= NUM_BASES) {
        // Don't even try out-of-range bases.
        return false;
    }

    // What will we extend to?
    FMDPosition extension;

    if(stack.size() == 0) {
        // Our "extension" is just starting with this base.
        extension = parent.getCharPosition(ALPHABETICAL_BASES[baseNumber]);
    }
    else {
        // Work out what we would select if we extended forwards with this
        // letter (i.e. appended it to the suffix).
        extension = stack.back().first;
        parent.extendFast(extension, ALPHABETICAL_BASES[baseNumber], false);
    }



    if(extension.isEmpty()) {
        // This would be a suffix that doesn't appear.
        return false;
    }

    // This would not be an empty place. Go there.
    // Add a stack frame
    stack.push_back(std::make_pair(extension, baseNumber));
    // And record the change to the pattern.
    pattern.push_back(ALPHABETICAL_BASES[baseNumber]);

    return true;
}

bool FMDIndexIterator::tryRecurse(size_t baseNumber) {
  // Try every base number after the given starting one until we either recurse
  // or run out.
  for(; baseNumber < NUM_BASES && !recurse(baseNumber); baseNumber++);
  
  // If we didn't run out, we must have successfully recursed.
  return baseNumber < NUM_BASES;
}

bool FMDIndexIterator::tryRecurseToDepth(size_t baseNumber) {
    // What depth should we not go above?
    size_t topDepth = stack.size();

    while(stack.size() < depth) {
        // Until we get to the right depth...

        // We never go into an empty thing.
        // So we must be at a shallower depth. Try to go deeper.

        if(tryRecurse(baseNumber)) {
            // We made it down. Reset baseNumber
            baseNumber = 0;

            if(reportDeadEnds && stack.size() < depth) {
                // If we lose some range here (i.e. some positions are followed
                // by end of text and don't show up in an extension with any
                // base), we have to stop so that we yield them. Those positions
                // will be the first ones in our range if they exist.
                
                // TODO: replace with asking the BWT about the stop character

                // See what we would get if we extended. TODO: Can we show that
                // this will always work? It seems to work when there are none
                // of the first base, and should always work when there is some
                // of the first base, but I'm not sure it won't break.
                FMDPosition extension = stack.back().first;
                parent.extendFast(extension, ALPHABETICAL_BASES[0], false);

                if(extension.getForwardStart() != 
                    stack.back().first.getForwardStart()) {
                
                    // There are some suffixes that come into this node and
                    // leave before the first real base. They must end the text.

                    INFO(
                        std::cout << "End of text: " << pattern << "$" << 
                            std::endl;
                        std::cout << extension << " vs. " << 
                        stack.back().first << std::endl;
                    )

                    // We got to a place we want to yield.

                    // Grab the FMDPosition we want to return
                    FMDPosition patternPosition = stack.back().first; 

                    // Fix it up to indicate only the part we aren't covering.
                    // Forward start is going to be the same, but it is going to
                    // run until the start of extension. Subtract 1 to keep this
                    // as an offset where 0 = a 1-base interval.
                    patternPosition.setEndOffset(extension.getForwardStart() - 
                        patternPosition.getForwardStart() - 1);

                    // TODO: The reverse start can't be moved sanely and is thus
                    // going to be invalid (we can't search anchored to text
                    // start).

                    // Grab the pattern and the FMDPosition from the top of the
                    // stack to go with it. Yield the pair.
                    yield(std::make_pair(pattern, patternPosition));

                    return true;
                }

            }

        } else {
            // We can't go to anything nonempty down. So we should try going up.

            if(stack.size() == topDepth) {
                // We don't want to go any higher than this. Give up on finding
                // a nonempty thing at the right depth in this subtree.
                return false;
            }

            // Otherwise, go up and continue on the next branch.
            baseNumber = pop().second + 1;
        }
    }

    // We got to the right depth, and we never go to empty things.

    // Grab the pattern and the FMDPosition from the top of the stack to go with
    // it. Yield the pair.
    yield(std::make_pair(pattern, stack.back().first));

    return true;
}

std::pair<FMDPosition, size_t> FMDIndexIterator::pop() {
    // Drop the character that this added to the pattern.
    pattern.erase(pattern.size() - 1, 1);

    // Grab the stack frame to return
    std::pair<FMDPosition, size_t> toReturn(stack.back());

    // Drop it from the stack.
    stack.pop_back();

    // Return it
    return toReturn;
}
