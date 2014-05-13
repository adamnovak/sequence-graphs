#include <string>
#include <vector>

#include "BWT.h"
#include "SampledSuffixArray.h"

#include "TextPosition.hpp"
#include "FMDIndexIterator.hpp"
#include "RangeVector.hpp"
#include "Mapping.hpp"
#include "MapAttemptResult.hpp"



/**
 * A class that encapsulates access to an FMDIndex, consisting of the underlying
 * indexing-library-specific index structure and some auxilliary structures
 * (cointig names and lengths) that we store along with it.
 *
 * Wraps underlying indexing library operations and provides an interface with
 * bi-directional search and suffix tree iteration.
 *
 *
 */
class FMDIndex {
public:
    /**
     * Load an FMD and metadata from the given basename.
     */
    FMDIndex(std::string basename);
    
    /***************************************************************************
     * Metadata Access Functions
     **************************************************************************/
    
    /**
     * Get the contig number (not the text number) from a (text, offset) pair.
     */
    size_t getContigNumber(TextPosition base) const;
    
    /**
     * Get the strand from a (text, offset) pair: either 0 or 1.
     */
    bool getStrand(TextPosition base) const;
    
    /**
     * Given a pair_type representing a base, and a vector of contig lengths by
     * number, determine the specified base's offset from the left of its
     * contig, 1-based.
     */
    size_t getOffset(TextPosition base) const;
    
    /**
     * Get a unique string name for a position.
     */
    std::string getName(TextPosition base) const;
    
    /**
     * Get the total length of all contigs, on both strands.
     */
    int64_t getTotalLength() const;
    
    /** 
     * Get the total number of contigs in the index.
     */
    size_t getContigs() const;
    
    /**
     * Get the name of the contig at the given index.
     */
    const std::string& getContigName(size_t index) const;
    
    /**
     * Get the length of the contig at the given index.
     */
    size_t getContigLength(size_t index) const;
    
    /***************************************************************************
     * Search Functions
     **************************************************************************/
     
    /**
     * Get an FMDPosition covering the whole SA.
     */
    FMDPosition getSAPosition() const;
    
    /**
     * Get an FMDPosition for the part of the BWT for things starting with the
     * given character.
     */
    FMDPosition getCharPosition(char c) const;
     
    /**
     * Extend a search by a character, either backward or forward. Ranges are in
     * BWT coordinates.
     */
    FMDPosition extend(FMDPosition range, char c, bool backward) const;
    
    /***************************************************************************
     * Location/Sampled Suffix Array Functions
     **************************************************************************/
 
    /**
     * Find the (text, offset) position for an index in the BWT.
     */
    TextPosition locate(int64_t index) const;
    
    /***************************************************************************
     * Retrieval Functions
     **************************************************************************/
    
    /**
     * Get the character at the given index in the BWT.
     */
    char display(int64_t index) const;
    
    /***************************************************************************
     * Mapping Functions
     **************************************************************************/
      
    /**
     * Attempt to map each base in the query string to a (text, position) pair.
     * The vector returned will have one entry for each character in the
     * selected range.
     * 
     * Optionally a start and length for the region to map can be specified. The
     * whole string will be used as context, but only that region will actually
     * be mapped. A length of -1 means to use the entire string after the start,
     * and is the default.
     */
    std::vector<Mapping> map(const std::string& query, int start = 0,
        int length = -1) const;
      
    /**
     * Try RIGHT-mapping each base in the query to one of the ranges represented
     * by the range vector. The range vector is in BWT space, and has a 1 in the
     * first position in each range, and is 0 everywhere else. So rank(k) on it
     * gives the number of the range containing position k, and we can easily
     * check if both the start and end of our (backwards) search interval are in
     * the same range.
     *
     * TODO: Unify semantics!
     *
     * The range starting points must be such that the ranges described are "bi-
     * ranges": each range has its reverse-complement range also present.
     *
     * Returns a vector of range numbers for left-mapping each base, or -1 if
     * the base did not map to a range.
     */
    std::vector<int64_t> map(const RangeVector& ranges,
        const std::string& query, int start = 0, int length = -1) const;
        
    /***************************************************************************
     * Iteration Functions
     **************************************************************************/
     
    /**
     * We have an iterator typedef, so we can get an iterator over the suffix
     * tree easily.
     */
    typedef FMDIndexIterator iterator;
    /**
     * const_iterator is the same as the normal iterator, since it would be
     * silly to try and modify the suffix tree we're iterating over.
     */
    typedef FMDIndexIterator const_iterator;
    
    /**
     * Get an iterator pointing to the first range in the suffix tree, iterating
     * down to the given depth. If reportDeadEnds is set, will yield shorter-
     * than-depth contexts that happen before end of texts.
     */
    iterator begin(size_t depth, bool reportDeadEnds = false) const;
     
    /**
     * Get a 1-past-the-end sentinel iterator for the given depth.
     * reportDeadEnds should match whatever was specified on the corresponding
     * begin.
     */
    iterator end(size_t depth, bool reportDeadEnds = false) const;
    
protected:
    
    /**
     * Holds the names of all the contigs.
     */
    std::vector<std::string> names;
    
    /**
     * Holds the lengths of all the contigs, in the same order.
     */
    std::vector<size_t> lengths;
    
    /**
     * Holds the actual underlying index.
     */
    BWT bwt;
    
    /**
     * Holds the sampled suffix array we use for locate queries.
     */
    SampledSuffixArray suffixArray;
    
    /**
     * Try left-mapping the given index in the given string, starting from
     * scratch. Start a backwards search at that index in the string and extend
     * left until we map to exactly one or zero places. Returns true or false
     * depending on whether we map, an FMDPosition (in BWT coordinates) that, if
     * nonempty, can be extended right to try and map the next base to the
     * right, and the number of characters in the pattern used to make that
     * FMDPosition.
     *
     * If the mapping succeeded, the FMDPosition returned has one thing in it,
     * which is the mapping upstream context.
     *
     * Index must be a valid character position in the string.
     */
    MapAttemptResult mapPosition(const std::string& pattern,
        size_t index) const;
      
    /**
     * Try RIGHT-mapping the given index in the given string to a unique forward-
     * strand range according to the bit vector of range start points, starting
     * from scratch. Start a backwards search at that index in the string and
     * extend left until we map to exactly one or zero ranges. Returns true or
     * false depending on whether we map, an FMDPosition (in BWT coordinates)
     * that, if nonempty, can be extended right to try and map the next base to
     * the right, and the number of characters in the pattern used to make that
     * FMDPosition.
     *
     * The range starting points must be such that the ranges described are "bi-
     * ranges": each range has its reverse-complement range also present.
     *
     * If the mapping succeeded, the FMDPosition returned is completely
     * contained within one range, which is the range to which the base has been
     * mapped.
     *
     * Index must be a valid character position in the string.
     */
    MapAttemptResult mapPosition(const RangeVector& ranges, 
      const std::string& pattern, size_t index) const;

};
