#ifndef FMDINDEX_HPP
#define FMDINDEX_HPP

#include <string>
#include <vector>
#include <stdint.h>
#include <map>
#include <unordered_map>
#include <mutex>

#include "BWT.h"
#include "SampledSuffixArray.h"
#include "SuffixArray.h"
#include "MarkovModel.hpp"

#include "TextPosition.hpp"
#include "FMDIndexIterator.hpp"
#include "GenericBitVector.hpp"
#include "Mapping.hpp"
#include "MapAttemptResult.hpp"
#include "MismatchResultSet.hpp"
#include "LCPArray.hpp"

// State that the test cases class exists, even though we can't see it.
class FMDIndexTests;

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

    // Let the test code test our protected methods.
    friend FMDIndexTests;

public:
    /**
     * Load an FMD and metadata from the given basename. Optionally, specify a
     * complete suffix array that the index can use. The index takes ownership
     * of that suffix array, and will free it on destruction. Also optionally
     * takes a MarkovModel to allow enforcing a minimum coding cost for mapping.
     * The index takes ownership of the Markov model too.
     */
    FMDIndex(std::string basename, SuffixArray* fullSuffixArray = NULL, 
        MarkovModel* markovModel = NULL);
    
    ~FMDIndex();
    
    /***************************************************************************
     * Index Component Update Functions
     **************************************************************************/
    
    /** 
     * Use the given Markov model to model the probabilities of query sequences.
     * Takes ownership of the model (which may be NULL);
     */
    void setMarkovModel(MarkovModel* model);
    
    
    /***************************************************************************
     * Metadata Access Functions
     **************************************************************************/
    
    // TODO: TextPositions can already do like all of these themsleves.
    
    /**
     * Get the contig number (not the text number) from a (text, offset) pair.
     */
    size_t getContigNumber(TextPosition base) const;
    
    /**
     * Get the strand from a (text, offset) pair: either 0 or 1.
     */
    bool getStrand(TextPosition base) const;
    
    /**
     * Given a TextPosition representing a base, determine the specified base's
     * offset from the left of its contig, 1-based. Not to be confuised with the
     * position's internal offset along its text, from base.getOffset().
     */
    size_t getContigOffset(TextPosition base) const;
    
    /**
     * Get a unique string name for a position.
     */
    std::string getName(TextPosition base) const;
    
    /**
     * Given a TextPosition representing a base as (text, offset), get the index
     * of that base on that contig out of all bases on all contigs.
     */
    size_t getBaseID(TextPosition base) const;
    
    /**
     * Convert a (contig, base, face) to a TextPosition.
     */
    TextPosition getTextPosition(
        std::pair<std::pair<size_t, size_t>, bool> base) const;
    
    /**
     * Get the total length of all contigs, on both strands.
     */
    int64_t getTotalLength() const;
    
    /**
     * Get the total length of the BWT, which ought to be the number of texts
     * (i.e. twice the number of contigs), plus the number of characters of
     * actual sequence.
     */
    int64_t getBWTLength() const;
    
    /** 
     * Get the total number of contigs in the index.
     */
    size_t getNumberOfContigs() const;
    
    /**
     * Get the name of the contig at the given index. This is the sequence name
     * from the FASTA header that the sequence which contained the contig
     * originally had, for an index made with FMDIndexBuilder.
     */
    const std::string& getContigName(size_t index) const;
    
    /**
     * Get the start position of the contig at the given index. This is the
     * start in the original FASTA sequence from which the contig was obtained,
     * 0-based, for an index made with FMDIndexBuilder.
     */
    size_t getContigStart(size_t index) const;
    
    /**
     * Get the length of the contig at the given index.
     */
    size_t getContigLength(size_t index) const;
    
    /**
     * Get the number of the genome to which the contig at the given index
     * belongs.
     */
    size_t getContigGenome(size_t index) const;
    
    /**
     * Get the total number of genomes in the index.
     */
    size_t getNumberOfGenomes() const;
    
    /**
     * Get the range of contig numbers [start, end) that belong to the given
     * genome. All contigs in a genome appear one after the other.
     */
    std::pair<size_t, size_t> getGenomeContigs(size_t genome) const;
    
    /**
     * Get the minimum length that the genome can be while holding all its
     * contigs. TODO: Assumes 1-scaffold genomes
     */
    size_t getGenomeLength(size_t genome);
    
    /**
     * Return whether the given BWT position is in the given genome.
     */
    bool isInGenome(int64_t bwtIndex, size_t genome) const;
    
    /**
     * Get the mask for the positions in the given genome.
     */
    const GenericBitVector& getGenomeMask(size_t genome) const;
    
    /***************************************************************************
     * Search Functions
     **************************************************************************/
     
    /**
     * Get an FMDPosition covering the whole BWT.
     */
    FMDPosition getCoveringPosition() const;
    
    /**
     * Get an FMDPosition for the part of the BWT for things starting with the
     * given character.
     */
    FMDPosition getCharPosition(char c) const;
     
    /**
     * Extend a search by a character, either backward or forward.
     */
    FMDPosition extend(FMDPosition range, char c, bool backward) const;
    
    /**
     * Extend a search by a character, either backward or forward, in an
     * optimized way. Modifies the range to be extended, replacing it with the
     * extended version.
     */
    void extendFast(FMDPosition& range, char c, bool backward) const;
    
    /**
     * Extend a search by a character on the left in a backwards-search-only
     * extension. Compatible with retract. Modifies the range to be extended,
     * replacing it with the extended version.
     */
    void extendLeftOnly(FMDPosition& range, char c) const;
    
    /**
     * Given a search range, retract characters on the right until it
     * corresponds to a search pattern of the given length. Modifies the range
     * in place.
     *
     * Only updates the interval used for backward search, so after this method
     * is used, only extendLeftOnly can be used to extend the search.
     *
     * May not be called on an interval with no results in it.
     */
    void retractRightOnly(FMDPosition& range, size_t newPatternLength) const;
    
    /**
     * Retract on the right to the parent suffix tree node. Returns the new
     * pattern length represented by the search query. Guaranteed to produce
     * more results, unless you are at the root already.
     */
    size_t retractRightOnly(FMDPosition& range) const;
    
    /**
     * Select all the occurrences of the given pattern, using FMD backwards
     * search.
     */
    FMDPosition count(std::string pattern) const;
    
    /***************************************************************************
     * Longest Common Prefix (LCP) functions
     **************************************************************************/
    
    // TODO: Privatize or elimintate these? They're super small.
    
    /**
     * Get the LCP value at a certain position. This is the length of the prefix
     * that the suffix corresponding to this BWT position shares with the
     * previous one.
     */
    size_t getLCP(size_t index) const;
    
    /**
     * Get the previous smaller value in the LCP array before the given index,
     * or 0 if no such value exists (or the index is out of range).
     */
    size_t getLCPPSV(size_t index) const;
    
    /**
     * Get the index of the next smaller value in the LCP array after the given
     * index, or 1 past the end of the BWT if no such value exists (or if the
     * index is out of range).
     */
    size_t getLCPNSV(size_t index) const;
    
    
    /***************************************************************************
     * Location/Sampled Suffix Array Functions
     **************************************************************************/
 
    /**
     * Find the (text, offset) position for an index in the BWT.
     */
    TextPosition locate(int64_t index) const;
    
    // Unfortunately, unlocate cannot be efficiently implemented with
    // libsuffixtools's SampledSuffixArray.
    
    /**
     * Find the endpoint of the given contig in the BWT.
     */
    int64_t getContigEndIndex(size_t contig) const;
    
    /***************************************************************************
     * Retrieval Functions
     **************************************************************************/
    
    /**
     * Get the character at the given index in the BWT.
     */
    char display(int64_t index) const;
    
    /**
     * Get the character at a given offset into the given contig. Offset is
     * 1-based. TODO: this is not very efficient at all.
     */
    char display(size_t offset, size_t contig) const;
    
    /**
     * Get the character at a given TextPosition. Efficient if characters are
     * accessed with locality.
     */
    char displayCached(const TextPosition& position) const;  
    
    /**
     * Get the character in the first column of the given row in the BWT matrix.
     */
    char displayFirst(int64_t index) const;
    
    /**
     * Extract and return the forward strand of the given contig.
     */
    std::string displayContig(size_t index) const;
    
    /**
     * Memoized interface to displayContig that pulls out texts when needed,
     * using a volatile cache. Thread safe.
     */
    const std::string& displayContigCached(size_t index) const;
    
    /**
     * Remove any entries for this object in the global thread-local contig
     * cache. This prevents a new object at the same memory address from picking
     * up cache entries from an old, destructed object at that address.
     *
     * A hack since the destructor can't go talk to all the threads.
     *
     * Should be called whenever a thread is done with an index, and might use
     * another. Automatically called by the destructor.
     */
    void clearThreadLocalCache() const;
    
    /**
     * Given an index in the BWT, do an LF-mapping to get where the character
     * that appears in the first column at this index shows up in the last
     * column.
     */
    int64_t getLF(int64_t index) const;
    
    /***************************************************************************
     * Mapping Functions
     **************************************************************************/
      
    /**
     * Attempt to RIGHT-map each base in the query string to a (text, position)
     * pair. The vector returned will have one entry for each character in the
     * selected range. Things that right-map to forward strands are on even
     * texts.
     * 
     * If a mask is non-NULL, only positions in the index with a 1 in the mask
     * will be counted for mapping purposes.
     *
     */
    std::vector<Mapping> mapRight(const std::string& query, 
        const GenericBitVector* mask, int minContext = 0) const;
        
    /**
     * RIGHT-map to a specific genome, or to all genomes if genome is -1. Same
     * semantics as the function above.
     */
    std::vector<Mapping> mapRight(const std::string& query, int64_t genome = -1, 
        int minContext = 0) const;
     
    /**
     * LEFT-map to a specific genome, or to all genomes if genome is -1. Same
     * semantics as the function above. Things that left-map to forward strands
     * are on even texts.
     */    
    std::vector<Mapping> mapLeft(const std::string& query,
        int64_t genome = -1, int minContext = 0) const;
    
    /**
     * Both left- and right-map the given string to the given genome (or all
     * genomes if genome is -1). Things that map to forward strands are on even
     * texts.
     */
    std::vector<Mapping> mapBoth(const std::string& query, int64_t genome = -1, 
        int minContext = 0) const;
        
      
    /**
     * Try RIGHT-mapping each base in the query to one of the ranges represented
     * by the range vector. The range vector is in BWT space, and has a 1 in the
     * first position in each range, and is 0 everywhere else. So rank(k) on it
     * gives the one-based number of the range containing position k, and we can
     * easily check if both the start and end of our (backwards) search interval
     * are in the same range.
     *
     * The range starting points must be such that the ranges described are "bi-
     * ranges": each range has its reverse-complement range also present.
     *
     * If a mask is non-NULL, only positions in the index with a 1 in the mask
     * will be counted for mapping purposes.
     *
     * Returns a vector of one-based range numbers for left-mapping each base,
     * or -1 if the base did not map to a range.
     */
    std::vector<int64_t> mapRight(const GenericBitVector& ranges,
        const std::string& query, const GenericBitVector* mask,
        int minContext = 0) const;
        
    /**
     * RIGHT-map to ranges using contexts from a specific genome, or all genomes
     * if genome is -1. Same semantics as the function above.
     */
    std::vector<int64_t> mapRight(const GenericBitVector& ranges,
        const std::string& query, int64_t genome = -1,
        int minContext = 0) const;
      
    /**
     * Exact left-map without ranges.
     */
    std::vector<Mapping> map(const std::string& query, 
        const GenericBitVector* mask = NULL, int minContext = 0, int start = 0, 
        int length = -1) const; 

    /**
     * Exact left-map with ranges ranges by original restart-based algorithm.
     */
    std::vector<std::pair<int64_t,size_t>> map(const GenericBitVector& ranges,
        const std::string& query, const GenericBitVector* mask, 
        int minContext = 0, int addContext = 0, double multContext = 0, 
        double minCodingCost = 0, int start = 0, int length = -1) const;
        
    /**
     * Exact left-map with ranges to a genome by original restart-based
     * algorithm.
     *
     * minContext is the minuimum number of bases of context needed to map,
     * including the base being mapped.
     *
     * addContext is the number of bases required after uniqueness in order to 
     * map.
     *
     * multContext is the minimum fraction of the context length at which the
     * context became unique that is required in order to map. So if it is 2 and
     * you were unique at 70 bases, you would need at least 140 bases to map.
     *
     * minCodingCost is the minimum coding cost under this FMIndex's Markov
     * model for potential input sequences required for a context to map a base,
     * in bits.
     */
    std::vector<std::pair<int64_t,size_t>> map(const GenericBitVector& ranges,
        const std::string& query, int64_t genome = -1, int minContext = 0, 
        int addContext = 0, double multContext = 0, double minCodingCost = 0,
        int start = 0, int length = -1) const;

    /**
     * Given a left mapping and a right mapping for a base, disambiguate them to
     * produce one mapping. Things that consistently left and right-mapped to a
     * forward strand will be on an even text.
     *
     * TODO: Should this be Mapping::disambiguate()?
     */
    Mapping disambiguate(const Mapping& left, const Mapping& right) const;
    
    /***************************************************************************
     * Mismatch
     **************************************************************************/
    
    /**
     * Given a query string, a number of mismatches, and a mask of elligible
     * positions (which may be null), find all the strings in the index that are
     * within the specified number of mismatches of the query.
     *
     * Returns a MismatchResultSet which contains ranges for each string
     * within the specified number of mismatches that has any results.
     *
     * ranges is a bit vector marking the ranges that belong to each position,
     * so that we can check for uniqueness and set the is_mapped flag
     * appropriately.
     *
     * pattern gives the query string to be searched for.
     *
     * z_max gives the maximum number of mismatches.
     *
     * mask indicates which BWT positions should be included by setting their
     * bits to 1. If it is null, all BWT positions are included.
     */
    MismatchResultSet mismatchCount(const GenericBitVector& ranges,
        const std::string& pattern, size_t z_max, 
        const GenericBitVector* mask = NULL) const;
    
    /**
     * Given a set of mismatch search results, extend each result in the set,
     * throwing out those which do not exist in the reference.
     * 
     * z_max is the maximum number of mismatches allowed
     * 
     * if startExtension is true, we only extend by the correct base.
     * 
     * if finishExtension is true we only extend by the incorrect bases.
     * 
     * These two options are used for right extension of left-context results,
     * respectively left extension of right-context results
     * 
     * If both bools are false, then we extend with everything. This is what
     * happens when left-extending a set of left context results.
     * 
     * For speed, this implementation does not sort search results.
     */ 
    MisMatchAttemptResults misMatchExtend(MisMatchAttemptResults& prevMisMatches,
        char c, bool backward, size_t z_max, const GenericBitVector* mask,
        bool startExtension = false, bool finishExtension = false) const;
        
    /**
     * An old implementation of the method described above, sorting the results
     * in order of number of mismatches. This will allow a speed-up if ever we
     * want to find the match with the fewest number of mismatches.
     */
    MisMatchAttemptResults sortedMisMatchExtend(MisMatchAttemptResults& prevMisMatches,
            char c, bool backward, size_t z_max, const GenericBitVector* mask) const;
        
    /**
     * A submethod of sortedMisMatchExtend to check for existence and then sort
     * extension results
     */
    void processMisMatchPositions(
        MisMatchAttemptResults& nextMisMatches,
        std::vector<std::pair<FMDPosition,size_t>>& waitingMatches,
        std::vector<std::pair<FMDPosition,size_t>>& waitingMisMatches,
        const GenericBitVector* mask) const;
                
    /**
     * Implementing mismatch search for Left-Right exact contexts. addContext is
     * minimum additional context after uniqueness required to map.
     *
     * keepIntermediates can be set to false to force a restart at every base,
     * so minContext will be accurate; there's no way to calculate it if we're
     * allowed to extend from a previous result.
     *
     * Does right mapping.
     *
     * Returns mappings to ranges, with their TextPositions set to TextPositions
     * that are guaranteed to be in those ranges.
     */
    std::vector<Mapping> misMatchMap(const GenericBitVector& ranges,
        const std::string& query, const GenericBitVector* mask, 
        int minContext = 0, int addContext = 0, double multContext = 0,
        size_t z_max = 0, bool keepIntermediates = true, int start = 0,
        int length = -1) const;
        
    std::vector<Mapping> misMatchMap(const GenericBitVector& ranges, 
        const std::string& query, int64_t genome = -1, int minContext = 0, 
        int addContext = 0, double multContext = 0, size_t z_max = 0,
        bool keepIntermediates = true, int start = 0, int length = -1) const;
        
    // We have to pass back the extra context after uniqueness through a
    // pointer, because this design is not really suitable for all the extra
    // options and constraints we want to be able to tack on.
    MisMatchAttemptResults misMatchMapPosition(const GenericBitVector& ranges, 
        const std::string& pattern, size_t index, size_t minContext, 
        size_t addContext, double multContext, 
        int64_t* extraContext, size_t z_max, 
        const GenericBitVector* mask = NULL, bool maxContext=true) const;
        
    // Here is a slightly better designed overload of the above which leaves min
    // context enforcement up to its caller, and can be run in two modes to get
    // min and max context so that that is actually practical. Now can also
    // consider mismatches on the base being mapped, to produce a properly
    // extensible intermediate result.
    MisMatchAttemptResults misMatchMapPosition(const GenericBitVector& ranges, 
        const std::string& pattern, size_t index, size_t z_max, bool maxContext,
        bool allowFirstMismatch, const GenericBitVector* mask = NULL) const;
    
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
     * Holds the sequence names of all the contigs.
     */
    std::vector<std::string> names;
    
    /**
     * Holds the starts of all the contigs, in the same order.
     */
    std::vector<size_t> starts;
    
    /**
     * Holds the lengths of all the contigs, in the same order.
     */
    std::vector<size_t> lengths;
    
    /**
     * Holds the partial sums of contig lengths, or the ID of the first base in
     * each contig.
     */
    std::vector<size_t> cumulativeLengths;
    
    /**
     * Holds the genome index to which each contig belongs.
     */
    std::vector<size_t> genomeAssignments;
    
    /**
     * Holds the index in the BWT of the last base in each contig.
     */
    std::vector<int64_t> endIndices;
    
    // TODO: Change all these indexed-by-contig things into one vector of Contig
    // objects.
    
    /**
     * Holds, for each genome, a [start, end) range of contig numbers that
     * belong to it.
     */
    std::vector<std::pair<size_t, size_t>> genomeRanges;
    
    /**
     * Holds the bit vector masks for the BWT positions belonging to each
     * genome. Note that we can't get genome by contig or contigs for genome.
     * Only genome by BWT position.
     *
     * If we had C++11 we would use a vector, but we don't and since
     * GenericBitVector is not copy constructable/assignable we can't put it in
     * a vector. So we put pointers in a vector.
     * TODO: We have C++11 now. Fix this.
     */
    std::vector<GenericBitVector*> genomeMasks;
    
    /**
     * Holds the actual underlying index.
     */
    BWT bwt;
    
    /**
     * Holds the sampled suffix array we use for locate queries.
     */
    SampledSuffixArray suffixArray;
    
    /**
     * Holds a full suffix array pointer that we can use for locate queries
     * instead. Owned by this object, if not null.
     */
    SuffixArray* fullSuffixArray;
    
    /**
     * Holds a MarkovModel for calculating the coding costs of search strings.
     */
    MarkovModel* markovModel;
    
    /**
     * Holds the longest common prefix array.
     */
    LCPArray lcpArray;
    
    /**
     * Holds a cache of contig strings we have had to reconstruct for mapping on
     * credit. TODO: throw these out eventually somehow.
     */
    mutable std::map<size_t, std::string> contigCache;
    
    /**
     * Holds a mutex on contigCache so we don't break it trying to save contigs
     * on multiple threads.
     */
    mutable std::mutex contigCacheMutex;
    
    /**
     * A thread-local mirror of contigCache, so each thread can look for a local
     * copy before locking the global cache. Holds iterators to entries in the
     * global cache.
     *
     * TODO: how to evict things from deallocated FMDIndexes?
     */
    static thread_local std::map<const FMDIndex*, std::map<size_t,
        std::map<size_t, std::string>::iterator>> threadContigCache;
    
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
     *
     * If a mask is specified, only positions in the index with a 1 in the mask
     * will be counted for mapping purposes.
     */
    MapAttemptResult mapPosition(const std::string& pattern,
        size_t index, const GenericBitVector* mask = NULL) const;
      
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
     *
     * If a mask is specified, only positions in the index with a 1 in the mask
     * will be counted for mapping purposes.
     */
    creditMapAttemptResult CmapPosition(const GenericBitVector& ranges, 
        const std::string& pattern, size_t index, 
        const GenericBitVector* mask = NULL) const;
        
    MapAttemptResult mapPosition(const GenericBitVector& ranges, 
        const std::string& pattern, size_t index, 
        const GenericBitVector* mask = NULL) const;
        
private:
    
    // No copy constructor.
    FMDIndex(const FMDIndex& other);
    // No assignment operator.
    FMDIndex& operator=(const FMDIndex& other);

};

#endif
