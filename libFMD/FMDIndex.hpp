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

#include "TextPosition.hpp"
#include "FMDIndexIterator.hpp"
#include "GenericBitVector.hpp"
#include "Mapping.hpp"
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
     * of that suffix array, and will free it on destruction.
     */
    FMDIndex(std::string basename, SuffixArray* fullSuffixArray = NULL);
    
    ~FMDIndex();
    
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
        
private:
    
    // No copy constructor.
    FMDIndex(const FMDIndex& other);
    // No assignment operator.
    FMDIndex& operator=(const FMDIndex& other);

};

#endif
