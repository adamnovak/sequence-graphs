#ifndef PINCHGRAPHUTIL_HPP
#define PINCHGRAPHUTIL_HPP

// Grab pinchesAndCacti dependency.
#include <stPinchGraphs.h>

#include <utility>
#include <vector>
#include <string>

#include <FMDIndex.hpp>
#include <TextPosition.hpp>

/**
 * pinchGraphUtil.hpp: utility functions for pinch graph.
 */
 

/**
 * Make a thread set with one thread representing each contig in the index.
 */
stPinchThreadSet* 
makeThreadSet(
    const FMDIndex& index
);

/**
 * Turn the given (contig number, 1-based offset from start, orientation)
 * position into the same sort of structure for the canonical base that
 * represents all the bases it has been pinched with.
 */
std::pair<std::pair<size_t, size_t>, bool>
canonicalize(
    stPinchThreadSet* threadSet, 
    size_t contigNumber,
    size_t offset,
    bool strand
);

/**
 * Turn the given (text, 0-based offset offset) pair into a canonical (contig
 * number, 1-based offset from contig start, orientation), using the given
 * FMDIndex and the given thread set. The orientation is which face of the
 * canonical base this (text, offset) pair means.
 */
std::pair<std::pair<size_t, size_t>, bool>
canonicalize(
    const FMDIndex& index, 
    stPinchThreadSet* threadSet, 
    TextPosition base
);

/**
 * Write the given threadSet on the contigs in the given index out as a multiple
 * alignment. Produces a cactus2hal (.c2h) file as described in
 * <https://github.com/benedictpaten/cactus/blob/development/hal/impl/hal.c>/
 *
 * Such a file looks like this.
 * 
 * s    'event1'   'seq1'    1
 * a    1 10   10
 * s    'event1'    'seq2'  0
 * a    20  10  1    0
 *
 * These files describe trees of sequences, with an alignment of each child onto
 * its parent. Each sequence has two ways of dividing it up into "segments", or
 * blocks: one which applies for alignments mapping upwards, and one which
 * appplies for alignments mapping downwards. 
 * 
 * The s lines define sequences. Sequences can be bottom or not. Each sequence
 * is defined by an event name, a contig name/FASTA header, and a bottom-or-not
 * flag. After each s line there are a series of a lines, which define alignment
 * segments on that sequence.
 * 
 * "Bottom" alignment lines define "bottom" segments, which are mapped down;
 * they consist only of a name (really a number), a start, and a length. Bottom
 * sequences only have bottom segments.
 *
 * "Top" alignment lines define "top" segments, which have a start, a length, a
 * name (which is numerical) of a bottom alignment block which they are aligned
 * to, and a flag for the relative orientation of the alignment (0 for same
 * direction, right for different directions). If the segment is an insertion
 * (i.e. unaligned) it has only a start and a length. Top sequences only have
 * top segments.
 *
 * All segments in a sequence are sorted by position. Each sequence appears only
 * once. Every base in a sequence is part of some segment. All positions are
 * 0-based. Ancestors must come before descendants.
 *
 * Returns the total number of bases in the made-up root node used to tie the
 * actual sequences together.
 */
size_t
writeAlignment(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    std::string filename
);

/**
 * Write a FASTA file that goes with the .c2h file written by writeAlignment, so
 * that the halAppendCactusSubtree tool can turn both into a HAL file.
 *
 * Strips out newlines so halAppendCactusSubtree will be happy with the
 * resulting FASTA.
 *
 * TODO: Drop this and just use the HAL API directly
 */
void
writeAlignmentFasta(
    std::vector<std::string> inputFastas,
    size_t rootBases,
    std::string filename
);

#endif
