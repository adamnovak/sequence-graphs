#ifndef PINCHGRAPHUTIL_HPP
#define PINCHGRAPHUTIL_HPP

// Grab pinchesAndCacti dependency.
#include <stPinchGraphs.h>

#include <utility>
#include <vector>
#include <string>
#include <set>

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
 * Turn the given TextPosition into the TextPosition for the canonical base that
 * represents all the bases it has been pinched with.
 */
TextPosition
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
    const std::string& filename
);

/*
 * Write the alignment represented by this pinch graph and index as a star tree
 * alignment where the given refrence sequence is used as the root. Saves in c2h
 * format, as described above, to the given file.
 */
void
writeAlignmentWithReference(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    const std::string& filename, 
    const size_t referenceGenomeNumber
);

/**
 * Write the bottom segments for a single-thread scaffold.
 */
void 
writeBottomSegments(
    std::ofstream& c2h, 
    stPinchThread* thread
);

/**
 * Write the top segments for a single-thread scaffold.
 */
void 
writeTopSegments(
    std::ofstream& c2h, 
    stPinchThread* thread
);

/*
 * Write the alignment represented by this pinch graph to the given file.
 *
 * Returns the total number of bases in the made-up root node used to tie the
 * actual sequences together.
 *
 * threadNames gives the sequence name for each pinch graph thread, and
 * threadEvents gives the event name for each thread.
 *
 * If eventsToKeep is not null, only outputs events with names in that set (plus
 * the automatically added rootSeq).
 *
 */
size_t
writeAlignment(
    stPinchThreadSet* threadSet, 
    std::vector<std::string> threadNames,
    std::vector<std::string> threadEvents,
    const std::string& filename,
    const std::set<std::string>* eventsToKeep = nullptr
);

/*
 * Write the alignment represented by this pinch graph alone as a star where the
 * given refrence thread is used as the root. Saves in c2h format, as
 * described above, to the given file.
 *
 * threadNames gives the sequence name for each pinch graph thread, and
 * threadEvents gives the event name for each thread.
 *
 * TODO: Multiple reference threads.
 */
void
writeAlignmentWithReference(
    stPinchThreadSet* threadSet, 
    std::vector<std::string> threadNames,
    std::vector<std::string> threadEvents,
    const std::string& filename, 
    const size_t referenceThreadNumber
);

/**
 * Write a FASTA file that goes with the .c2h file written by writeAlignment, so
 * that the halAppendCactusSubtree tool can turn both into a HAL file.
 *
 * Strips out newlines so halAppendCactusSubtree will be happy with the
 * resulting FASTA.
 *
 * If rootBases < 0, don't actually generate a rootSeq record.
 *
 * TODO: Drop this and just use the HAL API directly
 *
 * TODO: Move to another file since this doesn't touch pinch graphs.
 */
void
writeAlignmentFasta(
    std::vector<std::string> inputFastas,
    int64_t rootBases,
    std::string filename
);

/**
 * Write the given pinch graph in (badly non-standard) LastGraph format,
 * intended for import by Bandage.
 *
 * Takes the thread set, and a file to save to.
 *
 * The format begins with a tab-separated header line.
 *
 * <node count>\t<sequences (always 0)>\t<kmer length used (always 0)>\t<unknown
 * (always 0)>
 *
 * It then has NODE stanzas for each segment/node in the graph.
 *
 * NODE\t<ID>\t<length>\t<total bp merged into node>
 * <forward sequence of combined kmers for node>
 * <reverse sequence of combined kmers for node>
 *
 * Note that if k were nonzero, these sequences wouldn't be expected to be
 * proper reverse complements of each other. Also note that they are always
 * output as all Ns at the moment.
 *
 * Finally, there are ARC lines, between node IDs. The IDs are negated to
 * indicate the "wrong" sides of the relevant nodes (left side as the first
 * node, or right side as the second).
 *
 * ARC\t<(-)node 1>\t<(-)node 2>
 *
 * TODO: needs a way to get the sequence bases.
 */
void
writeLastGraph(
    stPinchThreadSet* threadSet,
    const std::string& filename
);

#endif
