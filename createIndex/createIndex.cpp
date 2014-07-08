#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <set>
#include <algorithm>
#include <utility>
#include <ctime>
#include <csignal>
#include <iterator>
#include <cstdint> 
#include <sys/resource.h>


#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// Grab pinchesAndCacti dependency.
#include <stPinchGraphs.h>

// Grab all the libFMD stuff.
#include <FMDIndexBuilder.hpp>
#include <BitVector.hpp>
#include <FMDIndexIterator.hpp>
#include <BitVector.hpp>
#include <TextPosition.hpp>
#include <util.hpp>
#include <Mapping.hpp>
#include <SmallSide.hpp>
#include <Log.hpp>

// Grab timers from libsuffixtools
#include <Timer.h>


#include "IDSource.hpp"
#include "ConcurrentQueue.hpp"
#include "OverlapMergeScheme.hpp"
#include "MappingMergeScheme.hpp"
#include "MergeApplier.hpp"


// TODO: replace with cppunit!

// Define a read to use for testing. 200 bp perfect match.
const std::string TEST_READ(std::string("GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGC") +
    "TGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGT" +
    "TGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGG" +
    "CCGGCGTATTAGAAGCGCG");
    
// 200 bp 5 changes
const std::string TEST_READ2(std::string("GACGGGACTCGCCGCCGCACAGCCGGGGTTCCCGC") +
    "TGGCGCAATTGAAAACTTTCGTCGACAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGT" +
    "TGGGGCAGTGACCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCCGAGAAAATGTCGATCGCCATTATGG" +
    "CCGGCGTATTCGAAGCGCG");

// 200 bp 10 changes
const std::string TEST_READ3(std::string("GACGGGAATCGCCGACGCCCAGCCGGGTTTCCGC") +
    "TGGCGCAATTTGAAAACTTTCGTCGATCAGGAATTTGCCCAAACAAAACATGTCCTGCATGGCGTTAGTTTGT" +
    "TGGGGCAGTGCCCGGTAGCATCAACGCTGCGCTGATTTGCCGTGGCGTAGAAAATGTCGATCGCCATTATGG" +
    "CCGGCGTATTCGAAGCGCG");
    
// 200bp unrelated
// >hg19_dna range=chr21:33031597-33031797 5'pad=0 3'pad=0 strand=+ repeatMasking=none
const std::string TEST_READ4(std::string("GCATCCATCTTGGGGCGTCCCAATTGCTGAGTAACAAATGAGACGCTGTG") +
    "GCCAAACTCAGTCATAACTAATGACATTTCTAGACAAAGTGACTTCAGAT" +
    "TTTCAAAGCGTACCCTGTTTACATCATTTTGCCAATTTCGCGTACTGCAA" +
    "CCGGCGGGCCACGCCCCCGTGAAAAGAAGGTTGTTTTCTCCACATTTCGG" +
    "G");

    
// How many times to try mapping this read?
const int TEST_ITERATIONS = 1000;

/**
 * Log current memory usage at INFO level.
 */
void logMemory() {
    
    // We have to interrogate /proc/self/status
    
    std::ifstream statusStream("/proc/self/status");
    
    std::string line;
    while(std::getline(statusStream, line)) {
        if(line.size() >= 2 && line[0] == 'V' && line[1] == 'm') {
            // This is a virtual memory line. Log it.
            Log::output() << line << std::endl;
        }
    }
    
    

}

/**
 * Signal handler to exit with exit(), preserving gprof profiling data if
 * applicable.
 */
void exitOnSignal(int signalNumber) {
    // Log the signal.
    Log::info() << "Exiting on signal " << signalNumber << std::endl;
    
    // Call exit and report the signal.
    exit(signalNumber);
}

/**
 * Look at the text of a base to see what arrow it needs to
 * have. Can also pass a base "strand" number.
 */
std::string 
getArrow(
    size_t textNumber
) {
    if(textNumber % 2 == 0) {
        // It's a left, so it gets a normal arrow
        return "normal";
    } else {
        // It's a right, so it gets that square thing.
        return "crow";
    }
}

/**
 * Start a new index in the given directory (by replacing it), and index the
 * given FASTAs for the bottom level FMD index. Optionally takes a suffix array
 * sample rate to use. Returns the basename of the FMD index that gets created.
 */
FMDIndex*
buildIndex(
    std::string indexDirectory,
    std::vector<std::string> fastas,
    int sampleRate = 128
) {

    // Make sure an empty indexDirectory exists.
    if(boost::filesystem::exists(indexDirectory)) {
        // Get rid of it if it exists already.
        boost::filesystem::remove_all(indexDirectory);
    }
    // Create it again.
    boost::filesystem::create_directory(indexDirectory);
    
    // Build the bottom-level index.
    
    // Work out its basename
    std::string basename(indexDirectory + "/index.basename");

    // Make a new builder
    FMDIndexBuilder builder(basename, sampleRate);
    for(std::vector<std::string>::iterator i = fastas.begin(); i < fastas.end();
        ++i) {
        
        Log::info() << "Adding FASTA " << *i << std::endl;
        
        // Add each FASTA file to the index.
        builder.add(*i);
    }
    
    Log::info() << "Finishing index..." << std::endl;    
    
    // Return the built index.
    return builder.build();
}

/**
 * Make a thread set with one thread representing each contig in the index.
 */
stPinchThreadSet* 
makeThreadSet(
    const FMDIndex& index
) {
    // Make a ThreadSet with one thread per contig.
    // Construct the thread set first.
    stPinchThreadSet* threadSet = stPinchThreadSet_construct();
    
    for(size_t i = 0; i < index.getNumberOfContigs(); i++) {
        // Add a thread for this contig number. The thread starts at base 1 and
        // has the appropriate length.
        stPinchThreadSet_addThread(threadSet, i, 1, index.getContigLength(i));
    }
    
    return threadSet;
}

/**
 * Create a new thread set from the given FMDIndex, and merge it down by the
 * overlap merging scheme, in parallel. Returns the pinched thread set.
 * 
 * If a context is specified, will not merge on fewer than that many bases of
 * context on a side, whether there is a unique mapping or not.
 */
stPinchThreadSet*
mergeOverlap(
    const FMDIndex& index,
    size_t context = 0
) {

    Log::info() << "Creating initial pinch thread set" << std::endl;
    
    // Make a thread set from our index.
    stPinchThreadSet* threadSet = makeThreadSet(index);
    
    // Make the merge scheme we want to use
    OverlapMergeScheme scheme(index, context);

    // Set it running and grab the queue where its results come out.
    ConcurrentQueue<Merge>& queue = scheme.run();
    
    // Make a merge applier to apply all those merges, and plug it in.
    MergeApplier applier(index, queue, threadSet);
    
    // Wait for these things to be done.
    scheme.join();
    applier.join();
    
    // Write a report before joining trivial boundaries.
    Log::output() << "Before joining boundaries:" << std::endl;
    Log::output() << "Pinch Blocks: " <<
        stPinchThreadSet_getTotalBlockNumber(threadSet) << std::endl;
    logMemory();
    
    // Now GC the boundaries in the pinch set
    Log::info() << "Joining trivial boundaries..." << std::endl;
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    
    // Write a similar report afterwards.
    Log::output() << "After joining boundaries:" << std::endl;
    Log::output() << "Pinch Blocks: " <<
        stPinchThreadSet_getTotalBlockNumber(threadSet) << std::endl;
    logMemory();
    
    // Now our thread set has been pinched. Return it.
    return threadSet;

}

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
) {
    
    Log::debug() << "Canonicalizing " << contigNumber << ":" << offset << 
        "." << strand << std::endl;
        
    // Now we need to look up what the pinch set says is the canonical
    // position for this base, and what orientation it should be in. All the
    // bases in this range have the same context and should thus all be
    // pointing to the canonical base's replacement.
    
    // Get the segment (using the 1-based position).
    stPinchSegment* segment = stPinchThreadSet_getSegment(threadSet, 
        contigNumber, offset);
        
    if(segment == NULL) {
        throw std::runtime_error("Found position in null segment!");
    }
        
    // How is it oriented in its block?
    bool segmentOrientation = stPinchSegment_getBlockOrientation(segment);
    // How far into the segment are we? 0-based, from subtracting 1-based
    // positions.
    size_t segmentOffset = offset - stPinchSegment_getStart(segment);
        
    // Get the first segment in the segment's block, or just this segment if
    // it isn't in a block.
    stPinchSegment* firstSegment = segment;
    
    // Get the block that that segment is in
    stPinchBlock* block = stPinchSegment_getBlock(segment);
    if(block != NULL) {
        // Put the first segment in the block as the canonical segment.
        firstSegment = stPinchBlock_getFirst(block);
    }
    
    // Work out what the official contig number for this block is (the name
    // of the first segment).
    size_t canonicalContig = stPinchSegment_getName(firstSegment);
    // How should it be oriented?
    bool canonicalOrientation = stPinchSegment_getBlockOrientation(
        firstSegment);
    
    // We need to calculate an offset into the canonical segment. If the
    // segments are pinched together backwards, this will count in opposite
    // directions on the two segments, so we'll need to flip the within-sement
    // offset around.
    size_t canonicalSegmentOffset = segmentOffset;
    if(segmentOrientation != canonicalOrientation) {
        // We really want this many bases in from the end of the segment, not
        // out from the start of the segment. Keep it 0-based.
        canonicalSegmentOffset = stPinchSegment_getLength(firstSegment) - 
            canonicalSegmentOffset - 1;
    }
    
    Log::debug() << "Canonicalized segment offset " << segmentOffset << 
        " to " << canonicalSegmentOffset << std::endl;
    
    // What's the offset into the canonical contig? 1-based because we add a
    // 0-based offset to a 1-based position.
    size_t canonicalOffset = stPinchSegment_getStart(firstSegment) + 
        canonicalSegmentOffset;
    
    Log::debug() << "Canonicalized contig " << contigNumber << " offset " <<
        offset << " to contig " << canonicalContig << " offset " << 
        canonicalOffset << std::endl;
    
    // Return all three values, and be sad about not having real tuples. What
    // orientation should we use?  Well, we have the canonical position's
    // orientation within the block, our position's orientation within the
    // block, and the orientation that this context attaches to the position in.
    // Flipping any of those will flip the orientation in which we need to map,
    // so we need to xor them all together, which for bools is done with !=.
    Log::debug() << "Canonical orientation: " << canonicalOrientation << 
            std::endl;
    Log::debug() << "Segment orientation: " << segmentOrientation << 
            std::endl;
    Log::debug() << "Strand: " << strand << std::endl;

    if(canonicalOffset <= 0 || 
        canonicalOffset > stPinchThread_getLength(
        stPinchSegment_getThread(firstSegment))) {
        
        // The answer is supposed to be a 1-based offset. Make sure that is
        // true.
        throw std::runtime_error("Canonical offset " + 
            std::to_string(canonicalOffset) + 
            " out of range for 1-based position");
    }

    return std::make_pair(std::make_pair(canonicalContig, canonicalOffset),
        canonicalOrientation != segmentOrientation != strand);
}

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
) {
    
    Log::debug() << "Canonicalizing 0-based " << base << std::endl;
    
    // What contig corresponds to that text?
    size_t contigNumber = index.getContigNumber(base);
    // And what strand corresponds to that text? This tells us what
    // orientation we're actually looking at the base in.
    bool strand = (bool) index.getStrand(base);
    // And what base position is that from the front of the contig? This is
    // 1-based.
    size_t offset = index.getOffset(base);
    
    if(offset <= 0 || offset > index.getContigLength(contigNumber)) {
        // Complain that we got an out of bounds offset from this thing.
        throw std::runtime_error("Tried to canonicalize text" +
            std::to_string(base.getText()) + " offset " + 
            std::to_string(base.getOffset()) + " which is out of bounds");
    }
    
    // Canonicalize that pinch thread set position.
    return canonicalize(threadSet, contigNumber, offset, strand);
    
}

/**
 * Given a pinch thread set and a filename, write the degree of each pinch block
 * or segment without a block (always degree 2, except at the ends) to the file,
 * one per line.
 *
 * TODO: Maybe restructure as one pass through all blocks, and one scan through
 * threads ignoring any blocks? Then we could eliminate the seen set for blocks.
 */
void
writeDegrees(
    stPinchThreadSet* threadSet, 
    std::string filename
) {

    Log::info() << "Saving pinch graph degrees to " << filename << std::endl;

    // Open up the file to write.
    std::ofstream degrees(filename.c_str());

    // We're going go through each thread one segment at a time, to make sure we
    // catch all the single segments.
    // TODO: When we get c++11, make this an unordered_set
    std::set<stPinchBlock*> seen;
    
    // Make an iterator over pinch threads
    stPinchThreadSetIt threadIterator = stPinchThreadSet_getIt(threadSet);
    
    // And a pointer to hold the current thread
    stPinchThread* thread;
    while((thread = stPinchThreadSetIt_getNext(&threadIterator)) != NULL) {
        // For each pinch thread
        
        // Go through all its pinch segments in order. There's no iterator so we
        // have to keep looking 3'
        
        // Get the first segment in the thread.
        stPinchSegment* segment = stPinchThread_getFirst(thread);
        while(segment != NULL) {
        
            stPinchBlock* block;
            if((block = stPinchSegment_getBlock(segment)) != NULL) {
                // Look at its block if it has one
                
                if(seen.count(block) == 0) {
                    // This is a new block we haven't reported the degree of
                    // yet.
                    
                    // Report the degree.
                    degrees << stPinchBlock_getDegree(block) << std::endl;
                    
                    // Remember we saw it.
                    seen.insert(block);
                }
            } else {
                // This is a bare segment; report degree depending on what it's
                // attached to on each end.
                
                // How many neighbors do we have. TODO: track this as we scan,
                // resulting in fewer queries.
                int degree = (int)(stPinchSegment_get3Prime(segment) != NULL) +
                    (int)(stPinchSegment_get5Prime(segment) != NULL);
                    
                degrees << degree << std::endl;
                
            }
            
            // Jump to the next 3' segment. This needs to return NULL if we go
            // off the end.
            segment = stPinchSegment_get3Prime(segment);
        }
        
    }
    
    // We wrote out all the degrees
    degrees.flush();
    degrees.close();

}

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
    std::string filename) 
{

    // We're going to lay out the segments in our root sequence with some
    // locality. We're going to scan through all threads, and put each new block
    // as it is encountered. So we need a set of the encountered blocks.
    // TODO: When we get c++11, make this an unordered_set
    std::set<stPinchBlock*> seen;
    
    // Open up the file to write.
    std::ofstream c2h(filename.c_str());
    
    // First, make a hacked-up consensus reference sequence to be the root of
    // the tree. It just has all the pinch blocks in some order.
    
    // Write the root sequence line. It is a bottom sequence since it has bottom
    // segments.
    c2h << "s\t'rootSeq'\t'rootSeq'\t1" << std::endl;
    
    // Keep track of the total root sequence space already used
    size_t nextBlockStart = 0;
    
    // Make an iterator over pinch threads
    stPinchThreadSetIt threadIterator = stPinchThreadSet_getIt(threadSet);
    
    // And a pointer to hold the current thread
    stPinchThread* thread;
    while((thread = stPinchThreadSetIt_getNext(&threadIterator)) != NULL) {
        // For each pinch thread
        
        // Go through all its pinch segments in order. There's no iterator so we
        // have to keep looking 3'
        
        // Get the first segment in the thread.
        stPinchSegment* segment = stPinchThread_getFirst(thread);
        while(segment != NULL) {
        
            stPinchBlock* block;
            if((block = stPinchSegment_getBlock(segment)) != NULL) {
                // Look at its block if it has one
                
                if(seen.count(block) == 0) {
                    // This is a new block we haven't allocated a range for yet.
                    
                    // For each block, put a bottom segment. Just name the
                    // segment after the block's address. We need block name
                    // (number), start, and length.
                    c2h << "a\t" << (uintptr_t)block << "\t" << 
                        nextBlockStart << "\t" << 
                        stPinchBlock_getLength(block) << std::endl;
                        
                    // Advance nextBlockStart.
                    nextBlockStart += stPinchBlock_getLength(block);
                    
                    // Record that we have seen this block now.
                    seen.insert(block);
                }
            }
            
            // Jump to the next 3' segment. This needs to return NULL if we go
            // off the end.
            segment = stPinchSegment_get3Prime(segment);
        }
        
    }
    
    // Keep a mapping from scaffold name to event name. Event name will be
    // either the contig scaffold name if all the contigs are from a single
    // scaffold, or "genome-<number>" if there are multiple scaffolds involved.
    std::map<std::string, std::string> eventNames;
    
    // Keep track of the original source sequence
    std::string sourceSequence;
    
    // Go through contigs in order. The index spec requires them to be grouped
    // by original source sequence.
    for(size_t contig = 0; contig < index.getNumberOfContigs(); contig++) {
        // Grab the sequence name that the contig is on
        std::string contigName = index.getContigName(contig);
        
        if(eventNames.count(contigName) == 0) {
            // We need to figure out the event name for this contig.
            
            // Start by naming the event after the contig.
            eventNames[contigName] = contigName;
            
            // What genome does it belong to?
            size_t genomeNumber = index.getContigGenome(contig);
            
            // What contigs are in that genome?
            auto genomeRange = index.getGenomeContigs(genomeNumber);
            
            for(size_t i = genomeRange.first; i < genomeRange.second; i++) {
                // For every contig in that genome
                if(index.getContigName(i) != contigName) {
                    // One of them doesn't match this contig's name; they are
                    // not all from the same scaffold. Re-name the event with a
                    // new generic name.
                    
                    eventNames[contigName] = "genome-" + 
                        std::to_string(index.getContigGenome(contig));
                        
                    // Now we're done
                    break;
                }
            }
        }
        
        if(contigName != sourceSequence || contig == 0) {
            // This is a new scaffold, not the same as the one we were on last.
            
            // Start a new sequence with a sequence line. The sequence is a top
            // sequence, since it is only connected up.
            c2h << "s\t'" << eventNames[contigName] << "'\t'" << contigName <<
                "'\t0" << std::endl;
            
            // Remember that we are on this sequence
            sourceSequence = contigName;
        } else {
            // This is the same sequence as before, but we need an unaligned
            // segment to cover the distance from the last contig to the start
            // of this one.
            
            // What's 1 base after the end of the last contig? Leave as 0-based.
            size_t prevContigEnd = index.getContigStart(contig - 1) + 
                index.getContigLength(contig - 1);
            
            // Unaligned top segments are "a <start> <length>" only. Leave as
            // 0-based.
            c2h << "a\t" << prevContigEnd << "\t" << 
                index.getContigStart(contig) - prevContigEnd << std::endl;
            
        }
        
        // Regardless of whether we changed sequence or not, we need an
        // alignment segment for every segment in this thread.
        
        // Each segment has to account for the offset of the contig on the
        // sequence.
        size_t contigStart = index.getContigStart(contig);
        
        // So first get the thread by name.
        stPinchThread* thread = stPinchThreadSet_getThread(threadSet, contig);
        
        // Go through all its pinch segments in order. There's no iterator so we
        // have to keep looking 3'
        
        // Get the first segment in the thread.
        stPinchSegment* segment = stPinchThread_getFirst(thread);
        while(segment != NULL) {
            
            // Start where the next segment starts. Convert from 1-based pinch
            // segments to 0-based HAL.
            size_t segmentStart = contigStart +
                stPinchSegment_getStart(segment) - 1;
            
            if(stPinchSegment_getBlock(segment) != NULL) {
                // It actually aligned
            
                // Write a top segment mapping to the segment named after the
                // address of the block this segment belongs to.
                c2h << "a\t" << segmentStart << "\t" << 
                    stPinchSegment_getLength(segment) << "\t" << 
                    (uintptr_t)stPinchSegment_getBlock(segment) << "\t" << 
                    stPinchSegment_getBlockOrientation(segment) << std::endl;
            } else {
                // Write a segment for the unaligned sequence.
                c2h << "a\t" << segmentStart << "\t" << 
                    stPinchSegment_getLength(segment) << std::endl;
            }
            
            // Jump to the next 3' segment. This needs to return NULL if we go
            // off the end.
            segment = stPinchSegment_get3Prime(segment);
        }
        
    }
    
    // Close up the finished file
    c2h.close();
    
    // Return the total length of the made-up root sequence. Later we will
    // probably need a FASTA with this many Ns in order to turn this c2h into a
    // HAL file.
    return nextBlockStart;
}

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
) {

    // Open the FASTA to write
    std::ofstream fasta(filename.c_str());
    
    Log::info() << "Generating " << rootBases << 
        " bases of root node sequence." << std::endl;
    
    // First we put the right number of Ns in a sequence named "rootSeq"
    fasta << ">rootSeq" << std::endl;
    for(size_t i = 0; i < rootBases; i++) {
        fasta << "N";
        // Entire sequence must be on one line.
    }
    
    for(std::vector<std::string>::iterator i = inputFastas.begin(); 
        i != inputFastas.end(); ++i) {
        
        Log::info() << "Copying over " << *i << std::endl;
        
        // Then we just copy all the other FASTAs in order.
        
        // Open the input FASTA
        std::ifstream inputFasta((*i).c_str());
        
        // Read it line by line. See <http://stackoverflow.com/a/7868998/402891>
        std::string line;
        while(std::getline(inputFasta, line)) {
            // For each line (without trailing newline)
        
            if(line.size() == 0) {
                // Drop blank lines
                continue;
            }
            
            if(line[0] == '>') {
                // Make sure there are newlines before and after header lines.
                fasta << std::endl << line << std::endl;
            } else {
                // Don't put any newlines.
                fasta << line;
            }
        }
            
        // Close up this input file and move to the next one.
        inputFasta.close();
    }
    
    // Insert a linebreak at the end of the file.
    fasta << std::endl;
    
    // Now we're done.
    fasta.close();

}

/**
 * Canonicalize each contigous run of positions mapping to the same canonical
 * base and face.
 *
 * Takes a pinched thread set, and the index on which it is defined.
 *
 * If a mask bit vector is specified, only positions which have a 1 in the mask
 * will be able to break up ranges.
 *
 * Returns a BitVector marking each such range with a 1 at the start, and a
 * vector of canonicalized positions. Each anonicalized position is a contig
 * number, a base number, and a face flag.
 *
 * All BWT positions must be represented in the pinch set.
 *
 * Scans through the entire BWT.
 */
std::pair<BitVector*, std::vector<std::pair<std::pair<size_t, size_t>, bool> > > 
identifyMergedRuns(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index,
    const BitVector* mask = NULL
) {
    
    // We need an iterator over the mask, if applicable.
    BitVectorIterator* maskIterator = (mask != NULL ? 
        new BitVectorIterator(*mask) : NULL);
    
    // We need to make bit vector denoting ranges, which we encode with this
    // encoder, which has 32 byte blocks.
    BitVectorEncoder encoder(32);
    
    // We also need to make a vector of canonical positions.
    std::vector<std::pair<std::pair<size_t, size_t>, bool> > mappings;
    
    Log::info() << "Building merged run index by scan..." << std::endl;
    
    // Do the thing where we locate each base and, when the canonical position
    // changes, add a 1 to start a new range and add a mapping.

    // Keep track of the ID and relative orientation for the last position we
    // canonicalized.
    // TODO: typedef this! It is getting silly.
    std::pair<std::pair<size_t, size_t>, bool> lastCanonicalized;
    
    for(int64_t j = index.getNumberOfContigs() * 2; j < index.getBWTLength();
        j++) {
        
        // For each base in the BWT (skipping over the 2*contigs stop
        // characters)... 
        
        // Note: stop characters aren't at the front in the last column, only
        // the first.
        
        if(maskIterator != NULL && !maskIterator->isSet(j)) {
            // This position is masked out. We don't allow it to break up
            // ranges, so no range needs to start here. Skip it and pretend it
            // doesn't exist.
            continue;
        }
        
        // Where is it located?
        TextPosition base = index.locate(j);
        
        // Canonicalize it. The second field here will be the relative
        // orientation and determine the Side.
        std::pair<std::pair<size_t, size_t>, bool> canonicalized = 
            canonicalize(index, threadSet, base);
            
        if(j == 0 || canonicalized != lastCanonicalized) {
            // We need to start a new range here, because this BWT base maps to
            // a different position than the last one.
            
            // Say this range is going to belong to the canonical base.
            mappings.push_back(canonicalized);
            
            if(j != index.getNumberOfContigs() * 2) {
                // Record a 1 in the vector at the start of every range except
                // the first. The first needs no 1 before it so it will be rank
                // 0 (and match up with mapping 0), and it's OK not to split it
                // off from the stop characters since they can't ever be
                // searched.
                encoder.addBit(j);
                Log::debug() << "Set bit " << j << std::endl;
            }
            
            
            // Remember what canonical base and face we're doing for this range.
            lastCanonicalized = canonicalized;
        }
        // Otherwise we had the same canonical base, so we want this in the same
        // range we already started.
    }
            
    // Set a bit after the end of the last range (i.e. at the end of the BWT).
    encoder.addBit(index.getBWTLength());
    
    // Finish the vector encoder into a vector of the right length, leaving room
    // for that trailing bit. Make sure to flush first.
    encoder.flush();
    BitVector* bitVector = new BitVector(encoder,
        index.getBWTLength() + 1);
    
    // If we made a mask iterator, get rid of it.
    if(maskIterator != NULL) {
        delete maskIterator;
    }
    
    // Return the bit vector and the canonicalized base vector
    return std::make_pair(bitVector, mappings);
}

/**
 * Make the range vector and list of matching Sides for the hierarchy level
 * implied by the given thread set in the given index. Gets IDs for created
 * positions from the given source.
 * 
 * Manages this by scanning the BWT from left to right, seeing what cannonical
 * position and orientation each base belongs to, and putting 1s in the range
 * vector every time a new one starts.
 *
 * Don't forget to delete the bit vector when done!
 */
std::pair<BitVector*, std::vector<SmallSide> > 
makeLevelIndexScanning(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    IDSource<long long int>& source
) {
    
    // First, canonicalize everything, yielding a bitvector of ranges and a
    // vector of canonicalized positions.
    auto mergedRuns = identifyMergedRuns(threadSet, index);
    
    // We also need to make a vector of SmallSides, which are the things that
    // get matched to by the corresponding ranges in the bit vector.
    std::vector<SmallSide> mappings;
    
    // We also need this map of position IDs (long long ints) by canonical
    // contig name (a size_t) and base index (also a size_t). Relative
    // orientation is always forward.
    std::map<std::pair<size_t, size_t>, long long int>
        idReservations;

    
    size_t runNumber = 0;

    Log::info() << "Building mapping data structure..." << std::endl;
    
    for(auto canonicalized : mergedRuns.second) {
        
        // For each canonical 1-based ((contig, base), face) corresponding to a
        // merged range...
        
        // What position is that?
        long long int positionCoordinate;
        if(idReservations.count(canonicalized.first) > 0) {
            // Load the previously chosen ID
            positionCoordinate = 
                idReservations[canonicalized.first];
        } else {
            // Allocate and remember a new ID.
            positionCoordinate = 
                idReservations[canonicalized.first] = 
                source.next();
        }
            
        // Say this range is going to belong to the ID we just looked up, on
        // the appropriate face.
        mappings.push_back(
            SmallSide(positionCoordinate, canonicalized.second));
    }
            
    // Return the bit vector and the Side vector
    return std::make_pair(mergedRuns.first, mappings);
}

/**
 * Save both parts of the given level index to files in the given directory,
 * which must not yet exist. Does not delete the bit vector from the level
 * index, so it can be reused.
 */
void saveLevelIndex(
    std::pair<BitVector*, std::vector<SmallSide> > levelIndex,
    std::string directory
) {
    
    Log::info() << "Saving index to disk..." << std::endl;
    
    // Make the directory
    boost::filesystem::create_directory(directory);
    
    // Save the bit vector to a file.
    std::ofstream vectorStream((directory + "/vector.bin").c_str(), 
        std::ios::binary);
    levelIndex.first->writeTo(vectorStream);
    vectorStream.close();
    
    // Open a file to save all the sides to
    std::ofstream sideStream((directory + "/mappings.bin").c_str(), 
        std::ios::binary);
    for(std::vector<SmallSide>::iterator i = levelIndex.second.begin(); 
        i != levelIndex.second.end(); ++i) {
        
        // For each side, write it. We can't write them all at once since we
        // might have more than an int's worth of bytes.
        (*i).write(sideStream);
    
    }
    sideStream.close();
}

/**
 * Create a new thread set from the given FMDIndex, and merge it down by the
 * greedy merging scheme, in parallel. Returns the pinched thread set.
 *
 * The greedy merging scheme is to take one genome, map the next genome to it,
 * merge at the mappings, compute the upper-level index, map the next genome to
 * that, merge at the mappings, compute the upper level index, and so on.
 * 
 * If a context is specified, will not merge on fewer than that many bases of
 * context on a side, whether there is a unique mapping or not.
 */
stPinchThreadSet*
mergeGreedy(
    const FMDIndex& index,
    size_t context = 0
) {

    Log::info() << "Creating initial pinch thread set" << std::endl;
    
    // Make a thread set from our index.
    stPinchThreadSet* threadSet = makeThreadSet(index);
    
    if(index.getNumberOfGenomes() == 0) {
        // Make sure we have at least 1 genome.
        throw std::runtime_error("Can't merge 0 genomes greedily!");
    }
    
    // Keep around a bit vector of all the positions that are in. This will
    // start with the very first genome, which we know exists.
    const BitVector* includedPositions = &index.getGenomeMask(0);
    
    // Canonicalize everything, yielding a bitvector (pointer) of ranges and a
    // vector of canonicalized positions.
    auto mergedRuns = identifyMergedRuns(threadSet, index);
    
    for(size_t genome = 1; genome < index.getNumberOfGenomes(); genome++) {
        // For each genome that we have to merge in...
        
        // Make the merge scheme we want to use. We choose a mapping-to-second-
        // level-based merge scheme, to which we need to feed the details of the
        // second level (range vector, representative positions to merge into)
        // and also the bitmask of what bottom-level things to count. We also
        // need to tell it what genome to map the contigs of.
        MappingMergeScheme scheme(index, *mergedRuns.first, mergedRuns.second,
            *includedPositions, genome, context);

        // Set it running and grab the queue where its results come out.
        ConcurrentQueue<Merge>& queue = scheme.run();
        
        // Make a merge applier to apply all those merges, and plug it in.
        MergeApplier applier(index, queue, threadSet);
        
        // Wait for these things to be done.
        scheme.join();
        applier.join();
        
        // Compute and log some statistics about the coverage of the alignments
        // used in the merge.
        
        // How many bases were aligned?
        auto lock = queue.lock();
        size_t basesAligned = queue.getThroughput(lock);
        
        // How many bases were alignable?
        size_t basesAlignable = 0;
        for(size_t i = index.getGenomeContigs(genome).first; 
            i < index.getGenomeContigs(genome).second; i++) {
                
            // Sum up the bases in every contig in the genome. Anything outside
            // a contig (i.e. Ns) could never be aligned adn should not have
            // been counted.
            basesAlignable += index.getContigLength(i);
        }
        
        // Print out the coverage obtained from this.
        Log::output() << "Coverage from alignment of genome " << genome << 
            ": " << basesAligned << " / " << basesAlignable << " = " <<
            ((double)basesAligned) / basesAlignable << std::endl;
        
        // Join any trivial boundaries.
        stPinchThreadSet_joinTrivialBoundaries(threadSet);
        
        // Merge the new genome into includedPositions, replacing the old
        // bitvector.
        BitVector* newIncludedPositions = includedPositions->createUnion(
            index.getGenomeMask(genome));
        if(genome > 1) {
            // If we already alocated a new BitVector that wasn't the one that
            // came when we loaded in the genomes, we need to delete it.
            delete includedPositions;
        }
        includedPositions = newIncludedPositions;
        
        // Delete the old merged runs bit vector and recalculate merged runs on
        // the newly updated thread set. Make sure to specify the new mask of
        // included positions, so that masked-out positions don't break ranges
        // that would otherwise be merged.
        delete mergedRuns.first;
        mergedRuns = identifyMergedRuns(threadSet, index, includedPositions);
        
    }
    
    // Delete the final merged run vector
    delete mergedRuns.first;
    
    if(index.getNumberOfGenomes() > 1) {
        // And, if we had to make any additional included position BitVectors,
        // get the last one of those too.
        delete includedPositions;
    }
    
    // Let the caller rebuild that from this pinch graph which we
    // return.
    return threadSet;
}

/**
 * Map a read repeatedly to the bottom level.
 */
void
testBottomMapping(
    const FMDIndex& index
) {
    // Start the timer
    clock_t start = clock();
    for(int i = 0; i < TEST_ITERATIONS; i++) {
        // Map repeatedly
        index.map(TEST_READ);
    }
    // Stop the timer
    clock_t end = clock();
    
    // Work out how many milliseconds each call took
    double msPerCall = ((double)(end - start)) / 
        (CLOCKS_PER_SEC / 1000.0) / TEST_ITERATIONS;
        
    Log::output() << "Mapping to bottom level: " << msPerCall << 
        " ms per call" << std::endl;
}

/**
 * Map a read repeatedly to the merged level.
 */
void
testMergedMapping(
    const FMDIndex& index, const BitVector* ranges
) {
    // Start the timer
    clock_t start = clock();
    for(int i = 0; i < TEST_ITERATIONS; i++) {
        // Map repeatedly
        index.map(*ranges, TEST_READ);
    }
    // Stop the timer
    clock_t end = clock();
    
    // Work out how many milliseconds each call took
    double msPerCall = ((double)(end - start)) / 
        (CLOCKS_PER_SEC / 1000.0) / TEST_ITERATIONS;
        
    Log::output() << "Mapping to merged level: " << msPerCall <<
        " ms per call" << std::endl;
}

/**
 * createIndex: command-line tool to create a multi-level reference structure.
 */
int 
main(
    int argc, 
    char** argv
) {

    // Register ctrl+c handler. See
    // <http://www.yolinux.com/TUTORIALS/C++Signals.html>
    signal(SIGINT, exitOnSignal);

    // Parse options with boost::programOptions. See
    // <http://www.radmangames.com/programming/how-to-use-boost-program_options>

    std::string appDescription = 
        std::string("Create a reference hierarchy for mapping to FASTAs.\n") + 
        "Usage: createIndex <index directory> <fasta> [<fasta> [<fasta> ...]]";

    // Make an options description for our program's options.
    boost::program_options::options_description description("Options");
    // Add all the options
    description.add_options() 
        ("help", "Print help messages") 
        ("test", "Run a mapping speed test")
        ("noMerge", "Don't compute merged level, only make lowest-level index")
        ("scheme", boost::program_options::value<std::string>()
            ->default_value("overlap"),
            "Merging scheme (\"overlap\" or \"greedy\")")
        ("alignment", boost::program_options::value<std::string>(), 
            "File to save .c2h-format alignment in")
        ("alignmentFasta", boost::program_options::value<std::string>(), 
            "File in which to save FASTA records for building HAL from .c2h")
        ("degrees", boost::program_options::value<std::string>(), 
            "File in which to save degrees of pinch graph nodes")
        ("context", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum required context length to merge on")
        ("sampleRate", boost::program_options::value<unsigned int>()
            ->default_value(64), 
            "Set the suffix array sample rate to use")
        // These next two options should be ->required(), but that's not in the
        // Boost version I can convince our cluster admins to install. From now
        // on I shall work exclusively in Docker containers or something.
        ("indexDirectory", boost::program_options::value<std::string>(), 
            "Directory to make the index in; will be deleted and replaced!")
        ("fastas", boost::program_options::value<std::vector<std::string> >()
            ->multitoken(),
            "FASTA files to load");
        
    // And set up our positional arguments
    boost::program_options::positional_options_description positionals;
    // One index directory
    positionals.add("indexDirectory", 1);
    // And an unknown number of FASTAs
    positionals.add("fastas", -1);
    
    // Add a variables map to hold option variables.
    boost::program_options::variables_map options;
    
    try {
        // Parse options into the variable map, or throw an error if there's
        // something wring with them.
        boost::program_options::store(
            // Build the command line parser.
            boost::program_options::command_line_parser(argc, argv)
                .options(description)
                .positional(positionals)
                .run(),
            options);
        boost::program_options::notify(options);
            
        if(options.count("help")) {
            // The help option was given. Print program help.
            std::cout << appDescription << std::endl;
            std::cout << description << std::endl;
            
            // Don't do the actual program.
            return 0; 
        }
        
        if(!options.count("indexDirectory") || !options.count("fastas")) {
            throw boost::program_options::error("Missing important arguments!");
        }
            
    } catch(boost::program_options::error& error) {
        // Something is bad about our options. Complain on stderr
        std::cerr << "Option parsing error: " << error.what() << std::endl;
        std::cerr << std::endl; 
        // Talk about our app.
        std::cerr << appDescription << std::endl;
        // Show all the actually available options.
        std::cerr << description << std::endl; 
        
        // Stop the program.
        return -1; 
    }
    
    // If we get here, we have the right arguments. Parse them.
    
    // This holds the directory for the reference structure to build.
    std::string indexDirectory(options["indexDirectory"].as<std::string>());
    
    // This holds a list of FASTA filenames to load and index.
    std::vector<std::string> fastas(options["fastas"]
        .as<std::vector<std::string> >());
        
    // Dump options.
    Log::output() << "Options:" << std::endl;
    Log::output() << "Store index in: " << indexDirectory << std::endl;
    
    for(std::vector<std::string>::iterator i = fastas.begin();
        i != fastas.end(); ++i) {
        
        Log::output() << "Index file: " << *i << std::endl;
    }
    
    // Index the bottom-level FASTAs. Use the
    // sample rate the user specified.
    FMDIndex* indexPointer = buildIndex(indexDirectory, fastas,
        options["sampleRate"].as<unsigned int>());
        
    // Make a reference out of the index pointer because we're not letting it
    // out of our scope.
    FMDIndex& index = *indexPointer;
    
    // Log memory usage with no pinch graph stuff having yet happened.
    Log::output() << "Memory usage with no merging:" << std::endl;
    logMemory();
    
    if(options.count("noMerge")) {
        // Skip merging any of the higher levels.
        return 0;
    }
    
    // Grab the merging scheme we are going to use to merge.
    std::string mergeScheme = options["scheme"].as<std::string>();
    
    // Make a pointer to hold the threadset pointer it will create.
    stPinchThreadSet* threadSet;
    
    // We want to time the merge code.
    Timer* mergeTimer = new Timer("Merging");
    
    if(mergeScheme == "overlap") {
        // Make a thread set that's all merged, with the given minimum merge
        // context.
        threadSet = mergeOverlap(index, options["context"].as<size_t>());
    } else if(mergeScheme == "greedy") {
        // Use the greedy merge instead.
        threadSet = mergeGreedy(index, options["context"].as<size_t>());
    } else {
        // Complain that's not a real merge scheme. TODO: Can we make the
        // options parser parse an enum or something instead of this?
        throw std::runtime_error(mergeScheme + 
            " is not an implemented merge scheme.");
    }
    
    // Now the merge is done. Stop timing.
    delete mergeTimer;
        
    if(options.count("degrees")) {
        // Save a dump of pinch graph node degrees (for both blocks and bare
        // segments).
        writeDegrees(threadSet, options["degrees"].as<std::string>());
    }
    
    if(options.count("alignment")) {
        // Save the alignment defined by the pinched pinch graph to the file the
        // user specified. Save the number of bases of root sequence that were
        // used in the center of the star tree.
        size_t rootBases = writeAlignment(threadSet, index,
            options["alignment"].as<std::string>());
            
        if(options.count("alignmentFasta")) {
            // Also save a FASTA with the sequences necessary to generate a HAL
            // from the above.
            writeAlignmentFasta(fastas, rootBases,
                options["alignmentFasta"].as<std::string>());
        }
    }
    
    // Make an IDSource to produce IDs not already claimed by contigs.
    IDSource<long long int> source(index.getTotalLength());
    
    // This will hold the computed level index of the merged level.
    std::pair<BitVector*, std::vector<SmallSide> > levelIndex;
    
    // We also want to time the merged level index building code
    Timer* levelIndexTimer = new Timer("Level Index Construction");
    
    // Use a scanning strategy for indexing.
    levelIndex = makeLevelIndexScanning(threadSet, index, source);
    
    delete levelIndexTimer;
        
    // Write it out, deleting the bit vector in the process
    saveLevelIndex(levelIndex, indexDirectory + "/level1");
    
    // Clean up the thread set
    stPinchThreadSet_destruct(threadSet);
    
    // Run the speed tests if we want to
    if(options.count("test")) {
        Log::output() << "Running performance tests..." << std::endl;
        testBottomMapping(index);
        testMergedMapping(index, levelIndex.first);
    }
    
    // Get rid of the range vector
    delete levelIndex.first;
    
    // Get rid of the index itself. Invalidates the index reference.
    delete indexPointer;

    Log::output() << "Final memory usage:" << std::endl;
    logMemory();

    // Now we're done!
    return 0;
}
