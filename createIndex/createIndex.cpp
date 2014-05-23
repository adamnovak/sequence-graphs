#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <set>
#include <algorithm>
#include <utility>
#include <ctime>
#include <sys/resource.h>


#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// Grab pinchesAndCacti dependency.
#include <stPinchGraphs.h>

// Grab all the libFMD stuff.
#include <FMDIndexBuilder.hpp>
#include <FMDIndex.hpp>
#include <FMDIndexIterator.hpp>
#include <RangeVector.hpp>
#include <TextPosition.hpp>
#include <util.hpp>
#include <Mapping.hpp>
#include <SmallSide.hpp>
#include <IDSource.hpp>
#include <Log.hpp>


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
    
    // We have to interogate /proc/self/status
    
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
 * nonsymmetric merging scheme to contexts of the given length (including the
 * base being matched itself). Returns the pinched thread set.
 *
 * If dumpFile is set, dumps debug graph data to that file.
 * If quiet is true, don't announce contexts.
 *
 * Note that due to the nature of this merging scheme, any two nodes that would
 * merge at a longer context length will also merge at a shorter context length,
 * so we can just directly calculate each upper level in turn.
 */
stPinchThreadSet*
mergeNonsymmetric(
    const FMDIndex& index,
    size_t contextLength,
    std::ostream* dumpFile = NULL,
    bool quiet = false
) {
    
    // Keep track of nodes we have already dumped.
    std::set<std::string> dumped;
    
    // Make a thread set from our index.
    stPinchThreadSet* threadSet = makeThreadSet(index);
    
    // To construct the non-symmetric merged graph with p context:
    // Traverse the suffix tree down to depth p + 1
    for(FMDIndex::iterator i = index.begin(contextLength); 
        i != index.end(contextLength); ++i) {
        // For each pair of suffix and position in the suffix tree
        
        // Unpack the iterator into pattern and FMDPosition at which it happens.
        std::string pattern = (*i).first;
        // This is in SA coordinates.
        FMDPosition range = (*i).second;
        
        
        if(!quiet) {
            // Dump the context and range.
            Log::debug() << pattern << " at " << range << std::endl;
        }
        
        if(range.getEndOffset() >= 1 || dumpFile != NULL) {
            // We only need to do any pinching if this context appears in two or
            // more places. And since a range with offset 0 has one thing in it,
            // we check to see if it's 1 or more.
            
            Log::debug() << "Locating " << range << std::endl;
            
            // For each base location, we need to work out the contig and base
            // and orientation, and pinch with the first.
            
            // Work out what text and base the first base is.
            TextPosition firstBase = index.locate(range.getForwardStart());
            
            Log::trace() << "First relative position: text " << 
                firstBase.getText() << " offset " << firstBase.getOffset() <<
                std::endl;
            
            // What contig corresponds to that text?
            size_t firstContigNumber = index.getContigNumber(firstBase);
            // And what strand corresponds to that text?
            size_t firstStrand = index.getStrand(firstBase);
            // And what base position is that?
            size_t firstOffset = index.getOffset(firstBase);
            
            // Grab the first pinch thread
            stPinchThread* firstThread = stPinchThreadSet_getThread(threadSet,
                firstContigNumber);
                
            if(firstThread == NULL) {
                throw std::runtime_error("First thread was NULL!");
            }
            
            if(dumpFile != NULL) {
                // Report the first position as existing.
                std::string firstName = index.getName(firstBase); 
                
                if(!dumped.count(firstName)) {
                    // Get its base
                    char baseChar = index.displayFirst(range.getForwardStart());
                    if(firstStrand) {
                        // Flip it around so we always see forward strand bases.
                        baseChar = complement(baseChar);
                    }
                
                    // Write a node for it
                    *dumpFile << firstName << "[shape=\"record\",label=\"{" << 
                        firstName << "|" << baseChar <<
                        "}\"];" << std::endl;
                    if(firstOffset > 1) {
                        // Link to previous Position
                        *dumpFile << index.getName(TextPosition(
                            // Hack to get the base actually before us on the
                            // contig.
                            firstContigNumber * 2, firstOffset - 2)) << 
                            " -> " << firstName << ";" << std::endl;
                    }
                    dumped.insert(firstName);
                }
            }
            
            for(int64_t j = 1; j < range.getEndOffset() + 1; j++) {
                // For each subsequent base
                
                // Locate the base
                TextPosition otherBase = index.locate(range.getForwardStart() +
                    j);
                
                Log::trace() << "Relative position: (" << 
                    otherBase.getText() << "," << otherBase.getOffset() << 
                    ")" << std::endl;   
                
                size_t otherContigNumber = index.getContigNumber(otherBase);
                size_t otherStrand = index.getStrand(otherBase);
                size_t otherOffset = index.getOffset(otherBase);
                
                // Grab the other pinch thread.
                stPinchThread* otherThread = stPinchThreadSet_getThread(
                    threadSet, otherContigNumber);
                    
                if(firstThread == NULL) {
                    throw std::runtime_error("Other thread was NULL!");
                }
            
                // What orientation should we use for the second strand, given
                // that we are pinching against the first strand in orientation
                // 1 (reverse)?
                bool orientation = firstStrand == otherStrand;
            
                // Pinch firstBase on firstNumber and otherBase on otherNumber
                // in the correct relative orientation.
                Log::trace() << "\tPinching #" << firstContigNumber << ":" <<
                    firstOffset << " strand " << firstStrand << " and #" << 
                    otherContigNumber << ":" << otherOffset << " strand " << 
                    otherStrand << " (orientation: " << orientation << ")" <<
                    std::endl;
                
                stPinchThread_pinch(firstThread, otherThread, firstOffset,
                    otherOffset, 1, orientation);
                    
                if(dumpFile != NULL) {
                    // Report the other position as existing.
                    std::string otherName = index.getName(otherBase); 
                    
                    if(!dumped.count(otherName)) {
                        // Get its base character
                        char baseChar = index.displayFirst(
                            range.getForwardStart() + j);
                            
                        if(otherStrand) {
                            // Flip it around so we always see forward strand
                            // bases.
                            baseChar = complement(baseChar);
                        }
                    
                        // Write a node for it
                        *dumpFile << otherName << 
                            "[shape=\"record\",label=\"{" << otherName << "|" <<
                            baseChar << "}\"];" << std::endl;
                        
                        if(otherOffset > 1) {
                            // Link previous position to us.
                            *dumpFile << index.getName(TextPosition(
                                // Hack to get the base actually before us on
                                // the contig.
                                otherContigNumber * 2, otherOffset - 2)) << 
                                " -> " << otherName << ";" << std::endl;
                        }
                        dumped.insert(otherName);
                    }
                }
                
            }
            
            if(!quiet) {
                // Say we merged some bases.
                Log::debug() << "Merged " << range.getEndOffset() + 1 <<  
                    " bases" << std::endl;
            }
            
        }
    }
    
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
    
    // Return the finished thread set
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

    return std::make_pair(std::make_pair(canonicalContig, canonicalOffset),
        canonicalOrientation != segmentOrientation != strand);
}

/**
 * Turn the given (text, offset) pair into a canonical (contig number, 1-based
 * offset from contig start, orientation), using the given FMDIndex and the
 * given thread set. The orientation is which face of the canonical base this
 * (text, offset) pair means.
 */
std::pair<std::pair<size_t, size_t>, bool>
canonicalize(
    const FMDIndex& index, 
    stPinchThreadSet* threadSet, 
    TextPosition base
) {
    
    // What contig corresponds to that text?
    size_t contigNumber = index.getContigNumber(base);
    // And what strand corresponds to that text? This tells us what
    // orientation we're actually looking at the base in.
    bool strand = (bool) index.getStrand(base);
    // And what base position is that from the front of the contig? This is
    // 1-based.
    size_t offset = index.getOffset(base);
    
    // Canonicalize that pinch thread set position.
    return canonicalize(threadSet, contigNumber, offset, strand);
    
}

/**
 * Write GraphViz edges for adjacencies between merged nodes, which are derived
 * by canonicalizing and looking up the IDs for successive positions, and
 * linking them together.
 *
 * Takes a pinch thread set that knows what all the contigs are and how they got
 * pinched, a map from canonical pinch thread bases to merged position IDs, and
 * the stream to send the GraphViz data to.
 */
void makeMergedAdjacencies(
    stPinchThreadSet* threadSet, 
    std::map<std::pair<size_t, size_t>, long long int> idReservations, 
    std::ofstream* dumpFile
) {

    // We need a way to keep track of what edges we have already made between
    // sides. So keep a set of merged base IDs and orientation flags. Make sure
    // the key pairs are sorted!
    typedef std::pair<long long int, bool> mergedSide;
    std::set<std::pair<mergedSide, mergedSide> > addedEdges;
    
    // Iterate through the pinch threads. We need to do a bit of an odd loop
    // since this is a bit of an odd iterator.
    stPinchThread* thread;
    for(stPinchThreadSetIt i = stPinchThreadSet_getIt(threadSet); 
        (thread = stPinchThreadSetIt_getNext(&i)) != NULL; ) {
        // For each thread...
        
        // Get its name
        size_t name = stPinchThread_getName(thread);
        for(size_t j = stPinchThread_getStart(thread); 
            j < stPinchThread_getStart(thread) + 
            stPinchThread_getLength(thread) - 1; j++) {
            
            // For each base in the thread that isn't the last...
            
            Log::debug() << "Linking merged base " << j << " on thread " <<
                name << std::endl;
            
            
            
            // Canonicalize its right side (strand 0 or false), which is the one
            // linked by this adjacency to the next base.
            std::pair<std::pair<size_t, size_t>, bool> canonicalBase =
                canonicalize(threadSet, name, j, false);
            
            // Also canonicalize the left side of the next base.
            std::pair<std::pair<size_t, size_t>, bool> canonicalNextBase =
                canonicalize(threadSet, name, j + 1, true);
                
            // What IDs and sides do these go to? Hold pairs of coordinate and
            // side.
            mergedSide side = 
                std::make_pair(idReservations[canonicalBase.first],
                canonicalBase.second);
            mergedSide nextSide = 
                std::make_pair(idReservations[canonicalNextBase.first],
                canonicalNextBase.second);
            
            
            if(nextSide < side) {
                // These are in the wrong order. All our edges need to go from
                // smallest to largest, so we don't repeat them.
                std::swap(side, nextSide);
            }
            
            // What edge would we make? Pair up the sides in order.
            std::pair<mergedSide, mergedSide> edge = std::make_pair(side,
                nextSide);
                
            if(addedEdges.count(edge) == 0) {
                // This is a new edge between merged sides. Print it out.
                *dumpFile << "M" << side.first << " -> " << "M" << 
                    nextSide.first << "[dir=both,arrowtail=" << 
                    getArrow(!side.second) << ",arrowhead=" << 
                    getArrow(!nextSide.second) << ",color=red];" << std::endl;
                    
                // Add the edge.
                addedEdges.insert(edge);
            }
        
        }  
        
    }

    
}

/**
 * Make the range vector and list of matching Sides for the hierarchy level
 * implied by the given thread set in the given index. Gets IDs for created
 * positions from the given source.
 * 
 * Manages this by traversing the suffix tree, and blocking out each range
 * (since we're sure a range of the correct length has all been merged into the
 * same base and face).
 *
 * If dumpFile is set, also writes a graphviz-format debug graph.
 * 
 * Don't forget to delete the bit vector when done!
 */
std::pair<RangeVector*, std::vector<SmallSide> > 
makeLevelIndex(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    size_t contextLength,
    IDSource<long long int>& source, 
    std::ofstream* dumpFile = NULL
) {
    
    // We need to make bit vector denoting ranges, which we encode with this
    // encoder, which has 32 byte blocks.
    RangeEncoder encoder(32);
    
    // We also need to make a vector of SmallSides, which are the things that
    // get matched to by the corresponding ranges in the bit vector.
    std::vector<SmallSide> mappings;
    
    // We also need this map of position IDs (long long ints) by canonical
    // contig name (a size_t) and base index (also a size_t)
    std::map<std::pair<size_t, size_t>, long long int>
        idReservations;

    Log::info() << "Building mapping data structure by tree traversal..." <<
        std::endl;
    
    for(FMDIndex::iterator i = index.begin(contextLength, true); 
        i != index.end(contextLength, true); ++i) {
        
        // Go through all the contexts, including shortened ones before end of
        // text.
        
        // Unpack the iterator into pattern and FMDPosition at which it happens,
        // which is in BWT coordinates.
        std::string context = (*i).first;
        FMDPosition range = (*i).second;
        
        Log::debug() << "===Context: " << context << " at range " << range <<
            "===" << std::endl;
        
        if(context.size() == contextLength) {
        
            // If it's a full-requested-length context (and therefore all the
            // bases in it have been pinched)...
            
            // Locate the first base in the context. It's already in SA
            // coordinates.
            TextPosition base = index.locate(range.getForwardStart());
            Log::debug() << "Text/Offset: (" << base.getText() << ", " << 
                base.getOffset() << ")" << std::endl;
            
            // Canonicalize it. The second field here will be the relative
            // orientation and determine the face.
            std::pair<std::pair<size_t, size_t>, bool> canonicalized = 
                canonicalize(index, threadSet, base);
            
            // Figure out what the ID of the canonical base is.
            long long int positionCoordinate;
            if(idReservations.count(canonicalized.first) > 0) {
                // Load the previously chosen ID
                positionCoordinate = idReservations[canonicalized.first];
            } else {
                // Allocate and remember a new ID.
                positionCoordinate = idReservations[canonicalized.first] = 
                    source.next();
            }
            
            // Say this range is going to belong to the ID we just looked up, on
            // the appropriate face.
            mappings.push_back(
                SmallSide(positionCoordinate, canonicalized.second));
            
            if(range.getForwardStart() != 0) {
                // Record a 1 in the vector at the start of every range except
                // the first, in BWT coordinates. The first needs no 1 before it
                // so it will be rank 0 (and match up with mapping 0), and it's
                // OK not to split it off from the stop characters since they
                // can't ever be searched.
                encoder.addBit(range.getForwardStart());
                Log::debug() << "Set bit " << range.getForwardStart() <<
                    std::endl;
            }
            
            // Put a 1 at the start of its interval, and add a mapping to the
            // given ID and Side.
            
            // TODO: Make edges by traversing the pinch graph.
                
        } else if(context.size() < contextLength) {
            // Otherwise if it's a shorter than requested context, it's
            // possible not everything has been merged. So do the old thing
            // where we locate each base and, when the canonical position
            // changes, add a 1 to start a new range and add a mapping.

            // Keep track of the ID and relative orientation for the last
            // position we canonicalized.
            std::pair<std::pair<size_t, size_t>, bool> lastCanonicalized;
            
            for(int64_t j = 0; j < range.getEndOffset() + 1; j++) {
                // For each base in the range...
                
                // Where is it located?
                TextPosition base = index.locate(range.getForwardStart() + j);
                
                // Canonicalize it. The second field here will be the relative
                // orientation and determine the Side.
                std::pair<std::pair<size_t, size_t>, bool> canonicalized = 
                    canonicalize(index, threadSet, base);
                    
                if(j == 0 || canonicalized != lastCanonicalized) {
                    // We need to start a new range here, because this BWT base
                    // maps to a different position than the last one.
                    
                    // What position is that? TODO: Unify with ID assignment
                    // above.
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
                    
                    // Say this range is going to belong to the ID we just
                    // looked up, on the appropriate face.
                    mappings.push_back(
                        SmallSide(positionCoordinate, canonicalized.second));
                    
                    // Add in the bit the same as above, only here we add the
                    // offset of i.
                    if(range.getForwardStart() + j != 0) {
                        // Record a 1 in the vector at the start of every range
                        // except the first, in BWT coordinates. The first needs
                        // no 1 before it so it will be rank 0 (and match up
                        // with mapping 0), and it's OK not to split it off from
                        // the stop characters since they can't ever be
                        // searched.
                        encoder.addBit(range.getForwardStart() + j);
                        Log::debug() << "Set bit " << 
                            range.getForwardStart() + j << std::endl;
                    }
                    
                    
                    // Remember what canonical base and face we're doing for
                    // this range.
                    lastCanonicalized = canonicalized;
                }
                // Otherwise we had the same canonical base, so we want this in
                // the same range we already started.
            }
        } else {
            // We got a context that's longer than we asked for?
            throw std::runtime_error("Iterated context longer than requested!");
        }
    }
    
    // Set a bit after the end of the last range.
    encoder.addBit(index.getTotalLength());
    
    // Finish the vector encoder into a vector of the right length.
    // This should always end in a 1!
    // Make sure to flush first.
    encoder.flush();
    RangeVector* bitVector = new RangeVector(encoder,
        index.getTotalLength() + 1);
    
    if(dumpFile != NULL) {
        // We need to add in the edges that connect merged positions together.
        // This requires walking all the contigs again, so we put it in its own
        // function.
        makeMergedAdjacencies(threadSet, idReservations, dumpFile);
    }
    
    // Return the bit vector and the Side vector
    return make_pair(bitVector, mappings);
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
 * If dumpFile is set, also writes a graphviz-format debug graph.
 * 
 * Don't forget to delete the bit vector when done!
 */
std::pair<RangeVector*, std::vector<SmallSide> > 
makeLevelIndexScanning(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    size_t contextLength,
    IDSource<long long int>& source, 
    std::ofstream* dumpFile = NULL
) {
    
    // We need to make bit vector denoting ranges, which we encode with this
    // encoder, which has 32 byte blocks.
    RangeEncoder encoder(32);
    
    // We also need to make a vector of SmallSides, which are the things that
    // get matched to by the corresponding ranges in the bit vector.
    std::vector<SmallSide> mappings;
    
    // We also need this map of position IDs (long long ints) by canonical
    // contig name (a size_t) and base index (also a size_t)
    std::map<std::pair<size_t, size_t>, long long int>
        idReservations;

    Log::info() << "Building mapping data structure by scan..." << std::endl;
    
    
    // Otherwise if it's a shorter than requested context, it's possible not
    // everything has been merged. So do the old thing where we locate each base
    // and, when the canonical position changes, add a 1 to start a new range
    // and add a mapping.

    // Keep track of the ID and relative orientation for the last position we
    // canonicalized.
    std::pair<std::pair<size_t, size_t>, bool> lastCanonicalized;
    
    for(int64_t j = index.getNumberOfContigs() * 2; j < index.getBWTLength();
        j++) {
        
        // For each base in the BWT (skipping over the 2*contigs stop
        // characters)...
        
        // Where is it located?
        TextPosition base = index.locate(j);
        
        // Canonicalize it. The second field here will be the relative
        // orientation and determine the Side.
        std::pair<std::pair<size_t, size_t>, bool> canonicalized = 
            canonicalize(index, threadSet, base);
            
        if(j == 0 || canonicalized != lastCanonicalized) {
            // We need to start a new range here, because this BWT base maps to
            // a different position than the last one.
            
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
            
    // Set a bit after the end of the last range.
    encoder.addBit(index.getTotalLength());
    
    // Finish the vector encoder into a vector of the right length.
    // This should always end in a 1!
    // Make sure to flush first.
    encoder.flush();
    RangeVector* bitVector = new RangeVector(encoder,
        index.getTotalLength() + 1);
    
    if(dumpFile != NULL) {
        // We need to add in the edges that connect merged positions together.
        // This requires walking all the contigs again, so we put it in its own
        // function.
        makeMergedAdjacencies(threadSet, idReservations, dumpFile);
    }
    
    // Return the bit vector and the Side vector
    return make_pair(bitVector, mappings);
}

/**
 * Save both parts of the given level index to files in the given directory,
 * which must not yet exist. Does not delete the bit vector from the level
 * index, so it can be reused.
 */
void saveLevelIndex(
    std::pair<RangeVector*, std::vector<SmallSide> > levelIndex,
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
    const FMDIndex& index, const RangeVector* ranges
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
        ("dump", "Dump GraphViz graphs")
        ("test", "Run a mapping speed test")
        ("quiet", "Don't print every context")
        ("noMerge", "Don't compute merged level, only make lowest-level index")
        ("scan", "Make merged mapping bit vector by scanning BWT")
        ("context", boost::program_options::value<unsigned int>()
            ->default_value(3), 
            "Set the context length to merge on")
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
        
    // This holds the length of context to use
    unsigned int contextLength = options["context"].as<unsigned int>();
    
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
    
    if(options.count("noMerge")) {
        // Skip merging any of the higher levels.
        return 0;
    }
    
    Log::info() << "Use " << contextLength << " bases of context." << std::endl;
    
    // Make an IDSource to produce IDs not already claimed by contigs.
    IDSource<long long int> source(index.getTotalLength());
    
    // This is the file we will dump our graph to, if needed.
    std::ofstream* dumpFile = NULL;
    if(options.count("dump")) {
        // Open it up.
        dumpFile = new std::ofstream("dump.dot");
        *dumpFile << "digraph dump {" << std::endl;
        
        // Start a cluster
        *dumpFile << "subgraph cluster_L0 {" << std::endl; 
        *dumpFile << "style=filled;" << std::endl;
        *dumpFile << "color=lightgrey;" << std::endl;
        *dumpFile << "label=\"Level 0\";" << std::endl;
        
        // Add per-contig rank constraints.
        for(size_t i = 0; i < index.getNumberOfContigs(); i++) {
            // For every contig, add rank edges.
            for(size_t j = 0; j < index.getContigLength(i) - 1; j++) {
                // For every base except the last...
                // Add an edge from this one to the next one, to enforce order.
                *dumpFile << "N" << i << "B" << j + 1 << 
                "->N" << i << "B" << j + 2 << "[style=invis];" << std::endl;
            }
        }
        
    }
    
    // Make a thread set for the context length we want.
    stPinchThreadSet* threadSet = mergeNonsymmetric(index, contextLength,
        dumpFile, options.count("quiet"));
        
    if(options.count("dump")) {
        // End the cluster and start a new one.
        *dumpFile << "}" << std::endl << "subgraph cluster_L1 {" << std::endl;
        *dumpFile << "style=filled;" << std::endl;
        *dumpFile << "color=lightgrey;" << std::endl;
        *dumpFile << "label=\"Level 1\";" << std::endl;
    }
    
    // Index it so we have a bit vector and SmallSides to write out.
    std::pair<RangeVector*, std::vector<SmallSide> > levelIndex;
    
    if(options.count("scan")) {
        // Use a scanning strategy for indexing.
        levelIndex = makeLevelIndexScanning(threadSet, index, contextLength,
            source, dumpFile);
     } else {
        // Use a range-based strategy for indexing.
        levelIndex = makeLevelIndex(threadSet, index, contextLength, source,
            dumpFile);
    }
        
    // Write it out, deleting the bit vector in the process
    saveLevelIndex(levelIndex, indexDirectory + "/level1");
    
    // Clean up the thread set
    stPinchThreadSet_destruct(threadSet);
    
    if(options.count("dump")) {
        // Finish and clean up the dump file
        *dumpFile << "}" << std::endl << "}" << std::endl;
        dumpFile->close();
        delete dumpFile;
    }
    
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

    // Now we're done!
    return 0;
}
