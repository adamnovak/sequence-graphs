#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <set>
#include <algorithm>
#include <utility>

#include <boost/filesystem.hpp>
#include "boost/program_options.hpp" 

#include <rlcsa/fmd.h>
#include <rlcsa/bits/rlevector.h>

// Grab pinchesAndCacti dependency.
#include "stPinchGraphs.h"

// Grab the Avro header for the Face/Side objects we need to dump out. Only use
// the most dependent one, since they all define their dependencies and are thus
// incompatible with the headers for their dependencies.
#include "schemas/Side.hpp"

// Get the variables Side_schema and Side_schema_len that give us the actual
// text of that schema, before code generation.
#include "schemas/Side_schema.hpp"

#include <avro/Encoder.hh>
#include <avro/DataFile.hh>
#include <avro/Schema.hh>
#include <avro/ValidSchema.hh>
#include <avro/Compiler.hh>

#include "FMDIndexBuilder.hpp"
#include "FMDIndex.hpp"
#include "IDSource.hpp"

/**
 * Define a macro for easily compiling in/out detailed debugging information.
 * Replaces the one we got from the fmd header.
 */
//#define DEBUG(op) op
#define DEBUG(op)


/**
 * Look at the text of a base to see what arrow it needs to
 * have. Can also pass a base "strand" number.
 */
std::string 
getArrow(
    CSA::usint textNumber
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
 * given FASTAs for the bottom level FMD index. Returns the basename of the FMD
 * index that gets created.
 */
std::string
buildIndex(
    std::string indexDirectory,
    std::vector<std::string> fastas
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
    FMDIndexBuilder builder(basename);
    for(std::vector<std::string>::iterator i = fastas.begin(); i < fastas.end();
        ++i) {
        
        std::cout << "Adding FASTA " << *i << std::endl;
        
        // Add each FASTA file to the index.
        builder.add(*i);
    }
    
    // Return the basename
    return basename;
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
    
    for(size_t i = 0; i < index.lengths.size(); i++) {
        // Add a thread for this contig number. The thread starts at base 1 and
        // has the appropriate length.
        stPinchThreadSet_addThread(threadSet, i, 1, index.lengths[i]);
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
    for(CSA::FMD::iterator i = index.fmd.begin(contextLength); 
        i != index.fmd.end(contextLength); ++i) {
        // For each pair of suffix and position in the suffix tree
        
        // Unpack the iterator into pattern and FMDPosition at which it happens.
        std::string pattern = (*i).first;
        // This is in SA coordinates.
        CSA::FMDPosition range = (*i).second;
        
        
        if(!quiet) {
            // Dump the context and range.
            std::cout << pattern << " at " << range << std::endl;
        }
        
        if(range.end_offset >= 1 || dumpFile != NULL) {
            // We only need to do any pinching if this context appears in two or
            // more places. And since a range with offset 0 has one thing in it,
            // we check to see if it's 1 or more.
            
            // Grab just the one-strand range, and construct the actual
            // endpoint.
            CSA::pair_type oneStrandRange = std::make_pair(range.forward_start, 
                range.forward_start + range.end_offset);
            
            DEBUG(std::cout << "Locating (" << oneStrandRange.first << "," <<
                oneStrandRange.second << ")" << std::endl;)
            
            // Locate the entire range. We'll have to free this eventually.
            // There are at least two of these.
            CSA::usint* locations = index.fmd.locate(oneStrandRange);
                
            if(locations == NULL) {
                throw std::runtime_error("Coud not locate pair!");
            }
            
            // For each base location, we need to work out the contig and base
            // and orientation, and pinch with the first.
            
            // Work out what text and base the first base is.
            CSA::pair_type firstBase = index.fmd.getRelativePosition(
                locations[0]);
            
            DEBUG(std::cout << "First relative position: text " << 
                firstBase.first << " offset " << firstBase.second << std::endl;)
            
            // What contig corresponds to that text?
            CSA::usint firstContigNumber = index.getContigNumber(firstBase);
            // And what strand corresponds to that text?
            CSA::usint firstStrand = index.getStrand(firstBase);
            // And what base position is that?
            CSA::usint firstOffset = index.getOffset(firstBase);
            
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
                    char baseChar = index.display(firstBase);
                    if(firstStrand) {
                        // Flip it around so we always see forward strand bases.
                        baseChar = CSA::reverse_complement(baseChar);
                    }
                
                    // Write a node for it
                    *dumpFile << firstName << "[shape=\"record\",label=\"{" << 
                        firstName << "|" << baseChar <<
                        "}\"];" << std::endl;
                    if(firstOffset > 1) {
                        // Link to previous Position
                        *dumpFile << index.getName(std::make_pair(
                            // Hack to get the base actually before us on the
                            // contig.
                            firstContigNumber * 2, firstOffset - 2)) << 
                            " -> " << firstName << ";" << std::endl;
                    }
                    dumped.insert(firstName);
                }
            }
            
            for(CSA::usint j = 1; j < range.end_offset + 1; j++) {
                // For each subsequent located base
                
                // Find and unpack it just like for the first base.
                CSA::pair_type otherBase = index.fmd.getRelativePosition(
                    locations[j]);
                DEBUG(
                    std::cout << "Relative position: (" << otherBase.first << 
                    "," << otherBase.second << ")" << std::endl;   
                )
                CSA::usint otherContigNumber = index.getContigNumber(otherBase);
                CSA::usint otherStrand = index.getStrand(otherBase);
                CSA::usint otherOffset = index.getOffset(otherBase);
                
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
                DEBUG(std::cout << "\tPinching #" << firstContigNumber << ":" <<
                    firstOffset << " strand " << firstStrand << " and #" << 
                    otherContigNumber << ":" << otherOffset << " strand " << 
                    otherStrand << " (orientation: " << orientation << ")" <<
                    std::endl;)
                
                stPinchThread_pinch(firstThread, otherThread, firstOffset,
                    otherOffset, 1, orientation);
                    
                if(dumpFile != NULL) {
                    // Report the other position as existing.
                    std::string otherName = index.getName(otherBase); 
                    
                    if(!dumped.count(otherName)) {
                        // Get its base
                        char baseChar = index.display(otherBase);
                        if(otherStrand) {
                            // Flip it around so we always see forward strand
                            // bases.
                            baseChar = CSA::reverse_complement(baseChar);
                        }
                    
                        // Write a node for it
                        *dumpFile << otherName << 
                            "[shape=\"record\",label=\"{" << otherName << "|" <<
                            baseChar << "}\"];" << std::endl;
                        
                        if(otherOffset > 1) {
                            // Link previous position to us.
                            *dumpFile << index.getName(std::make_pair(
                                // Hack to get the base actually before us on
                                // the contig.
                                otherContigNumber * 2, otherOffset - 2)) << 
                                " -> " << otherName << ";" << std::endl;
                        }
                        dumped.insert(otherName);
                    }
                }
                
            }
            
            // Free the location buffer
            free(locations);
            
            if(!quiet) {
                // Say we merged some bases.
                std::cout << "Merged " << range.end_offset + 1 <<  " bases" <<
                    std::endl;
            }
            
        }
    }
    
    // Now GC the boundaries in the pinch set
    std::cout << "Joining trivial boundaries..." << std::endl;
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    
    // Return the finished thread set
    return threadSet;
}

/**
 * Turn the given (contig number, 1-based offset from start, orientation)
 * position into the same sort of structure for the canonical base that
 * represents all the bases it has been pinched with.
 */
std::pair<std::pair<size_t, CSA::usint>, bool>
canonicalize(
    stPinchThreadSet* threadSet, 
    size_t contigNumber,
    CSA::usint offset,
    bool strand
) {
    
    DEBUG(std::cout << "Canonicalizing " << contigNumber << ":" << offset << 
        "." << strand << std::endl;)
    
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
    CSA::usint segmentOffset = offset - stPinchSegment_getStart(segment);
        
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
    CSA::usint canonicalSegmentOffset = segmentOffset;
    if(segmentOrientation != canonicalOrientation) {
        // We really want this many bases in from the end of the segment, not
        // out from the start of the segment. Keep it 0-based.
        canonicalSegmentOffset = stPinchSegment_getLength(firstSegment) - 
            canonicalSegmentOffset - 1;
    }
    
    DEBUG(std::cout << "Canonicalized segment offset " << segmentOffset << 
        " to " << canonicalSegmentOffset << std::endl;)
    
    // What's the offset into the canonical contig? 1-based because we add a
    // 0-based offset to a 1-based position.
    CSA::usint canonicalOffset = stPinchSegment_getStart(firstSegment) + 
        canonicalSegmentOffset;
    
    DEBUG(std::cout << "Canonicalized contig " << contigNumber << " offset " <<
        offset << " to contig " << canonicalContig << " offset " << 
        canonicalOffset << std::endl;)
    
    // Return all three values, and be sad about not having real tuples. What
    // orientation should we use?  Well, we have the canonical position's
    // orientation within the block, our position's orientation within the
    // block, and the orientation that this context attaches to the position in.
    // Flipping any of those will flip the orientation in which we need to map,
    // so we need to xor them all together, which for bools is done with !=.
    DEBUG(
        std::cout << "Canonical orientation: " << canonicalOrientation << 
            std::endl;
        std::cout << "Segment orientation: " << segmentOrientation << 
            std::endl;
        std::cout << "Strand: " << strand << std::endl;
    )
    return std::make_pair(std::make_pair(canonicalContig, canonicalOffset),
        canonicalOrientation != segmentOrientation != strand);
}

/**
 * Turn the given (text, offset) pair into a canonical (contig number, 1-based
 * offset from contig start, orientation), using the given FMDIndex and the
 * given thread set. The orientation is which face of the canonical base this
 * (text, offset) pair means.
 */
std::pair<std::pair<size_t, CSA::usint>, bool>
canonicalize(
    const FMDIndex& index, 
    stPinchThreadSet* threadSet, 
    CSA::pair_type base
) {
    
    // What contig corresponds to that text?
    CSA::usint contigNumber = index.getContigNumber(base);
    // And what strand corresponds to that text? This tells us what
    // orientation we're actually looking at the base in.
    bool strand = (bool) index.getStrand(base);
    // And what base position is that from the front of the contig? This is
    // 1-based.
    CSA::usint offset = index.getOffset(base);
    
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
    std::map<std::pair<size_t, CSA::usint>, long long int> idReservations, 
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
            
            DEBUG(std::cout << "Linking merged base " << j << " on thread " <<
                name << std::endl;)
            
            
            
            // Canonicalize its right side (strand 0 or false), which is the one
            // linked by this adjacency to the next base.
            std::pair<std::pair<size_t, CSA::usint>, bool> canonicalBase =
                canonicalize(threadSet, name, j, false);
            
            // Also canonicalize the left side of the next base.
            std::pair<std::pair<size_t, CSA::usint>, bool> canonicalNextBase =
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
 * Manages this by scanning the BWT from left to right, seeing what cannonical
 * position and orientation each base belongs to, and putting 1s in the range
 * vector every time a new one starts.
 *
 * If dumpFile is set, also writes a graphviz-format debug graph.
 * 
 * Don't forget to delete the bit vector when done!
 */
std::pair<CSA::RLEVector*, std::vector<Side> > 
makeLevelIndex(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    size_t contextLength,
    IDSource<long long int>& source, 
    std::ofstream* dumpFile = NULL
) {
    
    // We need to make bit vector denoting ranges, which we encode with this
    // encoder, which has 32 byte blocks.
    CSA::RLEEncoder encoder(32);
    
    // We also need to make a vector of Sides, which are (contig, base, face),
    // and are the things that get matched to by the corresponding ranges in the
    // bit vector.
    std::vector<Side> mappings;
    
    // We also need this map of position IDs (long long ints) by canonical
    // contig name (a size_t) and base index (a CSA::usint)
    std::map<std::pair<size_t, CSA::usint>, long long int>
        idReservations;
        
    // Work out where the valid positions in the BWT are.
    CSA::pair_type bwtBounds = index.fmd.getBWTRange();

    std::cout << "Building mapping data structure..." << std::endl;
    
    for(CSA::FMD::iterator i = index.fmd.begin(contextLength, true); 
        i != index.fmd.end(contextLength, true); ++i) {
        
        // Go through all the contexts, including shortened ones before end of
        // text.
        
        // Unpack the iterator into pattern and FMDPosition at which it happens.
        std::string context = (*i).first;
        // Keep in mind that this is in SA coordinates!
        CSA::FMDPosition range = (*i).second;
        
        // Work out where the forward range starts in BWT coordinates (again)
        CSA::usint bwtStart = range.forward_start;
        index.fmd.convertToBWTIndex(bwtStart);
        
        DEBUG(std::cout << "===Context: " << context << " at range " << range <<
            "===" << std::endl;)
        
        if(context.size() == contextLength) {
        
            // If it's a full-requested-length context (and therefore all the
            // bases in it have been pinched)...
            
            // Locate the first base in the context. It's already in SA
            // coordinates.
            CSA::pair_type base = index.fmd.getRelativePosition(
                index.fmd.locate(range.forward_start));
            DEBUG(std::cout << "Text/Offset: (" << base.first << ", " << 
                base.second << ")" << std::endl;)
            
            // Canonicalize it. The second field here will be the relative
            // orientation and determine the Side.
            std::pair<std::pair<size_t, CSA::usint>, bool> canonicalized = 
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
            
            DEBUG(
                char sourceChar = index.display(base);
                char destChar = index.display(canonicalized.first.first, 
                            canonicalized.first.second, canonicalized.second);
                std::cout << "Merged a " << sourceChar << " into a " << 
                    destChar << std::endl;
            )
            
            // Say this range is going to belong to the ID we just looked up, on
            // the appropriate face.
            Side mapping;
            mapping.coordinate = positionCoordinate;
            // What face?
            mapping.face = canonicalized.second ? LEFT : RIGHT;
            // Add the mapping
            mappings.push_back(mapping);
            
            if(bwtStart != bwtBounds.first) {
                // Record a 1 in the vector at the start of every range except
                // the first, in BWT coordinates. The first needs no 1 before it
                // so it will be rank 0 (and match up with mapping 0), and it's
                // OK not to split it off from the stop characters since they
                // can't ever be searched.
                encoder.addBit(bwtStart);
                DEBUG(std::cout << "Set bit " << bwtStart <<
                    std::endl;)
            }
            
            DEBUG(std::cout << "Mapping to #" << mapping.coordinate << 
                "." << mapping.face << std::endl;)
            
            // Put a 1 at the start of its interval, and add a mapping to the
            // given ID and Side.
            
            // TODO: Make edges by traversing the pinch graph.
                
        } else if(context.size() < contextLength) {
            // Otherwise if it's a shorter than requested context, it's
            // possible not everything has been merged. So do the old thing
            // where we locate each base and, when the canonical position
            // changes, add a 1 to start a new range and add a mapping.

            // Grab just the one-strand range, and construct the actual
            // endpoint.
            CSA::pair_type oneStrandRange = std::make_pair(range.forward_start, 
                range.forward_start + range.end_offset);
            
            // Locate the entire range. We'll have to free this eventually.
            // There are at least two of these.
            CSA::usint* locations = index.fmd.locate(oneStrandRange);
                
            if(locations == NULL) {
                throw std::runtime_error("Coud not locate short suffix!");
            }
            
            // Keep track of the ID and relative orientation for the last
            // position we canonicalized.
            std::pair<std::pair<size_t, CSA::usint>, bool> lastCanonicalized;
            
            for(CSA::usint j = 0; j < range.end_offset + 1; j++) {
                // For each base in the range...
                
                // Where was it located?
                CSA::pair_type base = index.fmd.getRelativePosition(
                    locations[j]);
                
                // Canonicalize it. The second field here will be the relative
                // orientation and determine the Side.
                std::pair<std::pair<size_t, CSA::usint>, bool> canonicalized = 
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
                    
                    DEBUG(
                        char sourceChar = index.display(base);
                        char destChar = index.display(canonicalized.first.first, 
                                    canonicalized.first.second,
                                    canonicalized.second);
                        std::cout << "Merged a " << sourceChar << " into a " << 
                            destChar << std::endl;
                    )
                    
                    // Say this range is going to belong to the ID we just
                    // looked up, on the appropriate face.
                    Side mapping;
                    mapping.coordinate = positionCoordinate;
                    // What face?
                    mapping.face = canonicalized.second ? LEFT : RIGHT;
                    // Add the mapping
                    mappings.push_back(mapping);
                    
                    // Add in the bit the same as above, only here we add the
                    // offset of i.
                    if(bwtStart + j != bwtBounds.first) {
                        // Record a 1 in the vector at the start of every range
                        // except the first, in BWT coordinates. The first needs
                        // no 1 before it so it will be rank 0 (and match up
                        // with mapping 0), and it's OK not to split it off from
                        // the stop characters since they can't ever be
                        // searched.
                        encoder.addBit(bwtStart + j);
                        DEBUG(std::cout << "Set bit " << bwtStart + j <<
                            std::endl;)
                    }
                    
                    DEBUG(std::cout << "Mapping to #" << 
                        mapping.coordinate << "." << mapping.face << std::endl;)
                    
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
    encoder.addBit(bwtBounds.second + 1);
    
    // Finish the vector encoder into a vector of the right length.
    // This should always end in a 1! Make sure to flush first.
    encoder.flush();
    CSA::RLEVector* bitVector = new CSA::RLEVector(encoder, 
        bwtBounds.second + 1);
    
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
 * which must not yet exist. Also deletes the bit vector, so don't use that
 * level index again.
 */
void saveLevelIndex(
    std::pair<CSA::RLEVector*, std::vector<Side> > levelIndex,
    std::string directory
) {
    
    std::cout << "Saving index to disk..." << std::endl;
    
    // Make the directory
    boost::filesystem::create_directory(directory);
    
    // Save the bit vector to a file.
    std::ofstream vectorStream((directory + "/vector.bin").c_str());
    levelIndex.first->writeTo(vectorStream);
    vectorStream.close();
    
    // Delete the bit vector since we're done with it.
    delete levelIndex.first;
    
    // Now write out all the merged positions to Avro, in an Avro-format file.
    // See <http://avro.apache.org/docs/1.7.6/api/cpp/html/index.html>. This is
    // complicated by the fact that we need to have the Avro schema JSON text to
    // write such a file, and the generated C++ classes provide no access to
    // that text.
    
    // We previously hacked the Side schema into globals Side_schema and
    // Side_schema_len with the xdd tool. Now we make it into an std::string.
    std::string sideSchema((char*) Side_schema, (size_t) Side_schema_len);
    std::istringstream sideSchemaStream(sideSchema);
    
    // This will hold the actual built schema
    avro::ValidSchema validSchema;
    
    // Now parse the schema. Assume it works.
    avro::compileJsonSchema(sideSchemaStream, validSchema);
    
    // Make a writer to write to the file we want, in the schema we want.
    avro::DataFileWriter<Side> writer((directory + 
        "/mappings.avro").c_str(), validSchema);
    
    for(std::vector<Side>::iterator i = levelIndex.second.begin(); 
        i != levelIndex.second.end(); ++i) {
        // Write each mapping to the file.
        writer.write(*i);
        
    }
    writer.close();
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
        ("quiet", "Don't print every context")
        ("context", boost::program_options::value<unsigned int>()
            ->default_value(3), 
            "Set the context length to merge on")
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
    std::cout << "Options:" << std::endl;
    std::cout << "Store index in: " << indexDirectory << std::endl;
    
    for(std::vector<std::string>::iterator i = fastas.begin();
        i != fastas.end(); ++i) {
        
        std::cout << "Index file: " << *i << std::endl;
    }
    
    std::cout << "Use " << contextLength << " bases of context." << std::endl;
    
    // Index the bottom-level FASTAs and get the basename they go into.
    std::string basename = buildIndex(indexDirectory, fastas);
    
    // Load the index and its metadata.
    FMDIndex index(basename);
    
    // Make an IDSource to produce IDs not already claimed by contigs.
    IDSource<long long int> source(index.getTotalLength());
    
    // This si the file we will dump our graph to, if needed.
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
        for(size_t i = 0; i < index.lengths.size(); i++) {
            // For every contig, add rank edges.
            for(size_t j = 0; j < index.lengths[i] - 1; j++) {
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
    
    // Index it so we have a bit vector and Sides to write out
    std::pair<CSA::RLEVector*, std::vector<Side> > levelIndex = makeLevelIndex(
        threadSet, index, contextLength, source, dumpFile);
        
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

    // Now we're done!
    return 0;
}
