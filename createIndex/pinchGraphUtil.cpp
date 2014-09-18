#include "pinchGraphUtil.hpp"

#include <Log.hpp>
#include <unordered_set>


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

std::pair<std::pair<size_t, size_t>, bool>
canonicalize(
    stPinchThreadSet* threadSet, 
    size_t contigNumber,
    size_t offset,
    bool strand
) {
    
    Log::trace() << "Canonicalizing " << contigNumber << ":" << offset << 
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
    
    Log::trace() << "Canonicalized segment offset " << segmentOffset << 
        " to " << canonicalSegmentOffset << std::endl;
    
    // What's the offset into the canonical contig? 1-based because we add a
    // 0-based offset to a 1-based position.
    size_t canonicalOffset = stPinchSegment_getStart(firstSegment) + 
        canonicalSegmentOffset;
    
    Log::trace() << "Canonicalized contig " << contigNumber << " offset " <<
        offset << " to contig " << canonicalContig << " offset " << 
        canonicalOffset << std::endl;
    
    // Return all three values, and be sad about not having real tuples. What
    // orientation should we use?  Well, we have the canonical position's
    // orientation within the block, our position's orientation within the
    // block, and the orientation that this context attaches to the position in.
    // Flipping any of those will flip the orientation in which we need to map,
    // so we need to xor them all together, which for bools is done with !=.
    Log::trace() << "Canonical orientation: " << canonicalOrientation << 
            std::endl;
    Log::trace() << "Segment orientation: " << segmentOrientation << 
            std::endl;
    Log::trace() << "Strand: " << strand << std::endl;

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

std::pair<std::pair<size_t, size_t>, bool>
canonicalize(
    const FMDIndex& index, 
    stPinchThreadSet* threadSet, 
    TextPosition base
) {
    
    Log::trace() << "Canonicalizing 0-based " << base << std::endl;
    
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

size_t
writeAlignment(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    const std::string& filename
) {

    // We're going to lay out the segments in our root sequence with some
    // locality. We're going to scan through all threads, and put each new block
    // as it is encountered. So we need a set of the encountered blocks.
    std::unordered_set<stPinchBlock*> seen;
    
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
                
            stPinchBlock* block = stPinchSegment_getBlock(segment);
            if(block != NULL) {
                // It actually aligned
            
                // Are we in the same orientation as the root?
                bool orientation = 
                    (stPinchSegment_getBlockOrientation(segment) == 
                    stPinchSegment_getBlockOrientation(
                    stPinchBlock_getFirst(block)));
            
                // Write a top segment mapping to the segment named after the
                // address of the block this segment belongs to.
                c2h << "a\t" << segmentStart << "\t" << 
                    stPinchSegment_getLength(segment) << "\t" << 
                    (uintptr_t)block << "\t" << orientation << std::endl;
                    
                Log::debug() << "Bottom segment orientation: " << 
                    stPinchSegment_getBlockOrientation(segment) << std::endl;
                
                Log::debug() << "First segment orientation: " << 
                    stPinchSegment_getBlockOrientation(
                    stPinchBlock_getFirst(block)) << std::endl;
                        
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

void
writeAlignmentWithReference(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index, 
    const std::string& filename, 
    const size_t referenceGenomeNumber
) {
    
    Log::info() << "Serializing star-tree alignment to " << filename <<
        std::endl;
    
    // Open up the file to write.
    std::ofstream c2h(filename.c_str());
    
    // Keep a mapping from genome number to event name. Event name will be
    // either the contig scaffold name if all the contigs in the genome are from
    // a single scaffold, or "genome-<number>" if there are multiple scaffolds
    // involved.
    std::map<size_t, std::string> eventNames;
    
    for(size_t genome = 0; genome < index.getNumberOfGenomes(); genome++) {
        // Go through genomes in order and assign them all event names.
        
        // What contigs are in this genome?
        auto genomeRange = index.getGenomeContigs(genome);
        
        // This is the event name we are going to use.
        std::string eventName;
        
        for(size_t i = genomeRange.first; i < genomeRange.second; i++) {
            if(i == genomeRange.first) {
                // Grab the name of the first contig
                eventName = index.getContigName(i);
            } else if(index.getContigName(i) != eventName) {
                // This contig doesn't match the first contig's name; they are
                // not all from the same scaffold. Re-name the event with a new
                // generic name.
                
                eventName = "genome-" + 
                    std::to_string(index.getContigGenome(i));
                    
                // Now we're done
                break;
            }
        }
        
        // We have now named the event for this genome. Save it
        eventNames[genome] = eventName;
        
        Log::info() << "Named genome " << genome << "/" << 
            index.getNumberOfGenomes() << " event " << eventName << std::endl;
    }

    // What contigs are in the reference genome?
    auto referenceRange = index.getGenomeContigs(referenceGenomeNumber);
    
    // What blocks have we seen the reference genome be part of?
    std::unordered_set<stPinchBlock*> seen;
    
    // Keep track of where the last segment ended, across contigs on a scaffold.
    size_t lastSegmentEnd;
    
    // And the sequence name for the last contig we did
    std::string lastContigName;
    
    for(size_t i = referenceRange.first; i < referenceRange.second; i++) {
        Log::info() << "Processing reference contig " << i << std::endl;
        
        // Go through all threads in the reference
        stPinchThread* thread = stPinchThreadSet_getThread(threadSet,
            i);
        
        if(i == referenceRange.first || 
            index.getContigName(i) != lastContigName) {
            
            // Start a new sequence for the reference with a sequence line. The
            // sequence is not a top sequence.
            c2h << "s\t'" << eventNames[referenceGenomeNumber] << "'\t'" << 
                index.getContigName(i) << "'\t1" << std::endl;
            
            // Start over at 0 on a new scaffold
            lastContigName = index.getContigName(i);
            lastSegmentEnd = 0;
        }
            
        // Get the start offset fom 0 on the scaffold for the first block in the
        // thread.
        size_t contigStart = index.getContigStart(i);
        
        // Get the first segment in the thread.
        stPinchSegment* segment = stPinchThread_getFirst(thread);
        while(segment != NULL) {
            // Go through all the segments in the thread
            
            // Work out where this segment starts in the reference. Convert to
            // 0-based HAL.
            size_t segmentStart = contigStart +
                stPinchSegment_getStart(segment) - 1;
            
            if(segmentStart > lastSegmentEnd) {
                // We need an unaligned segment padding out to here. Make it
                // named after 1 less than this segment's address, where no
                // other segment could possibly fit.

                c2h << "a\t" << ((size_t)(uintptr_t) segment) - 1 << "\t" << 
                    lastSegmentEnd << "\t" << segmentStart - lastSegmentEnd <<
                    std::endl;
                
            }
            
            // This holds the block that this segment is in, if any.
            stPinchBlock* block;
            if((block = stPinchSegment_getBlock(segment)) != NULL) {
                // Look at its block if it has one
                
                if(seen.count(block) == 0) {
                    // This is a new block we haven't seen before on this
                    // thread.
                    
                    // Put a bottom segment for the block. Just name the segment
                    // after the block's address. We need block name (number),
                    // start, and length.
                    c2h << "a\t" << (uintptr_t)block << "\t" << 
                        segmentStart << "\t" << 
                        stPinchBlock_getLength(block) << std::endl;
                        
                    // Record that we have seen this block now.
                    seen.insert(block);
                    
                    // Make sure the reference segments are first in their
                    // blocks.
                    stPinchSegment_putSegmentFirstInBlock(segment);
                    
                } else {
                    // The reference is going back though a block we already saw
                    // reading along the reference. We can't peoperly serialize
                    // as a star tree due to self-alignment in the reference.
                    
                    throw std::runtime_error("Can't make a star-tree due to "
                        "self-alignment in reference!");
                }
            } else {
                // Make a segment named after this segment's address instead.
                // Nothing will map to it.
                
                c2h << "a\t" << (uintptr_t) segment << "\t" << 
                    segmentStart << "\t" << 
                    stPinchSegment_getLength(segment) << std::endl;
            }
        
            // Keep track of where this segment ended
            lastSegmentEnd = segmentStart + stPinchSegment_getLength(segment);
        
            // Jump to the next 3' segment. This needs will return NULL if we go
            // off the end.
            segment = stPinchSegment_get3Prime(segment);
        }
    }
    
    // Now we've made all the bottom segments for the reference. For all the
    // other genomes it's just like what we do when we have no reference.
    
    for(size_t genome = 0; genome < index.getNumberOfGenomes(); genome++) {
        if(genome == referenceGenomeNumber) {
            // Don't do the reference since it is already the root.
            continue;
        }
        
        Log::info() << "Processing query genome " << genome << std::endl;
        
        // What contigs are in this genome? We assume they are grouped by source
        // scaffold and sorted in order of increasing scaffold position within a
        // scaffold.
        auto genomeRange = index.getGenomeContigs(genome);
        
        // What scaffold was the last contig in this genome on? We can only have
        // one sequence line and coordinate space per scaffold, remember.
        std::string lastContigScaffold;
        
        // Keep track of where the last segment ended, across contigs on a
        // scaffold.
        size_t lastSegmentEnd;
        
        for(size_t i = genomeRange.first; i < genomeRange.second; i++) {
        
            Log::info() << "Processing query contig " << i << std::endl;
        
            if(i == genomeRange.first || 
                index.getContigName(i) != lastContigScaffold) {
                
                // We are starting on the contigs for a new scaffold.
                
                // Start a new sequence with a sequence line. The sequence is a
                // top sequence, since it is only connected up.
                c2h << "s\t'" << eventNames[genome] << "'\t'" << 
                    index.getContigName(i) << "'\t0" << std::endl;
                    
                // Remember not to make a sequence line for this scaffold again.
                lastContigScaffold = index.getContigName(i);
                
                // Start over from 0 on the scaffold.
                lastSegmentEnd = 0;
            }
        
            // Each segment has to account for the offset of the contig on the
            // sequence.
            size_t contigStart = index.getContigStart(i);
            
            // So first get the thread by name.
            stPinchThread* thread = stPinchThreadSet_getThread(threadSet, i);
            
            // Go through all its pinch segments in order. There's no iterator
            // so we have to keep looking 3'
            
            // Get the first segment in the thread.
            stPinchSegment* segment = stPinchThread_getFirst(thread);
            while(segment != NULL) {
                
                // Start where the next segment starts. Convert from 1-based
                // pinch segments to 0-based HAL.
                size_t segmentStart = contigStart +
                    stPinchSegment_getStart(segment) - 1;
                    
                if(segmentStart > lastSegmentEnd) {
                    // We need an unaligned segment padding out to here.
                    c2h << "a\t" << lastSegmentEnd << "\t" << 
                        segmentStart - lastSegmentEnd << std::endl;
                }
                
                stPinchBlock* block = stPinchSegment_getBlock(segment);
                
                if(block != NULL) {
                    // It actually aligned
                
                    if(seen.count(block) == 0) {
                        // This aligned block wasn't in the reference. Complain
                        // we can't represent it.
                        throw std::runtime_error("Star tree can't represent "
                            "alignment of sequence not aligned to reference");
                    }
                
                    // Are we in the same orientation as the root?
                    bool orientation = 
                        (stPinchSegment_getBlockOrientation(segment) == 
                        stPinchSegment_getBlockOrientation(
                        stPinchBlock_getFirst(block)));
                        
                    Log::debug() << "Bottom segment (" << segment << 
                        ") orientation: " << 
                        stPinchSegment_getBlockOrientation(segment) << 
                        std::endl;
                        
                    Log::debug() << "First segment (" << 
                        stPinchBlock_getFirst(block) << ") orientation: " <<
                        stPinchSegment_getBlockOrientation(
                        stPinchBlock_getFirst(block)) << std::endl;
                
                    // Write a top segment mapping to the segment named after
                    // the address of the block this segment belongs to. Fields
                    // are a, start, length, aligned bottom segment,
                    // orientation.
                    c2h << "a\t" << segmentStart << "\t" << 
                        stPinchSegment_getLength(segment) << "\t" << 
                        (uintptr_t) block << "\t" << orientation << std::endl;
                        
                } else {
                    // Write a segment for the unaligned sequence.
                    c2h << "a\t" << segmentStart << "\t" << 
                        stPinchSegment_getLength(segment) << std::endl;
                }
                
                // Keep track of where this segment ended
                lastSegmentEnd = segmentStart + stPinchSegment_getLength(
                    segment);
                
                // Jump to the next 3' segment. This will return NULL if we
                // go off the end.
                segment = stPinchSegment_get3Prime(segment);
            }
                
        }
    
    }
    
    // Now that we did every segment of every contig of every genome, we are
    // done and can close up the file.
    c2h.close();
}

void 
writeBottomSegments(
    std::ofstream& c2h, 
    stPinchThread* thread
) {
    // Get the first segment in the thread.
    stPinchSegment* segment = stPinchThread_getFirst(thread);
    while(segment != NULL) {
        // Go through all the segments in the thread
        
        // Work out where this segment starts in the reference. Convert to
        // 0-based HAL.
        size_t segmentStart = stPinchSegment_getStart(segment) - 1;
        
        // No part of the thread is not covered by a segment.
        
        // This holds the block that this segment is in, if any.
        stPinchBlock* block = stPinchSegment_getBlock(segment);
        
        // This holds the name of the segment we're adding. Name it after the
        // block if we have one, and the unaligned segment otherwise.
        uintptr_t segmentName = (block != NULL) ? (uintptr_t) block : 
            (uintptr_t) segment;
            
        // Put a bottom segment in the c2h
        c2h << "a\t" << segmentName << "\t" << segmentStart << "\t" << 
            stPinchSegment_getLength(segment) << std::endl;
            
        if(block != NULL) {
            // Make sure these segments are first, so we can easily check the
            // orientations of the top segments against them.
            // TODO: Assumes these are the root.
            stPinchSegment_putSegmentFirstInBlock(segment);
        }
    
        // Jump to the next 3' segment. This needs will return NULL if we go
        // off the end.
        segment = stPinchSegment_get3Prime(segment);
    }
}

void 
writeTopSegments(
    std::ofstream& c2h, 
    stPinchThread* thread
) {
    // Get the first segment in the thread.
    stPinchSegment* segment = stPinchThread_getFirst(thread);
    while(segment != NULL) {
        // Go through all the segments in the thread
        
        // Work out where this segment starts. Convert to 0-based HAL.
        size_t segmentStart = stPinchSegment_getStart(segment) - 1;
        
        // No part of the thread is not covered by a segment.
        
        // This holds the block that this segment is in, if any.
        stPinchBlock* block = stPinchSegment_getBlock(segment);
        
        // Start with the common part between aligned and unaligned blocks.
        c2h << "a\t" << segmentStart << "\t" << 
            stPinchSegment_getLength(segment);
        
        if(block != NULL) {
            // Are we in the same orientation as the root?
            // TODO: This assumes root is the first segment.
            bool orientation = 
                (stPinchSegment_getBlockOrientation(segment) == 
                stPinchSegment_getBlockOrientation(
                stPinchBlock_getFirst(block)));
        
            // Write the bit for aligning
            c2h << "\t" << (uintptr_t) block << "\t" << orientation;
                
            Log::debug() << "Bottom segment orientation: " << 
                stPinchSegment_getBlockOrientation(segment) << std::endl;
                
            Log::debug() << "First segment orientation: " << 
                stPinchSegment_getBlockOrientation(
                stPinchBlock_getFirst(block)) << std::endl;
        }
        
        // Finish the line
        c2h << std::endl;
        
        // Jump to the next 3' segment. This needs will return NULL if we go
        // off the end.
        segment = stPinchSegment_get3Prime(segment);
    }
}

void
writeAlignmentWithReference(
    stPinchThreadSet* threadSet, 
    std::vector<std::string> threadNames,
    std::vector<std::string> threadEvents,
    const std::string& filename, 
    const size_t referenceThreadNumber
) {

    // The same as above but without all that code to manage contigs spanning
    // multiple threads or to calculate names for things.
    
    Log::info() << "Serializing star-tree alignment to " << filename <<
        std::endl;
    
    // Open up the file to write.
    std::ofstream c2h(filename.c_str());
    
    Log::info() << "Processing reference thread " << referenceThreadNumber <<
        std::endl;
        
    // Grab the reference thread.
    stPinchThread* reference = stPinchThreadSet_getThread(threadSet,
        referenceThreadNumber);
        
    if(reference == NULL) {
        // We may have gotten empty inputs or something.
        throw std::runtime_error(std::string("Reference ") + 
            std::to_string(referenceThreadNumber) + " does not exist");
    }
    
    // Start a new sequence for the reference with a sequence line. The
    // sequence is not a top sequence.
    c2h << "s\t'" << threadEvents[referenceThreadNumber] << "'\t'" << 
        threadNames[referenceThreadNumber] << "'\t1" << std::endl;
        
    // Write all the bottom segments for the reference.
    writeBottomSegments(c2h, reference);
    
    // Get an iterator over the thread set, since thread numbers may have gaps.
    stPinchThreadSetIt iterator = stPinchThreadSet_getIt(threadSet);
    
    // This holds each thread we're going to look at.
    stPinchThread* thread;
    while((thread = stPinchThreadSetIt_getNext(&iterator)) 
        != NULL) {
        
        // Go through all the threads.
        // TODO: assumes we hit them in increasing numerical order.
    
        // Get the name number of the thread.
        size_t threadNumber = stPinchThread_getName(thread);
        
        if(threadNumber == referenceThreadNumber) {
            // Already did this one; it's the reference
            continue;
        }
        
        Log::info() << "Processing query thread " << threadNumber <<
            std::endl;
        
        // Write the top sequence header
        c2h << "s\t'" << threadEvents[threadNumber] << "'\t'" << 
            threadNames[threadNumber] << "'\t0" << std::endl;        
        
        // Write the top segments for the thread
        writeTopSegments(c2h, stPinchThreadSet_getThread(threadSet,
            threadNumber));
        
    }
    
    // And we're done.
    c2h.close();
}

void
writeAlignmentFasta(
    std::vector<std::string> inputFastas,
    int64_t rootBases,
    std::string filename
) {

    // Open the FASTA to write
    std::ofstream fasta(filename.c_str());
    
    Log::info() << "Generating " << rootBases << 
        " bases of root node sequence." << std::endl;
    
    
    if(rootBases >= 0) {
        // First we put the right number of Ns in a sequence named "rootSeq", if
        // the caller didn't ask to not have it.
        fasta << ">rootSeq" << std::endl;
        for(size_t i = 0; i < rootBases; i++) {
            fasta << "N";
            // Entire sequence must be on one line.
        }
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
