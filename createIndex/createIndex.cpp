#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>

#include <boost/filesystem.hpp>

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

/**
 * Get the contig number (not the text number) from a (text, offset) pair.
 */
CSA::usint getContigNumber(CSA::pair_type base) {
    // What contig corresponds to that text? Contigs all have both strands.
    return base.first / 2;
}

/**
 * Get the strand from a (text, offset) pair: either 0 or 1.
 */
CSA::usint getStrand(CSA::pair_type base) {
    // What strand corresponds to that text? Strands are forward, then reverse.
    return base.first % 2;
}

/**
 * Given a pair_type representing a base, and a vector of contig lengths by
 * number, determine the specified base's offset from the left of its contig,
 * 1-based.
 */
CSA::usint getOffset(CSA::pair_type base,
    const std::vector<long long int>& contigLengths) {
    // What base offset, 1-based, from the left, corresponds to this pair_type,
    // which may be on either strand.
    if(getStrand(base) == 0) {
        // We're on the forward strand; just correct for the 1-basedness offset.
        return base.second + 1;
    } else {
        // We're on the reverse strand, so we measured from the end.
        return contigLengths[getContigNumber(base)] - base.second;
    }
}

/**
 * createIndex: command-line tool to create a multi-level reference structure.
 */
int main(int argc, char** argv) {
    if(argc < 3) {
        // They forgot their arguments.
        std::cout << "Usage: " << argv[0] << " <index directory> <fasta> "
            << "[<fasta> [<fasta> ...]]" << std::endl;
        std::cout << "The index directory will be deleted and recreated, if it "
            << "exists." << std::endl;
        return 1;
    }
    
    // TODO: define context for merging in a more reasonable way (argument?)
    int contextLength = 3;
    
    // If we get here, we have the right arguments.
    
    // This holds the directory for the reference structure to build.
    std::string indexDirectory(argv[1]);
    
    // This holds a list of FASTA filenames to load and index.
    std::vector<std::string> fastas;
    for(int i = 2; i < argc; i++) {
        // Put each filename in the vector.
        fastas.push_back(std::string(argv[i]));
    }
    
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
    
    
    // Load the index with RLCSA.
    CSA::FMD index(basename);
    
    // Open the contig name/length file for reading.
    std::ifstream contigFile((basename + ".chrom.sizes").c_str());
    
    // Read all the contigs. We'll keep them in this map of contig name to
    // length.
    std::map<std::string, long long int> contigLengths;
    // Also we need a vector of contig names to give them all numbers.
    std::vector<std::string> contigNames;
    // Also a vector of lengths by contig number.
    std::vector<long long int> contigLengthVector;
    
    // Also a variable for the next available ID, which comes after 2 * all the
    // contigs.
    long long int nextID = 0;
    
    // Also a string to hold each line in turn.
    std::string line;
    while(std::getline(contigFile, line)) {
        // For each <contig>\t<length> pair...
        
        // Find the \t
        size_t separator = line.find('\t');
        
        // Split into name and length
        std::string contigName = line.substr(0, separator - 1);
        std::string contigLength = line.substr(separator + 1, line.size() - 1);
        
        // Parse the length
        long long int lengthNumber = atoll(contigLength.c_str());
        
        // Add it to the map
        contigLengths[contigName] = lengthNumber;
        
        // And to the vector of names in number order.
        contigNames.push_back(contigName);
        
        // And the vector of sizes in number order
        contigLengthVector.push_back(lengthNumber);
        
        // Also keep track of the IDs this contig is held to have used up.
        nextID += 2 * lengthNumber;
    }
    // Close up the contig file. We read our map.
    contigFile.close();
    
    // Make a ThreadSet with one thread per contig.
    // Construct the thread set first.
    stPinchThreadSet* threadSet = stPinchThreadSet_construct();
    
    for(size_t i = 0; i < contigNames.size(); i++) {
        // Add a thread for this contig number, starting at 1 and having the
        // appropriate length.
        stPinchThreadSet_addThread(threadSet, i, 1, 
            contigLengths[contigNames[i]]);
    } 
    
    // To construct the non-symmetric merged graph with p context:
    // Traverse the suffix tree down to depth p + 1
    
    for(CSA::FMD::iterator i = index.begin(contextLength); 
        i != index.end(contextLength); ++i) {
        // For each pair of suffix and position in the suffix tree
        
        // Unpack the iterator into pattern and FMDPosition at which it happens.
        std::string pattern = (*i).first;
        CSA::FMDPosition range = (*i).second;
        
        // Dump the context and range.
        std::cout << pattern << " at " << range << std::endl;
        
        // Check by counting again
        CSA::pair_type count = index.count(pattern);
        std::cout << "Count: (" << count.first << "," << count.second << ")" <<
            std::endl;
        
        if(range.end_offset >= 1) {
            // We only need to do any pinching if this context appears in two or
            // more places. And since a range with offset 0 has one thing in it,
            // we check to see if it's 1 or more.
            
            // Grab just the one-strand range, and construct the actual
            // endpoint. Use the reverse strand so that we handle things in
            // increasing coordinate order.
            CSA::pair_type oneStrandRange = std::make_pair(range.reverse_start, 
                range.reverse_start + range.end_offset);
            
            std::cout << "Locating (" << oneStrandRange.first << "," <<
                oneStrandRange.second << ")" << std::endl;
            
            // Locate the entire range. We'll have to free this eventually.
            // There are at least two of these.
            CSA::usint* locations = index.locate(oneStrandRange);
                
            if(locations == NULL) {
                throw std::runtime_error("Coud not locate pair!");
            }
            
            // For each base location, we need to work out the contig and base
            // and orientation, and pinch with the first.
            
            // Work out what text and base the first base is.
            CSA::pair_type firstBase = index.getRelativePosition(locations[0]);
            
            std::cout << "Relative position: (" << firstBase.first << "," << 
                firstBase.second << ")" << std::endl;
            
            // What contig corresponds to that text?
            CSA::usint firstContigNumber = getContigNumber(firstBase);
            // And what strand corresponds to that text?
            CSA::usint firstStrand = getStrand(firstBase);
            // And what base position is that?
            CSA::usint firstOffset = getOffset(firstBase, contigLengthVector);
            
            // Grab the first pinch thread
            stPinchThread* firstThread = stPinchThreadSet_getThread(threadSet,
                firstContigNumber);
            
            for(CSA::usint j = 1; j < range.end_offset + 1; j++) {
                // For each subsequent located base
                
                // Find and unpack it just like for the first base.
                CSA::pair_type otherBase = index.getRelativePosition(
                    locations[j]);
                //std::cout << "Relative position: (" << otherBase.first << "," << 
                //    otherBase.second << ")" << std::endl;   
                CSA::usint otherContigNumber = getContigNumber(otherBase);
                CSA::usint otherStrand = getStrand(otherBase);
                CSA::usint otherOffset = getOffset(otherBase,
                    contigLengthVector);
                
                // Grab the other pinch thread.
                stPinchThread* otherThread = stPinchThreadSet_getThread(
                    threadSet, otherContigNumber);
            
                // Pinch firstBase on firstNumber and otherBase on otherNumber
                // in the correct relative orientation.
                //std::cout << "\tPinching #" << firstContigNumber << ":" << 
                //    firstOffset << " strand " << firstStrand << " and #" << 
                //    otherContigNumber << ":" << otherOffset << " strand " << 
                //    otherStrand << std::endl;
                
                stPinchThread_pinch(firstThread, otherThread, firstOffset,
                    otherOffset, 1, firstStrand != otherStrand);
                
            }
            
            // Free the location buffer
            free(locations);
            
            // Say we merged some bases.
            std::cout << "Merged " << range.end_offset + 1 <<  " bases" <<
                std::endl;
            
        }
        
        // Now GC the boundaries in the pinch set
        std::cout << "Joining trivial boundaries..." << std::endl;
        stPinchThreadSet_joinTrivialBoundaries(threadSet);
        
    }
    
    // What are we going to save?
    
    // A bit vector denoting ranges, which we encode with this encoder, which
    // has 32 byte blocks.
    CSA::RLEEncoder encoder(32);
    
    // A vector of Sides, which are (contig, base, face).
    std::vector<Side> mappings;
    
    // We also need this map of position IDs (long long ints) by canonical
    // contig name (a size_t) and base index (a CSA::usint)
    std::map<std::pair<size_t, CSA::usint>, long long int>
        idReservations;
    
    // Now go through all the contexts again.
    for(CSA::FMD::iterator i = index.begin(contextLength); 
        i != index.end(contextLength); ++i) {
        // For each pair of suffix and position in the suffix tree, in
        // increasing SA coordinate order.
        
        // Unpack the iterator into pattern and FMDPosition at which it happens.
        std::string pattern = (*i).first;
        CSA::FMDPosition range = (*i).second;
        
        // Dump the context and range.
        std::cout << "Reprocessing " << pattern << " at " << range << std::endl;
        
        if(range.end_offset == -1) {
            // Not even a single thing with this context. Skip to the next one.
            continue;
        }
        
        // Grab just the one-strand range, and construct the actual
        // endpoint. Use the reverse strand so that we handle things in
        // increasing coordinate order.
        CSA::pair_type oneStrandRange = std::make_pair(range.reverse_start, 
            range.reverse_start + range.end_offset);
            
        // Work out what text and base the first base is.
        CSA::pair_type base = index.getRelativePosition(index.locate(
            oneStrandRange.first));
        
        std::cout << "Relative position: (" << base.first << "," << 
            base.second << ")" << std::endl;
        
        // What contig corresponds to that text?
        CSA::usint contigNumber = getContigNumber(base);
        // And what strand corresponds to that text?
        CSA::usint strand = getStrand(base);
        // And what base position is that from the front of the contig?
        CSA::usint offset = getOffset(base, contigLengthVector);
     
        // Now we need to look up what the pinch set says is the canonical
        // position for this base, and what orientation it should be in.
        
        // Get the segment
        stPinchSegment* segment = stPinchThreadSet_getSegment(threadSet, 
            contigNumber, offset);
            
        if(segment == NULL) {
            throw std::runtime_error("Found position in null segment!");
        }
            
        // How is it oriented?
        bool segmentOrientation = stPinchSegment_getBlockOrientation(segment);
        // How far into the segment are we?
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
        
        // What's the offset into the canonical segment?
        CSA::usint canonicalSegmentOffset = segmentOffset;
        if(segmentOrientation != canonicalOrientation) {
            // We really want this many bases in from the end of the contig, not
            // out from the start.
            // TODO: Is this 0-based or 1-based?
            canonicalSegmentOffset = stPinchSegment_getLength(
                firstSegment) - canonicalSegmentOffset;
        }
        // What is the offset in the canonical sequence? TODO: needs to be
        // 1-based.
        CSA::usint canonicalOffset = canonicalSegmentOffset + 
            stPinchSegment_getStart(firstSegment);
        
        // So what's the contig-and-offset pair?
        std::pair<size_t, CSA::usint> contigAndOffset = std::make_pair(
            canonicalContig, canonicalOffset); 
        
        // What Position ID does this base get?
        long long int positionID;
        
        if(idReservations.count(contigAndOffset) > 0) {
            // Load the previously chosen ID
            positionID = idReservations[contigAndOffset];
        } else {
            // Allocate and remember a new ID.
            positionID = idReservations[contigAndOffset] = nextID++;
        }
        
        // Record a 1 in the vector at the end of this range, in BWT
        // coordinates.
        encoder.addBit(range.reverse_start + range.end_offset + 
            index.getNumberOfSequences());
          
        // Say this range is going to belong to the ID we just looked up, on the
        // appropriate face.
        Side mapping;
        mapping.coordinate = positionID;
        mapping.face = canonicalOrientation ? LEFT : RIGHT;
            
        // Add a string describing the mapping
        mappings.push_back(mapping);
        
        // Every range belongs to some ID. TODO: test this. Especially with Ns.
        
        std::cout << "Canonicalized to #" << mapping.coordinate << "." << 
            mapping.face << std::endl;
    }
    
    
    // Throw out the threadset
    stPinchThreadSet_destruct(threadSet);
    
    // Finish the vector encoder into a vector of the right length.
    CSA::RLEVector bitVector(encoder, index.getBWTRange().second);
    
    // Save it to a file.
    std::ofstream vectorStream((indexDirectory + "/vector.bin").c_str());
    bitVector.writeTo(vectorStream);
    vectorStream.close();
    
    // Now write out all the merged positions to Avro, in an Avro-format file.
    // See <http://avro.apache.org/docs/1.7.6/api/cpp/html/index.html>. This is
    // complicated by the fact that we need to have the Avro schema JSON text to
    // write such a file, and the generated C++ classes provide no access to
    // that text.
    
    // We previously hacked the Side schema into Side_schema and Side_schema_len
    // with the xdd tool. Now we make it into an std::string.
    std::string sideSchema((char*) Side_schema, (size_t) Side_schema_len);
    std::istringstream sideSchemaStream(sideSchema);
    
    // This will hold the actual built schema
    avro::ValidSchema validSchema;
    
    // Now parse the schema. Assume it works.
    avro::compileJsonSchema(sideSchemaStream, validSchema);
    
    // Make a writer to write to the file we want, in the schema we want.
    avro::DataFileWriter<Side> writer((indexDirectory + 
        "/mappings.avro").c_str(), validSchema);
    
    for(std::vector<Side>::iterator i = mappings.begin(); i != mappings.end();
        ++i) {
        
        // Write each mapping to the file.
        writer.write(*i);
        
    }
    writer.close();
    
    // Now we're done!
    
    
    return 0;
}
