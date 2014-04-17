#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>

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
 * Look at the text of a base to see what arrow it needs to
 * have. Can also pass a base "strand" number.
 */
std::string getArrow(CSA::usint textNumber) {
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
std::string buildIndex(std::string indexDirectory,
    std::vector<std::string> fastas) {

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
stPinchThreadSet* makeThreadSet(const FMDIndex& index) {
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
 *
 * Note that due to the nature of this merging scheme, any two nodes that would
 * merge at a longer context length will also merge at a shorter context length,
 * so we can just directly calculate each upper level in turn.
 */
stPinchThreadSet* mergeNonsymmetric(const FMDIndex& index,
    size_t contextLength, std::ostream* dumpFile = NULL) {
    
    // Keep track of nodes we have already dumped.
    std::map<std::string, bool> dumped;
    
    // Make a thread set from our index.
    stPinchThreadSet* threadSet = makeThreadSet(index);
    
    // To construct the non-symmetric merged graph with p context:
    // Traverse the suffix tree down to depth p + 1
    for(CSA::FMD::iterator i = index.fmd.begin(contextLength); 
        i != index.fmd.end(contextLength); ++i) {
        // For each pair of suffix and position in the suffix tree
        
        // Unpack the iterator into pattern and FMDPosition at which it happens.
        std::string pattern = (*i).first;
        CSA::FMDPosition range = (*i).second;
        
        // Dump the context and range.
        std::cout << pattern << " at " << range << std::endl;
        
        // Check by counting again
        CSA::pair_type count = index.fmd.count(pattern);
        std::cout << "Count: (" << count.first << "," << count.second << ")" <<
            std::endl;
        
        if(range.end_offset >= 1 || dumpFile != NULL) {
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
            CSA::usint* locations = index.fmd.locate(oneStrandRange);
                
            if(locations == NULL) {
                throw std::runtime_error("Coud not locate pair!");
            }
            
            // For each base location, we need to work out the contig and base
            // and orientation, and pinch with the first.
            
            // Work out what text and base the first base is.
            CSA::pair_type firstBase = index.fmd.getRelativePosition(
                locations[0]);
            
            std::cout << "First relative position: text " << firstBase.first << 
                " offset " << firstBase.second << std::endl;
            
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
                
                    std::cout << "Making first: " << firstName << std::endl;
                    
                    // Write a node for it
                    *dumpFile << firstName << "[shape=\"record\",label=\"" << 
                        firstName << "\"];" << std::endl;
                    if(firstOffset > 1) {
                        // Link to previous Position
                        *dumpFile << index.getName(std::make_pair(
                            // Hack to get the base actually before us on the
                            // contig.
                            firstContigNumber * 2, firstOffset - 2)) << " -> " <<
                            firstName << ";" << std::endl;
                    }
                    dumped[firstName] = true;
                }
            }
            
            for(CSA::usint j = 1; j < range.end_offset + 1; j++) {
                // For each subsequent located base
                
                // Find and unpack it just like for the first base.
                CSA::pair_type otherBase = index.fmd.getRelativePosition(
                    locations[j]);
                //std::cout << "Relative position: (" << otherBase.first << "," << 
                //    otherBase.second << ")" << std::endl;   
                CSA::usint otherContigNumber = index.getContigNumber(otherBase);
                CSA::usint otherStrand = index.getStrand(otherBase);
                CSA::usint otherOffset = index.getOffset(otherBase);
                
                // Grab the other pinch thread.
                stPinchThread* otherThread = stPinchThreadSet_getThread(
                    threadSet, otherContigNumber);
                    
                if(firstThread == NULL) {
                    throw std::runtime_error("Other thread was NULL!");
                }
            
                // Pinch firstBase on firstNumber and otherBase on otherNumber
                // in the correct relative orientation.
                //std::cout << "\tPinching #" << firstContigNumber << ":" << 
                //    firstOffset << " strand " << firstStrand << " and #" << 
                //    otherContigNumber << ":" << otherOffset << " strand " << 
                //    otherStrand << std::endl;
                
                stPinchThread_pinch(firstThread, otherThread, firstOffset,
                    otherOffset, 1, firstStrand != otherStrand);
                    
                if(dumpFile != NULL) {
                    // Report the other position as existing.
                    std::string otherName = index.getName(otherBase); 
                    
                    if(!dumped.count(otherName)) {
                    
                        std::cout << "Making other: " << otherName << std::endl;
                        
                        // Write a node for it
                        *dumpFile << otherName << 
                            "[shape=\"record\",label=\"" << otherName << 
                            "\"];" << std::endl;
                        
                        if(otherOffset > 1) {
                            // Link previous position to us.
                            *dumpFile << index.getName(std::make_pair(
                                // Hack to get the base actually before us on the
                                // contig.
                                otherContigNumber * 2, otherOffset - 2)) << 
                                " -> " << otherName << ";" << std::endl;
                        }
                        dumped[otherName] = true;
                    }
                }
                
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
    
    // Return the finished thread set
    return threadSet;
}

/**
 * Make the range vector and list of matching Sides for the hierarchy level
 * implied by the given thread set in the given index, assuming it was created
 * by going through all contexts of the given length and assigning them to
 * nodes. Gets IDs for created positions from the given source.
 * 
 * If dumpFile is set, also writes a graphviz-format debug graph.
 * 
 * Don't forget to delete the bit vector when done!
 */
std::pair<CSA::RLEVector*, std::vector<Side> > makeLevelIndex(
    stPinchThreadSet* threadSet, const FMDIndex& index, size_t contextLength,
    IDSource<long long int>& source, std::ofstream* dumpFile = NULL) {
    
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
    
    // Now go through all the contexts again.
    for(CSA::FMD::iterator i = index.fmd.begin(contextLength); 
        i != index.fmd.end(contextLength); ++i) {
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
        CSA::pair_type base = index.fmd.getRelativePosition(index.fmd.locate(
            oneStrandRange.first));
        
        std::cout << "First relative position: text " << base.first << 
                " offset " << base.second << std::endl;
        
        // What contig corresponds to that text?
        CSA::usint contigNumber = index.getContigNumber(base);
        // And what strand corresponds to that text?
        CSA::usint strand = index.getStrand(base);
        // And what base position is that from the front of the contig?
        CSA::usint offset = index.getOffset(base);
     
        // Now we need to look up what the pinch set says is the canonical
        // position for this base, and what orientation it should be in. All the
        // bases in this range have the same context and should thus all be
        // pointing to the canonical base's replacement.
        
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
            positionID = idReservations[contigAndOffset] = source.next();
            
            if(dumpFile != NULL) {
                // Put in this new position
                *dumpFile << positionID << "[shape=\"record\",label=\"" << 
                    positionID << "\"];" << std::endl;
            }
        }
        
        if(dumpFile != NULL) {
            // Locate every base in oneStrandRange
            CSA::usint* locations = index.fmd.locate(oneStrandRange);
            
            std::cout << "Located bases (" << oneStrandRange.first << "," <<
                oneStrandRange.second << ")" << std::endl;
            
            for(CSA::usint j = 0; j < range.end_offset + 1; j++) {
                // For each location
                
                // Get its text and offset.
                CSA::pair_type childBase = index.fmd.getRelativePosition(
                    locations[j]);
                
                std::cout << "\tAt text " << childBase.first << " offset " << 
                    childBase.second << std::endl;
                    
                // Link from the child base to the merged base it goes into.
                *dumpFile << index.getName(childBase) << " -> " << positionID <<
                    "[color=red,dir=both,arrowtail=" << 
                    getArrow(childBase.first) << ",arrowhead=" << 
                    getArrow(canonicalOrientation) << "];" << std::endl;
                
            }
            
            // Throw out the location data.
            // TODO: Did they use new?
            free(locations);
            
        }
        
        
        // Record a 1 in the vector at the end of this range, in BWT
        // coordinates.
        encoder.addBit(range.reverse_start + range.end_offset + 
            index.fmd.getNumberOfSequences());
          
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
    
    // Finish the vector encoder into a vector of the right length.
    CSA::RLEVector* bitVector = new CSA::RLEVector(encoder,
        index.fmd.getBWTRange().second);
    
    // Return the bit vector and the Side vector
    return make_pair(bitVector, mappings);
    
}

/**
 * Save both parts of the given level index to files in the given directory,
 * which must not yet exist. Also deletes the bit vector, so don't use that
 * level index again.
 */
void saveLevelIndex(std::pair<CSA::RLEVector*, std::vector<Side> > levelIndex,
    std::string directory) {
    
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
int main(int argc, char** argv) {
    // Parse options with boost::programOptions. See
    // <http://www.radmangames.com/programming/how-to-use-boost-program_options>

    std::string appDescription = 
        "Create a reference hierarchy for mapping to FASTAs";

    // Make an options description for our program's options.
    boost::program_options::options_description description("Options");
    // Add all the options
    description.add_options() 
        ("help", "Print help messages") 
        ("dump", "Dump GraphViz graphs")
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
            
    } catch(boost::program_options::error& error) {
        // Something is bad about our options. Complain on stderr
        std::cerr << "Option parsing error: " << error.what() << std::endl;
        std::cerr << std::endl; 
        std::cerr << description << std::endl; 
        
        // Stop the program.
        return -1; 
    }
    
    // TODO: define context for merging in a more reasonable way (argument?)
    int contextLength = 3;
    
    // If we get here, we have the right arguments. Parse them.
    
    // This holds the directory for the reference structure to build.
    std::string indexDirectory(options["indexDirectory"].as<std::string>());
    
    // This holds a list of FASTA filenames to load and index.
    std::vector<std::string> fastas(options["fastas"]
        .as<std::vector<std::string> >());
    
    // Dump options.
    std::cout << "Options:" << std::endl;
    std::cout << "Store index in: " << indexDirectory << std::endl;
    
    for(std::vector<std::string>::iterator i = fastas.begin();
        i != fastas.end(); ++i) {
        
        std::cout << "Index file: " << *i << std::endl;
    }
    
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
    }
    
    // Make a thread set for the context length we want.
    stPinchThreadSet* threadSet = mergeNonsymmetric(index, contextLength,
        dumpFile);
        
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
