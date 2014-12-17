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


#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// Grab pinchesAndCacti dependency.
#include <stPinchGraphs.h>

// Grab all the libFMD stuff.
#include <FMDIndex.hpp>
#include <FMDIndexBuilder.hpp>
#include <GenericBitVector.hpp>
#include <FMDIndexIterator.hpp>
#include <TextPosition.hpp>
#include <util.hpp>
#include <Mapping.hpp>
#include <SmallSide.hpp>
#include <Log.hpp>
#include <MappingScheme.hpp>
#include <LRMappingScheme.hpp>
#include <NaturalMappingScheme.hpp>

// Grab timers from libsuffixtools
#include <Timer.h>


#include "IDSource.hpp"
#include "ConcurrentQueue.hpp"
#include "OverlapMergeScheme.hpp"
#include "MappingMergeScheme.hpp"
#include "MergeApplier.hpp"

#include "unixUtil.hpp"
#include "adjacencyComponentUtil.hpp"
#include "pinchGraphUtil.hpp"

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
 * Canonicalize each contigous run of positions mapping to the same canonical
 * base and face.
 *
 * Takes a pinched thread set, and the index on which it is defined.
 *
 * If a mask bit vector is specified, only positions which have a 1 in the mask
 * will be able to break up ranges.
 *
 * Returns a GenericBitVector marking each such range with a 1 at the start, and
 * a vector of canonicalized positions. Each anonicalized position is a contig
 * number, a base number, and a face flag.
 *
 * All BWT positions must be represented in the pinch set.
 *
 * Scans through the entire BWT.
 */
std::pair<GenericBitVector*, 
    std::vector<std::pair<std::pair<size_t, size_t>, bool>>> 
identifyMergedRuns(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index,
    const GenericBitVector* mask = NULL
) {
    
    // We need to make bit vector denoting ranges, which we encode with this
    // encoder, which has 32 byte blocks.
    GenericBitVector* encoder = new GenericBitVector();
    
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
        
        if(mask != NULL && !mask->isSet(j)) {
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
                encoder->addBit(j);
                Log::trace() << "Set bit " << j << std::endl;
            }
            
            
            // Remember what canonical base and face we're doing for this range.
            lastCanonicalized = canonicalized;
        }
        // Otherwise we had the same canonical base, so we want this in the same
        // range we already started.
    }
            
    // Set a bit after the end of the last range (i.e. at the end of the BWT).
    encoder->addBit(index.getBWTLength());
    
    // Finish the vector at the right length, leaving room for that trailing
    // bit.
    encoder->finish(index.getBWTLength() + 1);
    
    // Return the bit vector and the canonicalized base vector
    return std::make_pair(encoder, mappings);
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
std::pair<GenericBitVector*, std::vector<SmallSide> > 
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
    std::pair<GenericBitVector*, std::vector<SmallSide> > levelIndex,
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
 * Create a new thread set from the given FMDIndex, and merge it down by the
 * greedy merging scheme, in parallel. Returns the pinched thread set.
 *
 * The greedy merging scheme is to take one genome, map the next genome to it,
 * merge at the mappings, compute the upper-level index, map the next genome to
 * that, merge at the mappings, compute the upper level index, and so on.
 * 
 * Takes a factory function that can allocate new MappingSchemes for an index
 * and the ranges and mask bitvectors.
 */
stPinchThreadSet*
mergeGreedy(
    const FMDIndex& index,
    std::function<MappingScheme*(const FMDIndex&,
        const GenericBitVector& ranges, const GenericBitVector* mask)> 
        mappingSchemeFactory
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
    const GenericBitVector* includedPositions = &index.getGenomeMask(0);
    
    // Canonicalize everything, yielding a bitvector (pointer) of ranges and a
    // vector of canonicalized positions.
    auto mergedRuns = identifyMergedRuns(threadSet, index);
    
    for(size_t genome = 1; genome < index.getNumberOfGenomes(); genome++) {
        // For each genome that we have to merge in...
        
        // Allocate a new MappingScheme, giving it the ranges and mask bitvectors.
        MappingScheme* mappingScheme = mappingSchemeFactory(index,
            *mergedRuns.first, includedPositions);
        
        // Make the merge scheme we want to use. We choose a mapping-to-second-
        // level-based merge scheme, to which we need to feed the details of the
        // second level (range vector, representative positions to merge into)
        // and also the bitmask of what bottom-level things to count. We also
        // need to tell it what genome to map the contigs of.
        MappingMergeScheme scheme(index, mappingScheme, genome);

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
        
        // Delete the mapping scheme, so we can delete the stuff it uses
        delete mappingScheme;
        
        // Merge the new genome into includedPositions, replacing the old
        // bitvector.
        GenericBitVector* newIncludedPositions = includedPositions->createUnion(
            index.getGenomeMask(genome));
        
        if(genome > 1) {
            // If we already alocated a new GenericBitVector that wasn't the one
            // that came when we loaded in the genomes, we need to delete it.
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
        // And, if we had to make any additional included position
        // GenericBitVectors, get the last one of those too.
        delete includedPositions;
    }
    
    // Let the caller rebuild that from this pinch graph which we
    // return.
    return threadSet;
}

/**
 * Save an adjacency component spectrum to a file as a <size>\t<count> TSV.
 */
void
writeAdjacencyComponentSpectrum(
    std::map<size_t, size_t> spectrum,
    std::string filename
) {

    Log::info() << "Saving adjacency component spectrum to " <<
        filename << std::endl;

    // Open up the file to write.
    std::ofstream file(filename.c_str());
    
    for(auto kv : spectrum) {
        file << kv.first << "\t" << kv.second << std::endl;
    }    
    
    file.close();
    
}

/**
 * Save a vector of numbers as a single-column TSV.
 */
template<typename T>
void
writeColumn(
    std::vector<T> numbers,
    std::string filename
) {

    // Open up the file to write.
    std::ofstream file(filename.c_str());
    
    for(auto number : numbers) {
        // Write each number on its own line
        file << number << std::endl;
    }    
    
    file.close();
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
    signal(SIGINT, stacktraceOnSignal);
    
    // Register segfaults with the stack trace handler
    signal(SIGSEGV, stacktraceOnSignal);
    
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
        ("spectrum", boost::program_options::value<std::string>(), 
            "File in which to save graph adjacency component size spectrum")
        ("indelLengths", boost::program_options::value<std::string>(), 
            "File in which to save indel lengths between a pair of genomes")
        ("tandemDuplications", boost::program_options::value<std::string>(), 
            "File in which to save the number of tandem duplications")
        ("nontrivialRearrangements", 
            boost::program_options::value<std::string>(), 
            "File in which to dump nontrivial rearrangements")
        ("context", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum required context length to merge on")
        ("addContext", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Extra context beyond that needed to be unique for greedy LR")
        ("multContext", boost::program_options::value<double>()
            ->default_value(0), 
            "Minimum context length as a fraction of uniqueness distance")
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
            "FASTA files to load")
        ("credit", "Enable mapping on credit")
        ("mapType", boost::program_options::value<std::string>()
            ->default_value("LR"),
            "Merging scheme (\"natural\" or \"LR\")")
        ("mismatches", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Maximum allowed number of mismatches")
        ("ignoreMatchesBelow", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Length below which to ignore maximal unique matches")
        ("minHammingBound", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum Hamming distance lower bound on a maximum unique match")
        ("maxHammingDistance", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Maximum Hamming distance from reference location");
        
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
            // These are required.
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
    
    // Grab the context types we are going to use to merge
    std::string mapType = options["mapType"].as<std::string>();
    
    if(mapType != "LR" && mapType != "natural") {
        // They gave us a map type we don't have probably.
        Log::critical() << "Invalid mapType: " << mapType << std::endl;
        return 1;
    }
    
    // Grab the merging scheme we are going to use to merge.
    std::string mergeScheme = options["scheme"].as<std::string>();
    
    // Make a pointer to hold the threadset pointer it will create.
    stPinchThreadSet* threadSet;
    
    // We want to time the merge code.
    Timer* mergeTimer = new Timer("Merging");
    
    // We want to flag whether we want mismatches
    
    if(mergeScheme == "overlap") {
        // Make a thread set that's all merged, with the given minimum merge
        // context.
        threadSet = mergeOverlap(index, options["context"].as<size_t>());
    } else if(mergeScheme == "greedy") {
        // Use the greedy merge instead.
        threadSet = mergeGreedy(index, [&](const FMDIndex& index,
            const GenericBitVector& ranges, const GenericBitVector* mask) {
        
            // Make a new MappingScheme for this step and return a pointer to
            // it. TODO: would it be better to just make one MappingScheme and
            // let the ranges and mask be updated? Or passed to the map method?
            // Hiding our parameters by sneaking an option struct into a closure
            // seems a bit odd...
        
            if(options["mapType"].as<std::string>() == "LR") {
                // We want an LRMappingScheme
                LRMappingScheme* scheme = new LRMappingScheme(index, ranges,
                    mask);
                    
                // Populate it
                scheme->minContext = options["context"].as<size_t>();
                scheme->addContext = options["addContext"].as<size_t>();
                scheme->multContext = options["multContext"].as<double>();
                scheme->credit = options.count("credit");
                scheme->z_max = options["mismatches"].as<size_t>();
                
                return (MappingScheme*) scheme;
            } else if(options["mapType"].as<std::string>() == "natural") {
                // We want a NaturalMappingScheme
                NaturalMappingScheme* scheme = new NaturalMappingScheme(index,
                    ranges, mask);
                    
                // Populate it
                scheme->credit = options.count("credit");
                scheme->minContext = options["context"].as<size_t>();
                scheme->z_max = options["mismatches"].as<size_t>();
                scheme->ignoreMatchesBelow = options[
                    "ignoreMatchesBelow"].as<size_t>();
                scheme->minHammingBound = options[
                    "minHammingBound"].as<size_t>();
                scheme->maxHammingDistance = options[
                    "maxHammingDistance"].as<size_t>();
                
                return (MappingScheme*) scheme;
            } else {
                // They asked for a mapping scheme we don't have.
                throw std::runtime_error("Invalid mapping scheme: " +
                    options["mapType"].as<std::string>());
            }
        });
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
        if(mergeScheme == "overlap") {
            // We can have self-alignment in the first sequence, so we need to
            // use this serializer.
        
            // Save the alignment defined by the pinched pinch graph to the file
            // the user specified. Save the number of bases of root sequence
            // that were used in the center of the star tree.
            size_t rootBases = writeAlignment(threadSet, index,
                options["alignment"].as<std::string>());
                
            if(options.count("alignmentFasta")) {
                // Also save a FASTA with the sequences necessary to generate a
                // HAL from the above.
                writeAlignmentFasta(fastas, rootBases,
                    options["alignmentFasta"].as<std::string>());
            }
        } else {
            // We are using the greedy scheme and can't have self-alignment in
            // the first genome, so we can use this star tree serializer.
            
            writeAlignmentWithReference(threadSet, index, 
                options["alignment"].as<std::string>(), 0);
                
            if(options.count("alignmentFasta")) {
                // Also save a FASTA, without the rootSeq sequence at all.
                writeAlignmentFasta(fastas, -1,
                    options["alignmentFasta"].as<std::string>());
            }
            
        }
    }
    
    // Make an IDSource to produce IDs not already claimed by contigs.
    IDSource<long long int> source(index.getTotalLength());
    
    // This will hold the computed level index of the merged level.
    std::pair<GenericBitVector*, std::vector<SmallSide> > levelIndex;
    
    // We also want to time the merged level index building code
    Timer* levelIndexTimer = new Timer("Level Index Construction");
    
    // Use a scanning strategy for indexing.
    levelIndex = makeLevelIndexScanning(threadSet, index, source);
    
    delete levelIndexTimer;
        
    // Write it out, deleting the bit vector in the process
    saveLevelIndex(levelIndex, indexDirectory + "/level1");
    
    // Now, while we still have the threadSet, we can work out the adjacency
    // components. 
    auto components = getAdjacencyComponents(threadSet);
    
    if(options.count("spectrum")) {
        // How many adjacency components are what size? Adjacency components of
        // size 2 are just SNPs or indels, while adjacency components of larger
        // sizes are generally more complex rearrangements.
        auto spectrum = getAdjacencyComponentSpectrum(components);
    
        // Save a dump of pinch graph adjacency component sizes
        writeAdjacencyComponentSpectrum(spectrum,
            options["spectrum"].as<std::string>());
    }
    
    if(options.count("indelLengths")) {
        // Get all the size-2 components, determine indel lengths for them, and
        // save them to the file the user wanted them in.
        writeColumn(getIndelLengths(filterComponentsBySize(components, 2)),
            options["indelLengths"].as<std::string>()); 
    }
    
    if(options.count("tandemDuplications")) {
        // Get all the size-4 components, and count tandem duplications.
        size_t tandemDuplications = countTandemDuplications(
            filterComponentsBySize(components, 4));
          
        // Hack the count into a 1-element vector and write it to the file.  
        std::vector<size_t> tandemDupeVector {tandemDuplications};
        writeColumn(tandemDupeVector,
            options["tandemDuplications"].as<std::string>()); 
    }
    
    if(options.count("nontrivialRearrangements")) {
        // Get all size>2 things and dump them.
        writeAdjacencyComponents(filterComponentsBySize(components, 2, 
            [](size_t a, size_t b) { return a > b; }), 
            options["nontrivialRearrangements"].as<std::string>());
    }
    
    // Clean up the thread set after we analyze everything about it.
    stPinchThreadSet_destruct(threadSet);
    
    // Get rid of the range vector
    delete levelIndex.first;
    
    // Get rid of the index itself. Invalidates the index reference.
    delete indexPointer;

    Log::output() << "Final memory usage:" << std::endl;
    logMemory();

    // Now we're done!
    return 0;
}
