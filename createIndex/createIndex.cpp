#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <unordered_set>
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
#include <NaturalMappingScheme.hpp>
#include <ZipMappingScheme.hpp>

// Grab timers from libsuffixtools
#include <Timer.h>


#include "IDSource.hpp"
#include "ConcurrentQueue.hpp"
#include "MappingMergeScheme.hpp"
#include "MergeApplier.hpp"

#include "unixUtil.hpp"
#include "adjacencyComponentUtil.hpp"
#include "pinchGraphUtil.hpp"

#include "indexUtil.hpp"

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
    std::unordered_set<stPinchBlock*> seen;
    
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
 * a vector of canonicalized TextPositions.
 *
 * All BWT positions must be represented in the pinch set.
 *
 * Scans through the entire BWT.
 */
std::pair<GenericBitVector*, std::vector<TextPosition>>
identifyMergedRuns(
    stPinchThreadSet* threadSet, 
    const FMDIndex& index,
    const GenericBitVector* mask = NULL
) {
    
    // We need to make bit vector denoting ranges, which we encode with this
    // encoder, which has 32 byte blocks.
    GenericBitVector* encoder = new GenericBitVector();
    
    // We also need to make a vector of canonical positions.
    std::vector<TextPosition> mappings;
    
    Log::info() << "Building merged run index by scan..." << std::endl;
    
    // Do the thing where we locate each base and, when the canonical position
    // changes, add a 1 to start a new range and add a mapping.

    // Keep track of the ID and relative orientation for the last position we
    // canonicalized.
    TextPosition lastCanonicalized;
    
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
        
        // Canonicalize it.
        TextPosition canonicalized = canonicalize(index, threadSet, base);
        
        if(j == 0 || canonicalized != lastCanonicalized) {
            // We need to start a new range here, because this BWT base maps to
            // a different position than the last one.
            
            // Say this range is going to belong to the canonical base.
            mappings.push_back(canonicalized);
            
            if(j >= index.getNumberOfContigs() * 2) {
                // Record a 1 in the vector at the start of every range,
                // including the first.
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
 * greedy merging scheme, in parallel. Returns the pinched thread set.
 *
 * The greedy merging scheme is to take one genome, map the next genome to it,
 * merge at the mappings, compute the upper-level index, map the next genome to
 * that, merge at the mappings, compute the upper level index, and so on.
 * 
 * Takes a factory function that can allocate new MappingSchemes for an index
 * and the ranges and mask bitvectors.
 *
 * If passed a StatTracker, will add in stats from every MappingScheme.
 */
stPinchThreadSet*
mergeGreedy(
    const FMDIndex& index,
    std::function<MappingScheme*(FMDIndexView&&)> mappingSchemeFactory,
    StatTracker* stats = nullptr
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
    // vector of canonicalized positions. Make sure only the selected positions
    // are included, because the graph mapping stuff doesn't work if you have
    // positions breaking up ranges that aren't masked in.
    auto mergedRuns = identifyMergedRuns(threadSet, index, includedPositions);
    
    for(size_t genome = 1; genome < index.getNumberOfGenomes(); genome++) {
        // For each genome that we have to merge in...
        
        // Make a map for all the merged positions, because that's what's
        // required by the view. TODO: Do this conversion earlier, and filter
        // out ranges that are just mapped to the expected places.
        std::map<size_t, TextPosition> rangesToPositions;
        for(size_t i = 0; i < mergedRuns.second.size(); i++) {
            
            // Enter every entry in the map. TODO: be selective, or just convert
            // to a full vector representation always.
            rangesToPositions[i] = mergedRuns.second[i];
        }
        
        // Make a new FMDIndexView, giving it the mask and ranges bitvectors.
        FMDIndexView view(index, includedPositions, mergedRuns.first,
            rangesToPositions);
        
        // Allocate a new MappingScheme using the view. 
        MappingScheme* mappingScheme = mappingSchemeFactory(std::move(view));
            
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
        
        if(stats != nullptr) {
            Log::info() << "Copying over stats after merge" << std::endl;
            // Save stats if applicable
            *(stats) += mappingScheme->getStats();
        }
        
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
    
    // Register aborts (like uncaught std::bad_alloc) with the stack trace
    // handler too
    signal(SIGABRT, stacktraceOnSignal);
    
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
        ("alignment", boost::program_options::value<std::string>(), 
            "File to save .c2h-format alignment in")
        ("alignmentFasta", boost::program_options::value<std::string>(), 
            "File in which to save FASTA records for building HAL from .c2h")
        ("lastGraph", boost::program_options::value<std::string>(),
            "File in which to dump the merged graph in LastGraph format")
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
        ("mapStats", boost::program_options::value<std::string>(),
            "File in which to save the mapping stats from merging all levels")
        ("sampleRate", boost::program_options::value<unsigned int>()
            ->default_value(64), 
            "Set the suffix array sample rate to use")
        ("indexDirectory", boost::program_options::value<std::string>()
            ->required(), 
            "Directory to make the index in; will be deleted and replaced!")
        ("fastas", boost::program_options::value<std::vector<std::string> >()
            ->required()
            ->multitoken(),
            "FASTA files to load")
        ("scheme", boost::program_options::value<std::string>()
            ->default_value("greedy"),
            "Merging scheme (\"greedy\" only)")
        ("mapType", boost::program_options::value<std::string>()
            ->default_value("natural"),
            "Merging scheme (\"natural\" or \"zip\")")
        ("context", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum required context length to map on")
        ("credit", "Enable mapping on credit")
        ("mismatches", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Maximum allowed number of mismatches for credit")
        ("ignoreMatchesBelow", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Length below which to ignore maximal unique matches")
        ("minEditBound", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum *edit* distance lower bound on a maximum unique match")
        ("maxEditDistance", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Maximum *edit* distance from reference location")
        ("unstable", "Allow unstable mapping for increased coverage")
        ("maxRangeCount", boost::program_options::value<size_t>()
            ->default_value(100),
            "Maximum number of merged ranges to visit in a mapping step")
        ("maxExtendThrough", boost::program_options::value<size_t>()
            ->default_value(100),
            "Maximum number of bases to try to extend through");
        
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
    
    // Dump our hostname
    logHostname();
    
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
    
    // Grab the merging scheme we are going to use to merge.
    std::string mergeScheme = options["scheme"].as<std::string>();
    
    // Make a pointer to hold the threadset pointer it will create.
    stPinchThreadSet* threadSet;
    
    // We want to time the merge code.
    Timer* mergeTimer = new Timer("Merging");
    
    // We want to know the stats that the MappingScheme(s) make while we merge
    StatTracker stats;
    
    // We want to flag whether we want mismatches
    
    if(mergeScheme == "greedy") {
        // Use the greedy merge instead.
        threadSet = mergeGreedy(index, [&](FMDIndexView&& view) {
        
            // Make a new MappingScheme for this step and return a pointer to
            // it. TODO: would it be better to just make one MappingScheme and
            // let the ranges and mask be updated? Or passed to the map method?
            
            // Hiding our parameters by sneaking an option struct into a closure
            // seems a bit odd...
        
            if(options["mapType"].as<std::string>() == "natural") {
                // We want a NaturalMappingScheme
                NaturalMappingScheme* scheme = new NaturalMappingScheme(
                    std::move(view));
                    
                // Populate it
                scheme->credit = options.count("credit");
                scheme->minContext = options["context"].as<size_t>();
                scheme->z_max = options["mismatches"].as<size_t>();
                scheme->ignoreMatchesBelow = options[
                    "ignoreMatchesBelow"].as<size_t>();
                scheme->minHammingBound = options[
                    "minEditBound"].as<size_t>();
                scheme->maxHammingDistance = options[
                    "maxEditDistance"].as<size_t>();
                scheme->unstable = options.count("unstable");
                
                return (MappingScheme*) scheme;
            } else if(options["mapType"].as<std::string>() == "zip") {
                // Make a ZipMappingScheme, which can handle graphs. But we need
                // the forward and reverse versions of merged ranges to agree on
                // what positions they are assigned when merging, but the view
                // takes care of that.
                ZipMappingScheme<FMDPosition>* scheme =
                    new ZipMappingScheme<FMDPosition>(std::move(view));
                
                // Set the parameters from the arguments
                scheme->minContextLength = options["context"].as<size_t>();
                scheme->maxRangeCount = options["maxRangeCount"].as<size_t>();
                scheme->maxExtendThrough =
                    options["maxExtendThrough"].as<size_t>();
                scheme->minUniqueStrings = options["minEditBound"].as<size_t>();
                
                // Set up credit
                scheme->credit.enabled = options.count("credit");
                scheme->credit.maxMismatches =
                    options["mismatches"].as<size_t>();
                    
                return (MappingScheme*) scheme;
            } else {
                // They asked for a mapping scheme we don't have.
                throw std::runtime_error("Invalid mapping scheme: " +
                    options["mapType"].as<std::string>());
            }
        }, &stats);
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
        
        // We are using the greedy scheme, but there may be alignment between
        // things not in the first genome, so we have to use this fake-root
        // output format.
        
        size_t rootSeqLength = writeAlignment(threadSet, index,
            options["alignment"].as<std::string>());
            
        if(options.count("alignmentFasta")) {
            // Also save a FASTA, with the root seq
            writeAlignmentFasta(fastas, rootSeqLength,
                options["alignmentFasta"].as<std::string>());
        }
        
    }
    
    if(options.count("lastGraph")) {
        // Save in as good an approximation of LastGraph format as we can get at
        // the moment.
        writeLastGraph(threadSet, options["lastGraph"].as<std::string>());
    }
    
    // Make an IDSource to produce IDs not already claimed by contigs.
    IDSource<long long int> source(index.getTotalLength());
    
    // This will hold the computed level index of the merged level.
    std::pair<GenericBitVector*, std::vector<SmallSide> > levelIndex;
    
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

    if(options.count("mapStats")) {
        // We should save the mapping scheme stats (combined) to a file.
        Log::output() << "Saving stats to " <<
            options["mapStats"].as<std::string>() << std::endl;
            
        stats.save(options["mapStats"].as<std::string>());
    }

    Log::output() << "Final memory usage:" << std::endl;
    logMemory();

    // Now we're done!
    return 0;
}
