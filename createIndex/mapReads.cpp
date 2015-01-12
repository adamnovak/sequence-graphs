// mapReads.cpp: program to map some reads (or other strings) to a string
// reference structure.
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
#include <Fasta.hpp>
#include <MappingScheme.hpp>
#include <LRMappingScheme.hpp>
#include <NaturalMappingScheme.hpp>
#include <OldNaturalMappingScheme.hpp>

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

#include "indexUtil.hpp"

/**
 * Load reads from the given FASTAs and queue them up in the given queue.
 * Returns total reads processed. Skips any reads with Ns.
 */
size_t
loadReads(
    const std::vector<std::string>& fastas,
    ConcurrentQueue<std::pair<std::string, std::string>>* recordsOut
) {

    size_t totalReads = 0;

    for(std::string readFasta : fastas) {
        // Now go through all the read FASTAs
        Fasta reader(readFasta);
        while(reader.hasNext()) {
            // Go through each FASTA record
            std::pair<std::string, std::string> headerAndSequence = 
                reader.getNextRecord();
            
            if(headerAndSequence.second.find('N') != std::string::npos) {
                // Skip this read
                Log::error() << "Skipping read with N: " << 
                    headerAndSequence.second << std::endl;
                continue;
            }
            
            totalReads++;
            Log::info() << "Mapping read " << totalReads << std::endl;
            
            // Put the read in the queue
            auto lock = recordsOut->lock();
            recordsOut->enqueue(headerAndSequence, lock);
        }
    }
    
    // Close the queue since we read everything.
    auto lock = recordsOut->lock();
    recordsOut->close(lock);
    
    return totalReads;
}

/**
 * Save each string in the queue to the given output stream. The strings already
 * have newlines.
 */
void
saveLines(
    ConcurrentQueue<std::string>* linesIn, std::ofstream& out
) {
    // Lock the record queue so we can maybe get a record    
    auto lineLock = linesIn->waitForNonemptyOrEnd();
    while(!linesIn->isEmpty(lineLock)) {
        
        // Dequeue and write every line
        out << linesIn->dequeue(lineLock);
    
        lineLock = linesIn->waitForNonemptyOrEnd();
    }
}

/**
 * Read FASTA sequence names and sequences from the input queue, map them to the
 * reference in the given index, according to the given mapping scheme, and send
 * lines of mapping TSV output to the output queue.
 *
 * Returns the total mappings made.
 */
size_t
mapSomeReads(
    ConcurrentQueue<std::pair<std::string, std::string>>* recordsIn, 
    const std::pair<std::string, std::string>& referenceRecord,
    const FMDIndex& index,
    const MappingScheme* mappingScheme,
    ConcurrentQueue<std::string>* linesOut
) {

    // We'll count all the mappings we make.
    size_t totalMappings = 0;

    // Lock the record queue so we can maybe get a record    
    auto recordLock = recordsIn->waitForNonemptyOrEnd();
    while(!recordsIn->isEmpty(recordLock)) {
        // We got a record to do. Dequeue it and unlock.
        std::pair<std::string, std::string> record = 
            recordsIn->dequeue(recordLock);
        
        // Parse out the record.
        std::string recordName = record.first;
        std::string sequence = record.second;
        
        // Map the sequence with the mapping scheme
        
        // Make a vector of mappings to populate.
        std::vector<Mapping> mappings(sequence.size());
        
        // Go and map
        mappingScheme->map(sequence, [&](size_t base, TextPosition mappedTo) {
            // For each TextPosition we get, make a Mapping
            mappings[base] = Mapping(mappedTo);
        });
        
        
        // Output each query base, noting which mapped and to where.
        for(size_t i = 0; i < sequence.size(); i++) {
            Mapping mapping = mappings[i];
            
            if(mapping.isMapped()) {
            
                // Did we map backwards or not?
                bool backwards = false;
            
                if(mapping.getLocation().getText() != 0) {
                    // Flip everything around to be on strand 0. Easy since
                    // there's exactly one contig indexed.
                    mapping = mapping.flip(index.getContigLength(0));
                    backwards = true;
                }
                
                // Get the character being mapped
                char mapped = sequence[i];
                // Get the character it is mapped to, from the forward strand.
                char mappedTo = referenceRecord.second[
                    mapping.getLocation().getOffset()];
                    
                if(backwards) {
                    // If we're mapping backwards, we have to match the reverse
                    // of that character.
                    mappedTo = complement(mappedTo);
                }
                    
                if(mapped != mappedTo) {
                    // Throw out this mapping, since it's placing a
                    // character on a character it doesn't match.
                    mapping = Mapping();
                    
                    Log::critical() << "Tried mapping " << recordName << ":" <<
                        i << " (" << mapped << ") to " << 
                        referenceRecord.first << ":" << 
                        mapping.getLocation().getOffset() << "." << backwards <<
                        " (" << mappedTo << ")" << std::endl;
                    
                    throw std::runtime_error("Non-matching mapping!");
                }
                
                
                // Now do the output for this line, because this position
                // mapped. Do it as reference, then query. And do it to this
                // stringstream.
                std::stringstream mappingStream;
                mappingStream << referenceRecord.first << "\t" << 
                    mapping.getLocation().getOffset() << "\t" << 
                    recordName << "\t" << 
                    i << "\t" << 
                    backwards << std::endl;
                
                // Make a string of the output
                std::string line(mappingStream.str());
                
                // Send it
                auto lineLock = linesOut->lock();
                linesOut->enqueue(line, lineLock);
                    
                // Count that we had a mapping
                totalMappings++;
                
            } else {
                // Report this query base as unaligned (just its contig and
                // base).
                
                std::stringstream mappingStream;
                mappingStream << recordName << "\t" << i << std::endl;
                
                // Make a string of the output
                std::string line(mappingStream.str());
                
                // Send it
                auto lineLock = linesOut->lock();
                linesOut->enqueue(line, lineLock);
                
            }
            
        }
        
        // Now wait for a new task, or for there to be no more records.
        recordLock = recordsIn->waitForNonemptyOrEnd();
        
    }
    
    // Close the output queue since we have run out of data.
    auto lineLock = linesOut->lock();
    linesOut->close(lineLock);
    
    // Give back the total mapping count.
    return totalMappings;

}

/**
 * mapReads: command-line tool to map strings to a string, using a reference
 * structure.
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
        std::string("Map strings to a string.\n") + 
        "Usage: mapReads <index directory> <reference> [<fasta> [<fasta> ...]] "
        "--alignment <alignment>";

    // Make an options description for our program's options.
    boost::program_options::options_description description("Options");
    // Add all the options
    description.add_options() 
        ("help", "Print help messages") 
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
        // These next three options should be ->required(), but that's not in
        // the Boost version I can convince our cluster admins to install. From
        // now on I shall work exclusively in Docker containers or something.
        ("indexDirectory", boost::program_options::value<std::string>(), 
            "Directory to make the index in; will be deleted and replaced!")
        ("reference", boost::program_options::value<std::string>(),
            "FASTA contining a single N-free reference sequence to map to")
        ("fastas", boost::program_options::value<std::vector<std::string> >()
            ->multitoken(),
            "FASTA files to map")
        ("alignment", boost::program_options::value<std::string>()
            ->required(), 
            "File to save alignment in, as a TSV of mappings")
        ("threads", boost::program_options::value<size_t>()
            ->default_value(16),
            "Number of mapping threads to run")
        ("credit", "Mapping on credit for greedy scheme")
        ("mismatches", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Maximum allowed number of mismatches")
        ("mapType", boost::program_options::value<std::string>()
            ->default_value("LR"),
            "Mapping scheme (\"natural\", \"old\", or \"LR\")")
        ("ignoreMatchesBelow", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Length below which to ignore maximal unique matches")
        ("minHammingBound", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum Hamming distance lower bound on a maximum unique match")
        ("maxHammingDistance", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Maximum Hamming distance from reference location")
        ("unstable", "Allow unstable mapping for increased coverage")
        ("stats", boost::program_options::value<std::string>(),
            "TSV file to save statistics to");
        
    // And set up our positional arguments
    boost::program_options::positional_options_description positionals;
    // One index directory
    positionals.add("indexDirectory", 1);
    // One reference
    positionals.add("reference", 1);
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
        
        if(!options.count("indexDirectory") || !options.count("reference") ||
            !options.count("fastas")) {
            
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
    
    // This holds the reference
    std::string reference(options["reference"].as<std::string>());
    
    // Read the reference contig. It must exist.
    std::pair<std::string, std::string> referenceRecord = 
        Fasta(reference).getNextRecord();
    
    // This holds a list of FASTA filenames to load and index.
    std::vector<std::string> fastas(options["fastas"]
        .as<std::vector<std::string> >());
        
    // Make a vector of just the reference.
    std::vector<std::string> referenceOnly = { reference };
        
    // Index the reference. Use the sample rate the user specified.
    FMDIndex* indexPointer = buildIndex(indexDirectory, referenceOnly,
        options["sampleRate"].as<unsigned int>());
        
    // Make a reference out of the index pointer because we're not letting it
    // out of our scope.
    FMDIndex& index = *indexPointer;
    
    // Make a vector that says every single position in the index is its own
    // range for mapping to.
    GenericBitVector ranges;
    for(size_t i = 0; i < index.getBWTLength(); i++) {
        ranges.addBit(i);
    }
    ranges.finish(index.getBWTLength());
    
    // Parse the number of threads to use
    size_t numThreads = options["threads"].as<size_t>();
    
    // Make a mapping scheme from the command-line options. TODO: unify with
    // createIndex's code for this.
    MappingScheme* mappingScheme;
    
    if(options["mapType"].as<std::string>() == "LR") {
        // We want an LRMappingScheme
        LRMappingScheme* scheme = new LRMappingScheme(index, ranges);
            
        // Populate it
        scheme->minContext = options["context"].as<size_t>();
        scheme->addContext = options["addContext"].as<size_t>();
        scheme->multContext = options["multContext"].as<double>();
        scheme->credit = options.count("credit");
        scheme->z_max = options["mismatches"].as<size_t>();
        
        mappingScheme = (MappingScheme*) scheme;
    } else if(options["mapType"].as<std::string>() == "natural") {
        // We want a NaturalMappingScheme
        NaturalMappingScheme* scheme = new NaturalMappingScheme(index, ranges);
            
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
        scheme->unstable = options.count("unstable");
        
        mappingScheme = (MappingScheme*) scheme;
    } else if(options["mapType"].as<std::string>() == "old") {
        // We want an OldNaturalMappingScheme
        OldNaturalMappingScheme* scheme = new OldNaturalMappingScheme(index,
            ranges);
            
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
        scheme->unstable = options.count("unstable");
        
        mappingScheme = (MappingScheme*) scheme;
    } else {
        // They asked for a mapping scheme we don't have.
        throw std::runtime_error("Invalid mapping scheme: " +
            options["mapType"].as<std::string>());
    }
    
    // Open the alignment file for writing
    std::ofstream alignment(options["alignment"].as<std::string>());
    
    // Now set up the parallel system we are going to use to map.
    
    // This holds records waiting to be mapped. Only one thread writes to it.
    ConcurrentQueue<std::pair<std::string, std::string>> recordQueue(1);
    
    // This holds output lines waiting to be written. All the mapping threads
    // write to it.
    ConcurrentQueue<std::string> lineQueue(numThreads);
    
    // This holds all our threads
    std::vector<Thread> threads;
    
    // Make a thread to load all the reads
    threads.push_back(Thread(&loadReads, std::ref(fastas), &recordQueue));
    
    for(size_t i = 0; i < numThreads; i++) {
        // Then some threads to process the reads
        threads.push_back(Thread(&mapSomeReads, &recordQueue, referenceRecord,
            std::ref(index), mappingScheme, &lineQueue));
    }
    
    // Then a thread to do the writing
    threads.push_back(Thread(&saveLines, &lineQueue, std::ref(alignment)));
    
    for(Thread& thread : threads) {
        // Wait for all the threads to be done
        thread.join();
    }
    
    if(options.count("stats")) {
        // Save statistics report to the specified file
        Log::info() << "Saving statistics to " <<
            options["stats"].as<std::string>() << std::endl;
        mappingScheme->stats.save(options["stats"].as<std::string>());
    }
    
    // Get rid of the mapping scheme now that everyone is done with it.
    delete mappingScheme;
    
    // Get rid of the index itself. Invalidates the index reference.
    delete indexPointer;

    // Now we're done!
    return 0;
}
