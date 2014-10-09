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
#include <MarkovModel.hpp>
#include <Fasta.hpp>
#include <CreditFilter.hpp>
#include <DisambiguateFilter.hpp>

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
 * Start a new index in the given directory (by replacing it), and index the
 * given FASTAs for the bottom level FMD index. Optionally takes a suffix array
 * sample rate to use. Returns the FMD index that gets created.
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
        // These next three options should be ->required(), but that's not in the
        // Boost version I can convince our cluster admins to install. From now
        // on I shall work exclusively in Docker containers or something.
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
        ("credit", "Mapping on credit for greedy scheme")
        ("mismatches", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Maximum allowed number of mismatches");
        
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
    
    // Grab a bool for whether we'll map on credit    
    bool credit = options.count("credit");
    
    // Parse out the rest of the mapping scheme parameters
    size_t minContext = options["context"].as<size_t>();
    size_t addContext = options["addContext"].as<size_t>();
    size_t mismatches = options["mismatches"].as<size_t>();
    double multContext = options["multContext"].as<double>();
    
    // Open the alignment fie for writing
    std::ofstream alignment(options["alignment"].as<std::string>());
    
    // How many reads are mapped total?
    size_t totalReads = 0;
    
    // And how many mappings?
    size_t totalMappings = 0;
    
    for(std::string readFasta : fastas) {
        // Now go through all the read FASTAs
        Fasta reader(readFasta);
        while(reader.hasNext()) {
            // Go through each FASTA record
            std::pair<std::string, std::string> headerAndSequence = 
                reader.getNextRecord();
                
            // Split it out
            std::string recordName = headerAndSequence.first;
            std::string sequence = headerAndSequence.second;
            
            if(sequence.find('N') != std::string::npos) {
                // Skip this read
                Log::error() << "Skipping read with N: " << sequence << 
                    std::endl;
                continue;
            }
            
            totalReads++;
            
            // Map the sequence with credit
            
            // Map it on the right.
            std::vector<Mapping> rightMappings = index.misMatchMap(ranges,
                sequence, -1, minContext, addContext, multContext, 0,
                mismatches);
            
            // Map it on the left
            std::vector<Mapping> leftMappings = index.misMatchMap(ranges, 
                reverseComplement(sequence), -1, minContext, 
                addContext, multContext, 0, mismatches);
            
            // Flip the left mappings back into the original order. They should stay
            // as other-side ranges.
            std::reverse(leftMappings.begin(), leftMappings.end());
            
            for(size_t i = 0; i < leftMappings.size(); i++) {
                // Convert left and right mappings from ranges to base positions.
                
                if(leftMappings[i].isMapped()) {
                    // Locate by the BWT index we infer from the range number.
                    // We can do this since we made each BWT position its own
                    // range.
                    leftMappings[i].setLocation(index.locate(
                        leftMappings[i].getRange() - 1));
                }

                if(rightMappings[i].isMapped()) {
                    // Locate by the BWT index we infer from the range number.
                    rightMappings[i].setLocation(index.locate(
                        rightMappings[i].getRange() - 1));
                }
            }

            for(size_t i = 0; i < leftMappings.size(); i++) {
                // Convert all the left mapping positions to right semantics
                
                if(leftMappings[i].isMapped()) {
                    // Flip anything that's mapped, using the length of the contig it
                    // mapped to.
                    leftMappings[i] = leftMappings[i].flip(index.getContigLength(
                        leftMappings[i].getLocation().getContigNumber()));
                }
            }
             
            // Run the mappings through a filter to disambiguate and possibly
            // apply credit.
            std::vector<Mapping> filteredMappings;    
            if(credit) {
                // Apply a credit filter to the mappings
                filteredMappings = CreditFilter(index).apply(leftMappings,
                    rightMappings);
            } else {
                // Apply only a disambiguate filter to the mappings
                filteredMappings = DisambiguateFilter(index).apply(leftMappings,
                    rightMappings);
            }
            
            // Do one final pass to remove mappings of the wrong base, and
            // output. TODO: Make a real filter.
            for(size_t i = 0; i < sequence.size(); i++) {
                Mapping mapping = filteredMappings[i];
                
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
                    // Get the character it is mapped to
                    char mappedTo = referenceRecord.second[
                        mapping.getLocation().getOffset()]; 
                        
                    if(mapped != mappedTo) {
                        // Throw out this mapping, since it's placing a
                        // character on a character it doesn't match.
                        mapping = Mapping();
                        
                        Log::info() << "Tried mapping " << recordName << ":" <<
                            i << " (" << mapped << ") to " << 
                            referenceRecord.first << ":" << 
                            mapping.getLocation().getOffset() << " (" << 
                            mappedTo << ")" << std::endl;
                        
                    } else {
                    
                        // Now do the output for this line, because this position
                        // mapped. Do it as reference, then query.
                        alignment << referenceRecord.first << "\t" << 
                            mapping.getLocation().getOffset() << "\t" << 
                            recordName << "\t" << 
                            i << "\t" << 
                            backwards << std::endl;
                            
                        // Count that we had a mapping
                        totalMappings++;
                            
                    }
                    
                }
                
            }
            
        }
    }
    
    Log::output() << "Mapped " << totalReads << " total reads with " << 
        totalMappings << " mapped positions" << std::endl;
    
    // Get rid of the index itself. Invalidates the index reference.
    delete indexPointer;

    // Now we're done!
    return 0;
}
