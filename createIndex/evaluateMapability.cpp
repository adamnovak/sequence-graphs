// evaluateMapability.cpp: program to see how much context is needed to map to a
// reference.
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
#include <regex>


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

// TODO: Matching stuff basically copied out of NaturalMappingScheme. Is there a
// good way to share it?

/**
 * Keep track of an occurrence of a unique string which matches a query
 * position to a reference position.
 */
struct Matching {
    /**
     * Where does it start in the query?
     */
    size_t start;
    /**
     * Where is it in the reference?
     */
    TextPosition location;
    /**
     * How long is it?
     */
    size_t length;
    
    /**
     * Make a new Matching.
     */
    inline Matching(size_t start, TextPosition location, size_t length): 
        start(start), location(location), length(length) {
    }
    
}

/**
 * Find all of the minimal unique matchings between query string characters
 * and the reference, in descending order by left endpoint.
 */
std::vector<Matching>
findMinMatchings(
    const FMDIndex& index,
    const std::string& query
) const {

    // What matchings have we found?
    std::vector<Matching> toReturn;
 
    // Start with everything selected.
    FMDPosition results = index.getCoveringPosition();
    
    // How many characters are currently searched?
    size_t patternLength = 0;

    // Flag that says whether we need to retract at least once before we can
    // find another minimal unique match, since we just found one ending at a
    // certain place. TODO: Find a cleaner way to do this.
    bool mustRetract = false;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--) {
        // For each position in the query from right to left, we're going to
        // consider any minimal unique matches with left endpoints here.
        
        // Retract on the right until we can successfully extend on the left
        // without running out of results.
        
        // We're going to extend backward with this new base.
        FMDPosition extended = results;
        index.extendLeftOnly(extended, query[i]);
        
        while(extended.isEmpty()) {
            // If you can't extend, retract until you can. TODO: Assumes we
            // can find at least one result for any character.
            
            // Retract the character
            FMDPosition retracted = results;
            // Make sure to drop characters from the total pattern length.
            index.retractRightOnly(retracted, --patternLength);
            mustRetract = false;
            
            // Try extending again
            extended = retracted;
            index.extendLeftOnly(extended, query[i]);
            
            // Say that last step we came from retracted.
            results = retracted;
        }
        
        // Extend on the left.
        results = extended;
        patternLength++;
        
        // Retract on the right until the next retraction would make us not
        // unique, and report a minimal unique match starting at this position.
        FMDPosition retracted = results;
        index.retractRightOnly(retracted, patternLength - 1);
        
        while(retracted.getLength() == 1) {
            // Retract until we would no longer be unique. Make sure to drop
            // characters from the total pattern length.
            results = retracted;
            index.retractRightOnly(retracted, (--patternLength) - 1);
            mustRetract = false;
        }
        
        if(results.getLength() == 1 && retracted.getLength() > 1 &&
            !mustRetract) {
            
            // We found a minimally unique match starting at this position and
            // ending patternLength right from here.
            toReturn.push_back(Matching(i,  index.locate(results.getResult()),
                patternLength));
                
            // We can't find another minimal match until we move the right
            // endpoint.
            mustRetract = true;
        }
    }
    
    // When we get here, we already know we retracted as much as we could for
    // the leftmost extension, so there are no more results to report.
    
    // This works: if there were a shorter unique match on the right starting at
    // this position, we would have retracted to find it. And if there were a
    // shorter unique match on the left ending at this position, we would have
    // already reported it and not reported any more until we retracted.
    
    // Results will also come out in descending order by left endpoint.
    
    return toReturn;
}

/**
 * evaluateMapability: command-line tool to evaluate how easy it is to map to
 * each base in a reference, in terms of minimum context length.
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
        std::string("Evaluate reference context lengths.\n") + 
        "Usage: evaluateMapability <index directory> <reference>";

    // Make an options description for our program's options.
    boost::program_options::options_description description("Options");
    // Add all the options
    description.add_options() 
        ("help", "Print help messages") 
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
        ("minHammingBound", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum Hamming distance lower bound on a maximum unique match");
        
    // And set up our positional arguments
    boost::program_options::positional_options_description positionals;
    // One index directory
    positionals.add("indexDirectory", 1);
    // One reference
    positionals.add("reference", 1);
    
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
        
        if(!options.count("indexDirectory") || !options.count("reference")) {
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
    
    // Make a vector of just the reference.
    std::vector<std::string> referenceOnly = { reference };
        
    // Index the reference. Use the sample rate the user specified.
    FMDIndex* indexPointer = buildIndex(indexDirectory, referenceOnly,
        options["sampleRate"].as<unsigned int>());
        
    // Make a reference out of the index pointer because we're not letting it
    // out of our scope.
    FMDIndex& index = *indexPointer;
    
    // Since we're mapping the reference to itself, there will be no mismatches.
    // All we have to do is load the reference, break it up on Ns, and for each
    // actual contig, find all the minimal unique matches. Then we need to do
    // another inchworm through each, counting up the minimum length required to
    // hit the needed alpha. Then we dump our histogram.
    
    // Make a Fasta to read the reference.
    Fasta referenceReader(reference);
    while(referenceReader.hasNext()) {
        // Get records until we run out
        std::string scaffold = referenceReader.getNext();
            
        std::for_each(std::sregex_token_iterator(scaffold.begin(),
            scaffold.end(), std::regex("N+"), -1), std::sregex_token_iterator(),
            [&](const std::string& contig) {
            
            // For every piece of sequence between runs of 1 or more N...
            
            // Find the minimum unique substrings between theis contig and the
            // whole reference.
            std::vector<Matching> minMatchings = findMinMatchings(index,
                contig);
                
            // Then we need to work out how far out we have to go to get a
            // certain minimum number of nonoverlapping ones.
            
        
        };
    }
    
    
    // Get rid of the mapping scheme now that everyone is done with it.
    delete mappingScheme;
    
    // Get rid of the index itself. Invalidates the index reference.
    delete indexPointer;

    // Now we're done!
    return 0;
}
