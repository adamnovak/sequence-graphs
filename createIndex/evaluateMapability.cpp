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
    
};

/**
 * Find all of the minimal unique matchings between query string characters
 * and the reference, in descending order by left endpoint.
 */
std::vector<Matching>
findMinMatchings(
    const FMDIndex& index,
    const std::string& query
) {

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
        "Usage: evaluateMapability <index directory> <reference> <output file>";

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
            "FASTA contining a single reference sequence to map to")
        ("outputFile", boost::program_options::value<std::string>(),
            "File to save context lengths to, one per line")
        ("minHammingBound", boost::program_options::value<size_t>()
            ->default_value(0), 
            "Minimum Hamming distance lower bound on a maximum unique match");
            
        
    // And set up our positional arguments
    boost::program_options::positional_options_description positionals;
    // One index directory
    positionals.add("indexDirectory", 1);
    // One reference
    positionals.add("reference", 1);
    // One output file
    positionals.add("outputFile", 1);
    
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
            !options.count("outputFile")) {
            
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
    
    // How many min matches do we need in a run
    size_t minHammingBound = options["minHammingBound"].as<size_t>();
    
    // This holds the directory for the reference structure to build.
    std::string indexDirectory(options["indexDirectory"].as<std::string>());
    
    // This holds the reference
    std::string reference(options["reference"].as<std::string>());
    
    // This holds the filename to write to
    std::string outputFilename = options["outputFile"].as<std::string>();
    // Open that file for writing.
    std::ofstream outputFile(outputFilename.c_str());
    
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
            
            Log::info() << "Contig: " << contig << std::endl;
            
            // Find the minimum unique substrings between theis contig and the
            // whole reference.
            std::vector<Matching> minMatchings = findMinMatchings(index,
                contig);
                
            // Flip them around to be in ascending order by left endpoint.
            std::reverse(minMatchings.begin(), minMatchings.end());

            // We're going to fill in this vector of runs, which are [start,
            // end) pairs, in ascending order by start position. All runs are
            // right-minimal, and all left-minimal runs will be present.
            std::vector<std::pair<size_t, size_t>> runs;
            
            // We're going to do a pretty dumb approach where we fill in this
            // vector repeatedly until it holds the right min context length for
            // each base. Fill it with the biggest size_t.
            std::vector<size_t> minContextLengths(contig.size(), (size_t) -1);

            for(size_t i = 0; i < minMatchings.size(); i++) {
                // For each match, find the right-minimal run of matches we need
                // to get sufficient alpha.
                
                Log::info() << "Start at matching " << minMatchings[i].start <<
                    " + " << minMatchings[i].length << std::endl;
                
                // How many of them have we found? We start with 1.
                size_t nonOverlappingFound = 1;
                
                // The last non-overlapping thing we found was this first one.
                size_t last = i;
                
                for(size_t j = i + 1; j < minMatchings.size() && 
                    nonOverlappingFound < minHammingBound; j++) {
                    
                    // For each subsequent matching until we run out or get
                    // enough, see if it overlaps the last one we grabbed.
                    
                    Log::info() << "\tMatching " << minMatchings[j].start <<
                        " + " << minMatchings[j].length << std::endl;
                    
                    if(minMatchings[j].start >= minMatchings[last].start +
                        minMatchings[last].length) {
                        
                        Log::info() << "\t+++ Non-overlapping" << std::endl;
                    
                        // We know it has to start after the last one starts,
                        // and we can see the last one also ends before it
                        // starts. So it's non-overlapping.
                        
                        // Since min matches can't contain each other, if this
                        // is the first match that begins after our old match
                        // ends, it necessarily ends before any other such
                        // matches. So we can just take it and we automatically
                        // have the greedy activity selection algorithm.
                        
                        // Say this is now our last match.
                        last = j;
                        // And that we found another non-overlapping match.
                        nonOverlappingFound++;
                    } else {
                        Log::info() << "\t--- Overlapping" << std::endl;
                    }
                }
                
                if(nonOverlappingFound >= minHammingBound) {
                    // We managed to find enough things for our run.
                    
                    // Where does it start and end, in bases? Start inclusive,
                    // end exclusive.
                    size_t runStart = minMatchings[i].start;
                    size_t runEnd = minMatchings[last].start +
                        minMatchings[last].length;
                        
                    Log::info() << "Successful run " << runStart << " - " <<
                        runEnd << std::endl;
                        
                    // Save this run.
                    runs.push_back(std::make_pair(runStart, runEnd));
                }
            }
            
            if(runs.size() == 0) {
                // Complain we found no runs.
                throw std::runtime_error(
                    "No valid runs found, so no bases can map at all");
            }
            
            for(auto run: runs) {
                // For each run, put its length if it is the shortest thing
                // overlapping a base.
                
                // How long is the run?
                size_t runLength = run.second - run.first;
                
                for(size_t i = run.first; i < run.second; i++) {
                    // For each base in the run, put the run length if it hasn't
                    // gotten anything smaller.
                    minContextLengths[i] = std::min(minContextLengths[i],
                        runLength);
                }
            }
            
            // Now we need to scan from side to side for things that are off the
            // ends of runs.
            
            // We can guarantee that the endpoints of runs in runs do not move
            // backwards. For them to move backwards, an earlier run would have
            // to completely contain a later run. But the first matching of the
            // first run would have chained to whatever the first matching of
            // the later run chained to in order to finish sooner, so that can
            // never happen.
            
            // What run starts last and finishes before this base?
            size_t latestStarting = 0;
            for(size_t i = runs[0].second; i < minContextLengths.size(); i++) {
            
                // For each base after the end of the first run...
                
                while(latestStarting + 1 < runs.size() &&
                    // TODO: This one will always be true because runs don't
                    // share left min matches and min matches can't share
                    // endpoints.
                    runs[latestStarting + 1].first >
                    runs[latestStarting].first && 
                    runs[latestStarting + 1].second <= i) {
                    
                    // We can advance to a later-starting run that still doesn't
                    // cover this base. Do that.
                    latestStarting++;
                }
                
                // How long a string do we need to get out to the left end of
                // the latest starting range left of us?
                size_t contextLength = i - runs[latestStarting].first + 1;
                
                Log::info() << "Latest starting run before " << i << " is " <<
                    runs[latestStarting].first << " - " <<
                    runs[latestStarting].second << " with context length " <<
                    contextLength << " vs. " << minContextLengths[i] <<
                    std::endl;
                
                // Adopt this context length if it is minimal.
                minContextLengths[i] = std::min(minContextLengths[i],
                        contextLength);
            }
            
            // What run ends last and starts after this base?
            size_t earliestEnding = runs.size() - 1;
            for(size_t i = runs[runs.size() - 1].first - 1; i != (size_t) -1;
                i--) {
                
                // For each base before the beginning of the last run...
                
                while(earliestEnding - 1 != (size_t) -1 &&
                    // TODO: this one will always be true because run endpoints
                    // are nondecreasing left to right.
                    runs[earliestEnding - 1].second <=
                    runs[earliestEnding].second &&
                    runs[earliestEnding - 1].first > i) {
                    
                    // We can move to a lefter run while still not overlapping
                    // base i. Do that.
                    earliestEnding--;
                }
                
                // How long a string do we need to get out to the right end of
                // the earliest ending range right of us?
                size_t contextLength = runs[earliestEnding].second - i;
                
                Log::info() << "Earliest ending run after " << i << " is " <<
                    runs[earliestEnding].first << " - " <<
                    runs[earliestEnding].second << " with context length " <<
                    contextLength << " vs. " << minContextLengths[i] <<
                    std::endl;
                
                // Adopt this context length if it is minimal.
                minContextLengths[i] = std::min(minContextLengths[i],
                        contextLength);
                
            }
            
            Log::output() << minContextLengths.size() << " lengths" << std::endl;
            
            for(size_t contextLength : minContextLengths) {
                // Output the final length at each base, one per line.
                outputFile << contextLength << std::endl;
            }
        
        });
    }
    
    outputFile.close();
    
    // Get rid of the index itself. Invalidates the index reference.
    delete indexPointer;

    // Now we're done!
    return 0;
}
