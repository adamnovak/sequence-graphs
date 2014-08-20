/**
 * speedMap: tool for mapping things to a reference structure and seeing how
 * fast that is. Meant to be profiled.
 */

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
#include <sys/resource.h>
 
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

// Grab timers from libsuffixtools
#include <Timer.h>

#include <FMDIndex.hpp>
 
#include "Fasta.hpp"
 
/**
 * createIndex: command-line tool to create a multi-level reference structure.
 */
int main(int argc, char** argv) {

    // Parse options with boost::programOptions. See
    // <http://www.radmangames.com/programming/how-to-use-boost-program_options>

    std::string appDescription = 
        std::string("Map a sequence to a reference structure.\n") + 
        "Usage: speedMap <index directory> <fasta>";

    // Make an options description for our program's options.
    boost::program_options::options_description description("Options");
    // Add all the options
    description.add_options() 
        ("help", "Print help messages") 
        ("indexDirectory", boost::program_options::value<std::string>(), 
            "Directory to load the index from")
        ("fasta", boost::program_options::value<std::string>(),
            "FASTA file to map")
        ("repeat", boost::program_options::value<size_t>()
            ->default_value(1),
            "number of times to repeat each mapping");
            
        
    // And set up our positional arguments
    boost::program_options::positional_options_description positionals;
    // One index directory
    positionals.add("indexDirectory", 1);
    // And an one FASTA
    positionals.add("fasta", 1);
    
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
        
        if(!options.count("indexDirectory") || !options.count("fasta")) {
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
    
    // This holds the directory for the reference structure to read.
    std::string indexDirectory(options["indexDirectory"].as<std::string>());

    // This holds the FASTA we need to read.
    std::string fastaName(options["fasta"].as<std::string>());
    
    // This holds the number of repetitions to do
    size_t repetitions(options["repeat"].as<size_t>());
    
    // Load the index.
    FMDIndex index(indexDirectory + "/index.basename");
    
    // We want copies of all the mappings from the new method.
    std::vector<std::vector<Mapping>> newMappings;
    
    // And the old one
    std::vector<std::vector<Mapping>> oldMappings;
    
    // Load the FASTA
    Fasta fasta(fastaName);
    
    // Load all the records
    std::vector<std::string> nonemptyContigs;
    while(fasta.hasNext()) {
        std::string fastaRecord = fasta.getNext();
    
        // We're going to split each record into contigs.
        std::vector<std::string> recordContigs;
        
        // We compress runs of Ns thanks to Boost
        boost::algorithm::split(recordContigs, fastaRecord,
            boost::is_any_of("Nn"), boost::token_compress_on);
        
        // Put all the new contigs in the vector of contigs.
        nonemptyContigs.insert(nonemptyContigs.end(), recordContigs.begin(),
            recordContigs.end());
    }
    
    size_t originalContigs = nonemptyContigs.size();
    for(size_t i = 0; i < originalContigs; i++) {
        // Make sure each contig is in forwards and backwards, to eliminate
        // directional bias. TODO: not really needed.
        nonemptyContigs.push_back(reverseComplement(nonemptyContigs[i]));
    }
    
    // This holds reverse-complemented contigs, so each mapper can map in its
    // prefered direction.
    std::vector<std::string> rcContigs;
    for(auto record: nonemptyContigs) {
        rcContigs.push_back(reverseComplement(record));
    }
    
    // Get a stream so we can print without a new message every time
    auto output = Log::output();    
    
    // We want to time the mapping code.
    Timer* timer = new Timer("Mapping New Way");
    for(auto record : rcContigs) {
        for(size_t i = 0; i < repetitions; i++) {
            // Map the record a bunch
            output << ".";
            // TODO: This flips, so compare internal mapRight all with the
            // original map. This still introduces possibly a bias depending on
            // sequence orientation.
            newMappings.push_back(index.mapRight(record));
        }
    }
    output << std::endl;
    
    // Stop the timer
    delete timer;
    
    // Reverse all the new mappings (off the clock), so the will match up with
    // the right mappings
    for(size_t i = 0; i < newMappings.size(); i++) {
        // Put them in proper base order.
        std::reverse(newMappings[i].begin(), newMappings[i].end());
        
        for(size_t j = 0; j < newMappings[i].size(); j++) {
            // Go through all the mappings
            
            if(newMappings[i][j].is_mapped) {
                // Flip the mapping onto the correct text for left semantics.
                size_t contigLength = index.getContigLength(
                    index.getContigNumber(newMappings[i][j].location));
                    
                newMappings[i][j].location.setText(
                    newMappings[i][j].location.getText() ^ 1);
                newMappings[i][j].location.setOffset(contigLength - 
                    newMappings[i][j].location.getOffset() - 1);
            }
        }
    }
    
    
    // Get a new message
    output = Log::output();
    
    // And the old mapping code
    timer = new Timer("Mapping Old Way");
    
    for(auto record : nonemptyContigs) {
        for(size_t i = 0; i < repetitions; i++) {
            // Map the record a bunch
            output << ".";
            oldMappings.push_back(index.map(record));
        }
    }
    output << std::endl;
    
    // Stop the timer
    delete timer;
    
    // Keep some simple mapping stats
    size_t mapped = 0;
    size_t unmapped = 0;
    
    for(size_t i = 0; i < oldMappings.size(); i++) {
        // Check the answers to make sure they match and they did the same work.
        if(oldMappings[i].size() != newMappings[i].size()) {
            Log::critical() << "Mapping " << i << " has a length mismatch (" <<
            oldMappings[i].size() << " vs " << newMappings[i].size() << ")" <<
            std::endl;
            throw std::runtime_error("Mapping length mismatch");
        }
        
        for(size_t j = 0; j < oldMappings[i].size(); j++) {
            // Complain if they got different answers
            if(oldMappings[i][j] != newMappings[i][j]) {
                Log::error() << "Mapping mismatch: " << oldMappings[i][j] << 
                    " vs. " << newMappings[i][j] << std::endl;
            }
            
            if(oldMappings[i][j].is_mapped) {
                mapped++;
            } else {
                unmapped++;
            }
        }
        
    }
    
    Log::info() << mapped << " total mapped, " << unmapped << 
        " total unmapped" << std::endl;
    
}
    
    
    
    
    
    
    
    
    
    
    
