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
            "FASTA file to map");
        
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
    
    // Load the index.
    FMDIndex index(indexDirectory + "/index.basename");
    
    // Load the FASTA
    Fasta fasta(fastaName);
    
    while(fasta.hasNext()) {
        // Grab each FASTA record
        std::string fastaRecord(fasta.getNext());
        
        // Only map to the whole index (bottom level).
        index.mapBoth(fastaRecord);
        
        Log::output() << "Mapped FASTA record:" << std::endl;
        Log::output() << fastaRecord << std::endl;
    }
    
    
    
    
    
    
    
    
    
    
    
    
}
    
    
    
    
    
    
    
    
    
    
    
