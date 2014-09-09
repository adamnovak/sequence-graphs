#include <csignal>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "unixUtil.hpp"
#include "pinchGraphUtil.hpp"

/**
 * Represents a merge to be executed later after we build the pinch graph and
 * know how long all the threads are.
 */
struct C2hMerge {
    size_t sequence1;
    size_t start1;
    size_t sequence2;
    size_t start2;
    size_t length;
    bool orientation;
}

/**
 * cactusMerge.cpp: merge two pairs of c2h and FASTA files into one pair. The
 * files must be star trees with the same root sequence.
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
        std::string("Merge c2h/FASTA file pairs.\n" 
        "Usage: cactusMerge <c2h1> <fasta1> <c2h2> <fasta2> "
        "<c2hOut> <fastaOut>");

    // Make an options description for our program's options.
    boost::program_options::options_description description("Options");
    // Add all the options
    description.add_options() 
        ("help", "Print help messages")
        ("c2h1", boost::program_options::value<std::string>()->required(),
            "First c2h file to merge")
        ("fasta1", boost::program_options::value<std::string>()->required(),
            "FASTA for first c2h file")
        ("c2h2", boost::program_options::value<std::string>()->required(),
            "Second c2h file to merge")
        ("fasta2", boost::program_options::value<std::string>()->required(),
            "FASTA for second c2h file")
        ("c2hOut", boost::program_options::value<std::string>()->required(), 
            "File to save .c2h-format alignment in")
        ("fastaOut", boost::program_options::value<std::string>()->required(), 
            "File in which to save FASTA records for building HAL from .c2h");
        
        
    // And set up our positional arguments
    boost::program_options::positional_options_description positionals;
    positionals.add("c2h1", 1);
    positionals.add("fasta1", 1);
    positionals.add("c2h2", 1);
    positionals.add("fasta2", 1);
    positionals.add("c2hOut", 1);
    positionals.add("fastaOut", 1);
    
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
    
    // If we get here, we have the right arguments.
    
    // Make a list of the c2h files to use
    std::vector<std::string> c2hFiles {
        options["c2h1"].as<std::string>(),
        options["c2h1"].as<std::string>()
    };
    
    // This will hold the event names for the two c2h files in order
    std::vector<std::string> eventNames;
    
    // And this will hold the sequence names
    std::vector<std::string> sequenceNames;
    
    // This will hold bottom (1) and top (0) flags for each sequence.
    std::vector<bool> isBottom;
    
    // And this will hold the sequence lengths
    std::vector<size_t> sequenceLengths;
    
    // This will hold the first bottom sequence for any event and sequence name
    std::map<std::pair<std::string, sts::string>, size_t> firstBottomSequence;
    
    // Holds Merge structs to be executed later.
    std::vector<C2hMerge> merges;
    
    for(auto c2hFile : c2hFiles) {
        // Scan through the c2h files to get the event, sequence, and length of
        // each thread, and to collect merges.
        
        // Open the file
        std::istream c2h(c2hFile);
        
        // This maps block name to (sequence number, start location) pairs for
        // this file. We use it to compose merges for our list in global
        // sequence number space.
        std::map<size_t, std::pair<size_t, size_t>> nameMap;
        
        for(std::string line; std::getline(c2h, line)) {
            // This is a new sequence. Split it up on \t.
            std::vector<std::string> parts;
            boost::split(parts, line, boost::is_any_of("\t"));
        
            if(parts.size() < 1) {
                // Skip lines that have nothing on them.
                continue;
            }
        
            // For each line
            if(parts[0] == "s") {
                
                // It's a sequence line. Start a new squence.
                
                if(parts.size() != 4) {
                    // Not the right number of fields.
                    throw std::runtime_error(
                        std::string("Invalid field count in ") + line);
                }
                
                // Grab the parts
                std::string eventName = parts[1];
                std::string sequenceName = parts[2];
                bool bottomFlag = std::stoi(parts[3]);
                
                // Save the names
                eventNames.push_back(eventName);
                sequenceNames.push_back(sequenceName);
                
                // Save the bottomness flag.
                isBottom.push_back(bottomFlag);
                
                // Initialize the total length to 0
                sequenceLengths.push_back(0);
                
                if(bottomFlag) {
                    // It's a bottom sequence and we have to worry about
                    // merging.
                
                    auto namePair = std::make_pair(eventName, sequenceName);
                    if(!firstBottomSequence.count(namePair)) {
                        // This is the first time we have seen a bottom sequence
                        // for this event and sequence name. Everything should
                        // merge against this one thread and not make more
                        // threads.
                        
                        // Later instances of this event and sequence name as a
                        // bottom sequence should redirect here.
                        firstBottomSequence[namePair] = isBottom.size() - 1;
                    }
                    
                    // Otherwise we will count up length and stuff as normal,
                    // but not merge against this thread. We still need to look
                    // at all the blocks so we can fugure out where on the first
                    // thread we should merge when people ask for them.
                    
                }
                
            } else if(parts[0] == "a") {
                // This is an alignment block
                
                if(isBottom.size() == 0) {
                    throw std::runtime_error(
                        "Found alignmet block before sequence");
                }
                
                // Which sequence are we working on?
                size_t sequenceNumber = isBottom.size() - 1;
                
                if(isBottom[sequenceNumber]) {
                    // Parse it as a bottom block: "a" name start length
                    
                    if(parts.size() != 4) {
                        // Not the right number of fields.
                        throw std::runtime_error(
                            std::string("Invalid field count in ") + line);
                    }
                    
                    std::string blockName = parts[1];
                    size_t blockStart = std::stoi(parts[2]);
                    size_t blockLength = std::stoi(parts[3]);
                    
                    // Look up the sequence number we actually want to merge
                    // against when we come get this block.
                    auto namePair = std::make_pair(eventName, sequenceName);
                    size_t mergeSequenceNumber = firstBottomSequence[namePair];
                    
                    // We need to associate the block name with the thread
                    // number for the sequence we want it to merge into, and the
                    // start location it specifies, for merging later.
                    nameMap[blockName] = std::make_pair(mergeSequenceNumber,
                        blockStart);
                    
                    // Also record the additional length on this sequence
                    sequenceLengths[sequenceNumber] += blockLength;
                    
                } else {
                    // Parse it as a top block: 
                    // "a" start length [name orientation]
                    
                    if(parts.size() < 3) {
                        // Not the right number of fields.
                        throw std::runtime_error(
                            std::string("Invalid field count in ") + line);
                    }
                    
                    // Parse out the start and length
                    size_t segmentStart = std::stoi(parts[2]);
                    size_t segmentLength = std::stoi(parts[3]);
                    
                    // Add in the length
                    sequenceLengths[sequenceNumber] += segmentLength;
                    
                    if(parts.size() == 5) {
                        // If it has a name and orientation, remember a merge.
                        
                        size_t blockName = std::stoi(parts[4]);
                        bool orientation = std::stoi(parts[5]);
                        
                        // Make a merge and populate it with everything we can
                        // get from this segment.
                        C2hMerge merge;
                        merge.sequence1 = sequenceNumber;
                        merge.start1 = segmentStart;
                        merge.length = segmentLength;
                        // TODO: error-check length
                        merge.orientation = orientation;
                        
                        // Grab the info from the bottom segment we are talking
                        // about earlier in this file.
                        merge.sequence2 = nameMap[blockName].first;
                        merge.start2 = nameMap[blockName].first;
                        
                        // Save the merge for doing later.
                        merges.push_back(merge);
                        
                    }
                    
                }
            }
        }
    }
    
    // Make a thread set with all those threads
    stPinchThreadSet* threadSet = stPinchThreadSet_construct();
    
    // Make all the threads. Be 1-based internally since the serialization code
    // wants that.
    for(size_t i = 0; i < sequenceLengths.size(); i++) {
        if(isBottom[i]) {
            auto namePair = std::make_pair(eventNames[i], sequenceNames[i]);
            if(firstBottomSequence[namePair] != i) {
                // This bottom sequence is not the first; it is getting merged
                // into another one. Don't make a thread for it.
                continue;
            }
        }
        
        // Make threads for all the top sequences and the first bottom sequence
        // for every event and sequence name pair.
        stPinchThreadSet_addThread(threadSet, i, 1, sequenceLengths[i]); 
    }
    
    for(auto merge : merges) {
        // Apply all the merges, converting merges to 1-based
        stPinchThread_pinch(stPinchThreadSet_getThread(merge.sequence1), 
            stPinchThreadSet_getThread(merge.sequence2), merge.start1 + 1, 
            merge.start2 + 1, merge.length, merge.orientation);
    }
    
    // Write out a new c2h file, with the 0th thread as the reference.
    writeAlignmentWithReference(threadSet, sequenceNames, eventNames, 
        options["c2hOut"].as<std::string>(), 0);
    
    // Merge the FASTAs, applying any renaming that needs to happen.
    
    // Now we're done!
    return 0;
}
