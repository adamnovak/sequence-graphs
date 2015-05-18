#include <csignal>
#include <unordered_set>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "unixUtil.hpp"
#include "pinchGraphUtil.hpp"

#include <Fasta.hpp>
#include <Log.hpp>

#include <string>
#include <vector>
#include <set>

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
};

/**
 * Remove single quotes quoting the given string, if any. Only sees quotes if
 * they are the first and last character.
 */
std::string unquote(const std::string& input) {

    if(input.size() > 0 && input[0] == '\'' && 
        input[input.size() - 1] == '\'') {
        
        // It's quoted.
        return input.substr(1, input.size() - 2);
    }
    
    return input;

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
        "Usage: cactusMerge <c2hOut> <fastaOut> --c2h <c2h files...> "
            "--fasta <fasta files...> --suffix <suffixes...>");

    // Make an options description for our program's options.
    boost::program_options::options_description description("Options");
    // Add all the options
    description.add_options() 
        ("help", "Print help messages")
        ("c2h", boost::program_options::value<std::vector<std::string>>(),
            "List of c2h files to merge")
        ("fasta", boost::program_options::value<std::vector<std::string>>(),
            "List of FASTA files for the given c2h files")
        ("suffix", boost::program_options::value<std::vector<std::string>>(),
            "List of suffixes to add on to event names")
        ("mergeOn", boost::program_options::value<std::string>()->required(), 
            "An event on which to merge the files")
        ("c2hOut", boost::program_options::value<std::string>()->required(), 
            "File to save .c2h-format alignment in")
        ("fastaOut", boost::program_options::value<std::string>()->required(), 
            "File in which to save FASTA records for building HAL from .c2h");
        
        
        
    // And set up our positional arguments
    boost::program_options::positional_options_description positionals;
    positionals.add("mergeOn", 1);
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
        
        if(!options.count("mergeOn") || !options.count("c2h") ||
            !options.count("fasta")) {
            
            // We need both of these
            throw boost::program_options::error("Missing important arguments!");
        }
        
        if(options["c2h"].as<std::vector<std::string>>().size() != 
            options["fasta"].as<std::vector<std::string>>().size()) {
            
            // Counts need to match up here, because these are pairs
            throw boost::program_options::error(
                "c2h/fasta counts don't match!");
        }
        
        if(options.count("suffix") && 
            options["c2h"].as<std::vector<std::string>>().size() != 
            options["suffix"].as<std::vector<std::string>>().size()) {
        
            // If we have any suffixes we must have the right number
            throw boost::program_options::error(
                "c2h/suffix counts don't match!");
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
    std::vector<std::string> c2hFiles(
        options["c2h"].as<std::vector<std::string>>());
    
    // This holds the suffix applied to all the top sequences and events in each
    // file.
    std::vector<std::string> suffixes(
        options["suffix"].as<std::vector<std::string>>());
        
    // Make a list of the FASTA files to use
    std::vector<std::string> fastaFiles(
        options["fasta"].as<std::vector<std::string>>());
    
    // This will hold all of the renames that have to happen for each file.
    // These are generated when we go through the file by renaming top and
    // bottom sequences with suffixes.
    std::vector<std::map<std::string, std::string>> renames;
    
    for(size_t i = 0; i < c2hFiles.size(); i++) {
        // Make sure it has an empty map of renames for each file.
        renames.push_back(std::map<std::string, std::string>());
    }
    
    // This will hold the event names for the c2h files in order
    std::vector<std::string> eventNames;
    
    // And this will hold the sequence names
    std::vector<std::string> sequenceNames;
    
    // This will hold bottom (1) and top (0) flags for each sequence.
    std::vector<bool> isBottom;
    
    // And this will hold the sequence lengths
    std::vector<size_t> sequenceLengths;
    
    // This will hold the first sequence number for any event and sequence name
    std::map<std::pair<std::string, std::string>, size_t> firstSequenceNumber;
    
    // Holds Merge structs to be executed later.
    std::vector<C2hMerge> merges;
    
    // We're going to throw out all of the events that are old rootSeqs, and
    // just keep the actual leaves. This holds the list of renamed event names
    // we are keeping.
    std::set<std::string> eventsToKeep;
    
    for(size_t fileIndex = 0; fileIndex < c2hFiles.size(); fileIndex++) {
        // Scan through the c2h files to get the event, sequence, and length of
        // each thread, and to collect merges.
        
        Log::output() << "Reading alignment " << c2hFiles[fileIndex] << 
            std::endl;
        
        // Open the file
        std::ifstream c2h(c2hFiles[fileIndex]);
        
        // This maps block name to (sequence number, start location) pairs for
        // this file. We use it to compose merges for our list in global
        // sequence number space.
        std::map<size_t, std::pair<size_t, size_t>> nameMap;
        
        for(std::string line; std::getline(c2h, line);) {
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
                std::string eventName = unquote(parts[1]);
                std::string sequenceName = unquote(parts[2]);
                bool bottomFlag = std::stoi(parts[3]);
                
                Log::info() << "Read sequence " << eventName << "." << 
                    sequenceName << (bottomFlag ? " (bottom)" : " (top)") << 
                    std::endl;
                
                if(eventName != options["mergeOn"].as<std::string>()) {
                    // We aren't merging on this sequence, so we may have to
                    // apply a suffix.
                
                    // We need to rename this event (possibly to the same thing)
                    renames[fileIndex][eventName] = eventName + 
                        suffixes[fileIndex];
                        
                    if(bottomFlag) {
                        // All the bottom events (that aren't being merged
                        // on) need to be renamed apart manually since the names
                        // may be reused.
                        renames[fileIndex][eventName] += "-" + 
                            std::to_string(fileIndex);
                    }
                    eventName = renames[fileIndex][eventName];
                    
                    if(!bottomFlag) {
                        // Keep this event when we do our final output.
                        eventsToKeep.insert(eventName);
                    }
                    
                    // And the sequence
                    renames[fileIndex][sequenceName] = sequenceName + 
                        suffixes[fileIndex];
                        
                    if(bottomFlag) {
                        // All the bottom sequences (that aren't being merged
                        // on) need to be renamed apart manually since the names
                        // may be reused.
                        renames[fileIndex][sequenceName] += "-" + 
                            std::to_string(fileIndex);
                    }
                    sequenceName = renames[fileIndex][sequenceName];
                    
                    Log::info() << "Canonical name: " << eventName << "." << 
                        sequenceName << std::endl;
                    
                } else {
                    // If we are going to merge on it, we keep its name the same
                    // and then later we just make one thread for that name. We
                    // do definitely need it in the output though.
                    eventsToKeep.insert(eventName);
                }
                
                
                // Save the names
                eventNames.push_back(eventName);
                sequenceNames.push_back(sequenceName);
                
                // Save the bottomness flag.
                isBottom.push_back(bottomFlag);
                
                // Initialize the total length to 0
                sequenceLengths.push_back(0);
                
                
                auto namePair = std::make_pair(eventName, sequenceName);
                if(!firstSequenceNumber.count(namePair)) {
                    // This is the first time we have seen a sequence
                    // for this event and sequence name. Everything should
                    // merge against this one thread and not make more
                    // threads.
                    
                    // If this is the mergeOn event, we'll only make this
                    // once across all the files.
                    
                    // Later instances of this event and sequence name should
                    // redirect here.
                    firstSequenceNumber[namePair] = sequenceNames.size() - 1;
                    
                    Log::info() << "This is the first time we have seen "
                        "this sequence." << std::endl;
                }
                
            } else if(parts[0] == "a") {
                // This is an alignment block
                
                if(sequenceNames.size() == 0) {
                    throw std::runtime_error(
                        "Found alignmet block before sequence");
                }
                
                // Which sequence are we working on?
                size_t sequenceNumber = sequenceNames.size() - 1;
                
                if(isBottom[sequenceNumber]) {
                    // Parse it as a bottom block: "a" name start length
                    
                    if(parts.size() != 4) {
                        // Not the right number of fields.
                        throw std::runtime_error(
                            std::string("Invalid field count in ") + line);
                    }
                    
                    size_t blockName = std::stoll(parts[1]);
                    size_t blockStart = std::stoll(parts[2]);
                    size_t blockLength = std::stoll(parts[3]);
                    
                    // Look up the sequence number we actually want to merge
                    // against when we come get this block.
                    auto namePair = std::make_pair(eventNames[sequenceNumber], 
                        sequenceNames[sequenceNumber]);
                    size_t mergeSequenceNumber = firstSequenceNumber[namePair];
                    
                    // We need to associate the block name with the thread
                    // number for the sequence we want it to merge into, and the
                    // start location it specifies, for merging later.
                    nameMap[blockName] = std::make_pair(mergeSequenceNumber,
                        blockStart);
                    
                    Log::debug() << "Bottom block " << blockName << " is " << 
                        blockStart << " on sequence " << mergeSequenceNumber << 
                        std::endl;
                    
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
                    size_t segmentStart = std::stoll(parts[1]);
                    size_t segmentLength = std::stoll(parts[2]);
                    
                    // Add in the length
                    sequenceLengths[sequenceNumber] += segmentLength;
                    
                    if(parts.size() == 5) {
                        // If it has a name and orientation, remember a merge.
                        
                        size_t blockName = std::stoll(parts[3]);
                        bool orientation = std::stoi(parts[4]);
                        
                        // Get the sequence number that canonically represents
                        // all sequences with this event/sequence name
                        // combination.
                        auto namePair = std::make_pair(
                            eventNames[sequenceNumber], 
                            sequenceNames[sequenceNumber]);
                        size_t mergeSequenceNumber = firstSequenceNumber[
                            namePair];
                        
                        // Make a merge and populate it with everything we can
                        // get from this segment.
                        C2hMerge merge;
                        merge.sequence1 = mergeSequenceNumber;
                        merge.start1 = segmentStart;
                        merge.length = segmentLength;
                        // TODO: error-check length
                        merge.orientation = orientation;
                        
                        // Grab the info from the bottom segment we are talking
                        // about earlier in this file.
                        merge.sequence2 = nameMap[blockName].first;
                        merge.start2 = nameMap[blockName].second;
                        
                        Log::debug() << "Going to merge " << segmentStart << 
                            " length " << segmentLength << " to " <<
                            blockName << " orientation " << orientation << 
                            std::endl;
                        
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
        
        auto namePair = std::make_pair(eventNames[i], sequenceNames[i]);
        if(firstSequenceNumber[namePair] != i) {
            // This sequence is not the first; it is getting merged into another
            // one that is the same length and structure and name (i.e. it
            // appears in two files). Don't make a thread for it.
            continue;
        }
        
        
        // Make threads for all the top sequences and the first bottom sequence
        // for every event and sequence name pair.
        stPinchThreadSet_addThread(threadSet, i, 1, sequenceLengths[i]); 
    }
    
    for(auto merge : merges) {
        // Apply all the merges, converting merges to 1-based
        stPinchThread_pinch(
            stPinchThreadSet_getThread(threadSet, merge.sequence1),
            stPinchThreadSet_getThread(threadSet, merge.sequence2),
            merge.start1 + 1, merge.start2 + 1, merge.length, 
            merge.orientation);
            
        Log::trace() << "Applied merge between threads " << merge.sequence1 << 
            ":" << merge.start1 << "-" << merge.start1 + merge.length << 
            " and " << merge.sequence2 << ":" << merge.start2 << "-" << 
            merge.start2 + merge.length << " orientation " << 
            merge.orientation << std::endl;
    }
    
    // Write out a new c2h file, with a new rootSeq.
    size_t newRootLength = writeAlignment(threadSet, sequenceNames, eventNames, 
        options["c2hOut"].as<std::string>(), &eventsToKeep);
        
    // Clean up thread set.
    stPinchThreadSet_destruct(threadSet);
    
    // Merge the FASTAs, applying any renaming that needs to happen.
    
    // We'll do the FASTA output ourselves. Open the file.
    std::ofstream fastaOut(options["fastaOut"].as<std::string>());
    
    // Write the newly synthesized rootSeq. TODO: unify with writeAlignmentFasta
    // by moving support for renames over there.
    fastaOut << ">rootSeq" << std::endl;
    for(size_t i = 0; i < newRootLength; i++) {
        // Write an n for every base
        fastaOut << "N";
    }
    fastaOut << std::endl;
    
    // This holds the IDs of all the sequences we already wrote. Only write
    // sequences if they aren't duplicates after renaming (which is how we
    // deduplicate the shared root)
    std::unordered_set<std::string> alreadyWritten;
    
    for(size_t fileIndex = 0; fileIndex < fastaFiles.size(); fileIndex++) {
        // Open up the FASTA for reading
        Fasta fasta(fastaFiles[fileIndex]);
        
        Log::info() << "Copying over FASTA records from " <<
            fastaFiles[fileIndex] << std::endl;
        
        while(fasta.hasNext()) {
            // Go through all the FASTA records.
            // TODO: assumes FASTA headers have nothing but IDs.
            std::pair<std::string, std::string> record = fasta.getNextRecord();
            
            if(renames[fileIndex].count(record.first)) {
                // Rename them if necessary
                record.first = renames[fileIndex][record.first];
            }
            
            if(!eventsToKeep.count(record.first)) {
                // This event wasn't on the list of events to actually output,
                // so don't output it.
                Log::info() << "Skipped event " << record.first << std::endl;
                continue;
            }
            
            if(!alreadyWritten.count(record.first)) {
            
                // Save the record to the output FASTA file.
                fastaOut << ">" << record.first << std::endl << record.second <<
                    std::endl;
                
                // Remember that we have written a record by this name.
                alreadyWritten.insert(record.first);
            }
            
        }
    }
    
    fastaOut.close();
    
    // Now we're done!
    return 0;
}
