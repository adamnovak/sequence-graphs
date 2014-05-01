#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include "kseq.h"
#include "util.hpp"

#include "FMDIndexBuilder.hpp"

/**
 * Utility function to report last error and kill the program.
 */
void report_error(const std::string message) {

    // Something went wrong. We need to report this error.
    int errorNumber = errno;
    char* error = strerror(errno);
    
    // Complain about the error
    std::cerr << message << " (" << errorNumber << "): " <<
        error << std::endl;
    
    // Don't finish our program.
    throw std::runtime_error(message);
    
}

// Tell kseq what files are (handle numbers) and that you read them with read.
// Don't hook in .gz support. See <http://stackoverflow.com/a/19390915/402891>
KSEQ_INIT(int, read)

FMDIndexBuilder::FMDIndexBuilder(const std::string& basename, int sampleRate):
    basename(basename), sampleRate(sampleRate), input() {
    // Nothing to do. Already initialized everything.
    
}

void FMDIndexBuilder::add(const std::string& filename) {
    
    // Open the main index contig size list for appending
    std::ofstream contigStream;
    contigStream.open((basename + ".chrom.sizes").c_str(), std::ofstream::out |
        std::ofstream::app);
        
    // Open the FASTA for reading.
    FILE* fasta = fopen(filename.c_str(), "r");
    
    if(fasta == NULL) {
        report_error("Failed to open FASTA " + filename);
    }
    
    int fileNumber = fileno(fasta);
    
    kseq_t* seq = kseq_init(fileNumber); // Start up the parser
    while (kseq_read(seq) >= 0) { // Read sequences until we run out.
        // Stringify the sequence name
        std::string name(seq->name.s);
        
        std::cout << "Adding contig " << name << std::endl;
        
        // And the sequence sequence
        std::string sequence(seq->seq.s);
        
        // Upper-case all the letters. TODO: complain now if any not-base
        // characters are in the string.
        boost::to_upper(sequence);
        
        // Add the forward strand to the concatenated input. See
        // <http://stackoverflow.com/a/2551785/402891>
        input.insert(input.end(), sequence.begin(), sequence.end());
        
        // Null-terminate it.
        input.push_back(0);
        
        // Take the reverse complement and do the same
        std::string reverseComplement = reverse_complement(sequence);
        input.insert(input.end(), reverseComplement.begin(),
            reverseComplement.end());
        input.push_back(0);
        
        // Write the sequence ID and size to the contig list.
        contigStream << name << '\t' << sequence.size() << std::endl;
    }  
    kseq_destroy(seq); // Close down the parser.
    
    // Close up streams
    contigStream.flush();
    contigStream.close();

    

}

void FMDIndexBuilder::close() {
    
    std::cout << "Creating final index of " << input.size() << 
        " byte input..." << std::endl;
    
    // Make an RLCSA from out whole concatenated input. Don't try and delete the
    // data.
    CSA::RLCSA rlcsa((CSA::uchar*)&input[0], input.size(),
        CSA::RLCSA_BLOCK_SIZE.second, sampleRate, THREADS, false);
        
    // Throw out our input data. It made a copy.
    input.clear();
    
    // Pull out the RLCSA and save it.
    
    if(!(rlcsa.isOk())) {
        // Complain if it's broken.
        throw std::runtime_error("RLCSA integrity check failed!");
    }
    std::cout << "Saving RLCSA to " << basename << std::endl;
    
    // Dump its info.
    rlcsa.printInfo();
    rlcsa.reportSize(true);
    
    // Save it
    rlcsa.writeTo(basename);
}


















