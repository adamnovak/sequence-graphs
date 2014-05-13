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

#include <SampledSuffixArray.h>
#include <SuffixArray.h>
#include <ReadInfoTable.h>
#include <ReadTable.h>
#include <BWT.h>

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
    basename(basename), tempDir(make_tempdir()), 
    tempFastaName(tempDir + "/temp.fa"), tempFasta(tempFastaName.c_str()), 
    contigFile((basename + ".chrom.sizes").c_str()), sampleRate(sampleRate) {

    // Nothing to do, already made everything.
    
}

void FMDIndexBuilder::add(const std::string& filename) {
    
    std::cout << "Extracting contiguous runs from " << filename << std::endl;
        
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
        
        // And the sequence sequence
        std::string sequence(seq->seq.s);
        
        // Upper-case all the letters. TODO: complain now if any not-base
        // characters are in the string.
        boost::to_upper(sequence);

        // Where does the next run of not-N characters start?        
        size_t runStart = 0;
        // Iterate over contiguous runs of not-N
        for(size_t i = 0; i <= sequence.size(); i++) {
            if(i == sequence.size() || sequence[i] == 'N') {
                // This position is after the end of a run of not-N characters.
                
                if(i > runStart) {
                    // That run is nonempty. Process it.
                    
                    std::cout << "Run of " << (i - runStart) <<
                        " characters on " << name << std::endl;
                    
                    // Pull it out
                    std::string run = sequence.substr(runStart, i - runStart);
                    
                    // Add the forward strand to the contig FASTA
                    tempFasta << ">" << name << "-" << runStart << "F" <<
                        std::endl;
                    tempFasta << run << std::endl;
                    
                    // And the contig file
                    contigFile << name << "-" << runStart << "\t" <<
                        (i - runStart) << std::endl;
                    
                    // And the reverse strand    
                    tempFasta << ">" << name << "-" << runStart << "R" <<
                        std::endl;
                    std::string reverseStrand = reverseComplement(run);
                    tempFasta << reverseStrand << std::endl;
                    
                }
                
                // The next run must start after here (or later).
                runStart = i + 1;
            }
        }
    }  
    kseq_destroy(seq); // Close down the parser.
}

void FMDIndexBuilder::close() {
    // Close up the temp file
    tempFasta.close();
    
    // And the contig sizes file
    contigFile.close();
    
    // Compute what we want to save: BWT and sampled suffix array
    std::string bwtFile = basename + ".bwt";
    std::string ssaFile = basename + ".ssa";    

    std::cout << "Loading reads..." << std::endl;
    
    // Produce the index of the temp file
    // Load all the sequences into memory (again).
    // TODO: Just keep them there
    ReadTable* readTable = new ReadTable(tempFastaName);
    
    std::cout << "Computing index..." << std::endl;
    
    // Compute the suffix array (which computes the BWT)
    SuffixArray* suffixArray = new SuffixArray(readTable, NUM_THREADS);
    
    std::cout << "Saving BWT to " << bwtFile << std::endl;
    
    // Write the BWT to disk
    suffixArray->writeBWT(bwtFile, readTable);
    
    // Delete everything we no longer need.
    delete suffixArray;
    delete readTable;
    
    std::cout << "Loading BWT..." << std::endl;
    
    // Load the BWT back in (instead of re-calculating it).
    // TODO: Add ability to save a calculated BWT object with a BWTWriter.
    BWT bwt(bwtFile);
    
    std::cout << "Scanning reads..." << std::endl;
    
    // Load all the sequence lengths from the original file.
    // TODO: just store these as we write the file.
    ReadInfoTable infoTable(tempFastaName);
    
    std::cout << "Sampling suffix array..." << std::endl;
    
    // Make a sampled suffix array
    SampledSuffixArray sampled;
    
    // Build it from the BWT and read info, with the specified sample rate
    sampled.build(&bwt, &infoTable, sampleRate);
    
    std::cout << "Saving sampled suffix array to " << ssaFile << std::endl;

    // Save it to disk    
    sampled.writeSSA(ssaFile);
    
    // Get rid of the temporary FASTA directory
    boost::filesystem::remove_all(tempDir);
    
}


















