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
#include "Log.hpp"

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
    contigFile((basename + ".contigs").c_str()), genomeAssignments(),
    sampleRate(sampleRate) {

    // Nothing to do, already made everything.
    
}

void FMDIndexBuilder::add(const std::string& filename) {
    
    // Open the FASTA for reading.
    FILE* fasta = fopen(filename.c_str(), "r");
    
    if(fasta == NULL) {
        report_error("Failed to open FASTA " + filename);
    }
    
    // Work out what genome number this file gets. Either 0, or 1 more than the
    // last one used.
    size_t genomeNumber = (genomeAssignments.size() == 0) ? 0 : 
        genomeAssignments.back() + 1;
        
    // Start up the parser
    int fileNumber = fileno(fasta);
    kseq_t* seq = kseq_init(fileNumber); 
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
                    
                    // Pull it out
                    std::string run = sequence.substr(runStart, i - runStart);
                    
                    // Add the forward strand to the contig FASTA
                    tempFasta << ">" << name << "-" << runStart << "F" <<
                        std::endl;
                    tempFasta << run << std::endl;
                    
                    // And the reverse strand    
                    tempFasta << ">" << name << "-" << runStart << "R" <<
                        std::endl;
                    std::string reverseStrand = reverseComplement(run);
                    tempFasta << reverseStrand << std::endl;
                    
                    // Add the contig to the contig file (where we store FASTA
                    // record name, start, length, and genome). All contigs will
                    // be ordered by FASTA record they are from, and FASTA
                    // records from the same genome all appear together.
                    contigFile << name << "\t" << runStart << "\t" <<
                        (i - runStart) << "\t" << genomeNumber << std::endl;
                    
                    // Record that this sequence belongs to this genome.
                    genomeAssignments.push_back(genomeNumber);
                    
                }
                
                // The next run must start after here (or later).
                runStart = i + 1;
            }
        }
    }  
    kseq_destroy(seq); // Close down the parser.
}

FMDIndex* FMDIndexBuilder::build() {
    // TODO: Quiet this procedure down, or get logging down into this library or
    // something.

    // Close up the temp file
    tempFasta.close();
    
    // And the contig sizes file
    contigFile.close();
    
    // Compute what we want to save: BWT, sampled suffix array, and per-genome
    // BitVector masks.
    std::string bwtFile = basename + ".bwt";
    std::string ssaFile = basename + ".ssa";
    std::string bitmaskFile = basename + ".msk";

    // Produce the index of the temp file
    // Load all the sequences into memory (again).
    // TODO: Just keep them there
    ReadTable* readTable = new ReadTable(tempFastaName);
    
    Log::info() << "Computing index of " << tempFastaName << std::endl;
    
    // Compute the suffix array (which computes the BWT)
    SuffixArray* suffixArray = new SuffixArray(readTable, NUM_THREADS, true);
    
    Log::info() << "Saving BWT to " << bwtFile << std::endl;
    
    // Write the BWT to disk
    suffixArray->writeBWT(bwtFile, readTable);
    
    // Delete the read table since we no lonfger need it. Keep the suffix array
    // around because the FMDIndex we return can cheat off it.
    delete readTable;
    
    // How many genomes are there?
    size_t numGenomes = (genomeAssignments.size() == 0) ? 0 :
        genomeAssignments.back() + 1;
        
    Log::info() << "Creating " << numGenomes << " genome bitmasks..." <<
        std::endl;
    
    // Holds a bit vector encoder for each genome. TODO: Make this a C++11
    // vector with emplace_back to work around non-copy-constructability of
    // encoders. TODO: Before that, make this pointers in a vector.
    BitVectorEncoder** encoders = new BitVectorEncoder*[numGenomes];
    for(size_t i = 0; i < numGenomes; i++) {
        // Make each of the individual encoders.
        encoders[i] = new BitVectorEncoder(32);
    }
    
   
    for(size_t i = 0; i < suffixArray->getSize(); i++) {
        // Scan the suffix array, and make a 1 in the correct place in each bit
        // vector for each genome.
    
        // For each place in the SA
        
        // Get the SA element which stores text and offset.
        SAElem element = suffixArray->get(i);
        
        // Get the text ID
        size_t text = element.getID();
        
        // Compute the contig from the text. Easy since each contig has exactly
        // 2 texts.
        size_t contig = text / 2;
        
        // Look up the genome, and set this bit in the appropriate encoder.
        encoders[genomeAssignments[contig]]->addBit(i);
    }
    
    // Open the bitmask file
    std::ofstream bitmaskStream(bitmaskFile.c_str(), std::ios::binary);
    
    for(size_t i = 0; i < numGenomes; i++) {
        // Save all the bit vectors to the bitmask file.
        
        // Finish encoding to an actual BitVector
        encoders[i]->flush();
        BitVector bitVector(*encoders[i], suffixArray->getSize() + 1);
        
        // Save the BitVector
        bitVector.writeTo(bitmaskStream);
        
        // All the bitvectors can go in the same file. When reading them just
        // see if there are any bytes left. TODO: This is probably true.
        
        // Delete the encoder, sicne we've already encoded with it.
        delete encoders[i];
    }

    // Delete the encoders array altogether.
    delete[] encoders;
    
    // Finish up the bitmask file.
    bitmaskStream.flush();
    bitmaskStream.close();
    
    Log::info() << "Re-loading BWT..." << std::endl;
    
    // Load the BWT back in (instead of re-calculating it).
    // TODO: Add ability to save a calculated BWT object with a BWTWriter.
    BWT bwt(bwtFile);
    
    // Load all the sequence lengths from the original file.
    // TODO: just store these as we write the file.
    ReadInfoTable infoTable(tempFastaName);
    
    Log::info() << "Sampling suffix array..." << std::endl;
    
    // Make a sampled suffix array
    SampledSuffixArray sampled;
    
    // Build it from the BWT and read info, with the specified sample rate
    sampled.build(&bwt, &infoTable, sampleRate);
    
    Log::info() << "Saving sampled suffix array to " << ssaFile << std::endl;

    // Save it to disk    
    sampled.writeSSA(ssaFile);
    
    // Get rid of the temporary FASTA directory
    boost::filesystem::remove_all(tempDir);
    
    // Hand our SuffixArray off to an FMDIndex.
    return new FMDIndex(basename, suffixArray);
    
}


















