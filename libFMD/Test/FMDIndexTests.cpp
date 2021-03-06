// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "FMDIndexTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( FMDIndexTests );

// Define constants
const std::string FMDIndexTests::filename = "Test/haplotypes.fa";

FMDIndexTests::FMDIndexTests() {
}

FMDIndexTests::~FMDIndexTests() {
}

void FMDIndexTests::setUp() {

    // We need a built index as a fixture. It has to get built for every test.

    // Set up a temporary directory to put the index in.
    tempDir = make_tempdir();
    
    // Make a new index builder.
    FMDIndexBuilder builder(tempDir + "/index.basename");
    
    // Add the haplotypes file
    builder.add(filename);
    
    // Finish the index.
    FMDIndex* tmpIndex = builder.build();
    
    // Don't leak it.
    delete tmpIndex;
    
    // Save a pointer to a new index that we just load (so we don't have the
    // full SA).
    index = new FMDIndex(tempDir + "/index.basename");
    
}


void FMDIndexTests::tearDown() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

/**
 * Test index metadata.
 */
void FMDIndexTests::testMetadata() {
    
    // Make sure it has the right number of characters.
    CPPUNIT_ASSERT(index->getTotalLength() == 35 * 2 * 2);
    
    // Make sure it has the right number of BWT positions (characters + texts).
    CPPUNIT_ASSERT(index->getBWTLength() == index->getTotalLength() + 4);
    
    // Make sure it has the right number of contigs
    CPPUNIT_ASSERT(index->getNumberOfContigs() == 2);
    
    // Make sure it has the right number of genomes.
    CPPUNIT_ASSERT(index->getNumberOfGenomes() == 1);
    
    // Make sure the genome contains the contigs
    CPPUNIT_ASSERT(index->getGenomeContigs(0).first == 0);
    CPPUNIT_ASSERT(index->getGenomeContigs(0).second == 2);
    
    // Make sure the contigs belong to the genome
    CPPUNIT_ASSERT(index->getContigGenome(0) == 0);
    CPPUNIT_ASSERT(index->getContigGenome(1) == 0);
    
    // Make sure all the positions belong to the one genome
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        CPPUNIT_ASSERT_MESSAGE("Genome claims position " + std::to_string(i), 
            index->getGenomeMask(0).isSet(i));
    }
}

/**
 * Test the LF mapping.
 */
void FMDIndexTests::testLF() {

    // Make sure LF behaves sanely and maps the first $ to the first $.
    std::cout << "Index 21 maps LF to " << index->getLF(21) << std::endl;
    for(int i = 0; i < 21; i++) {
        CPPUNIT_ASSERT(index->display(i) != '$');
    }
    CPPUNIT_ASSERT(index->display(21) == '$');
    CPPUNIT_ASSERT(index->displayFirst(0) == '$');
    CPPUNIT_ASSERT(index->getLF(21) == 0);
    
    // Make sure it is consistent and always maps at least to an instance of the
    // correct character.
    for(int i = 0; i < index->getBWTLength(); i++) {
        // Get the first character of where we go, and compare it to the last
        // character of where we are.
        CPPUNIT_ASSERT(index->displayFirst(index->getLF(i)) == 
            index->display(i));
    }

}

/**
 * Test dumping an FMD index's BWT.
 */
void FMDIndexTests::testDump() {
    
    // Dump the entire BWT
    for(int i = 0; i < index->getBWTLength(); i++) {
        // Reconstruct the string.
        std::string reconstruction;
        
        
        // Start at the last character of this row in the BWT.
        int64_t curIndex = i;
        
        do {
            // Put the current character on the reconstruction.
            reconstruction.push_back(index->display(curIndex));
            
            // Go to the previous character in the BWT row.
            curIndex = index->getLF(curIndex);
        } while(curIndex != i);      
    
        // Flip the string around forwards. See
        // <http://www.cplusplus.com/forum/beginner/11633/>
        std::string forwardReconstruction(reconstruction.rbegin(),
            reconstruction.rend());
    
        std::cout << i << ": " << index->displayFirst(i) << " " << 
            forwardReconstruction << " " << index->display(i) << " " << 
            index->isInGenome(i, 0) << " " << index->getLCP(i) << std::endl;
    }
    
}

/**
 * Test pulling out a contig.
 */
void FMDIndexTests::testDisplay() {
    
    // Make sure the first contig comes out right
    CPPUNIT_ASSERT(index->displayContig(0) == 
        "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT");
        
    // And the second
    CPPUNIT_ASSERT(index->displayContig(1) == 
        "CGGGCGCATCGCTATTATTTCTTTCTCTTTTCACA");
        
    std::string contig1("CGGGCGCATCGCTATTATTTCTTTCTCTTTTCACA");
        
    for(size_t i = 0; i < contig1.size(); i++) {
        CPPUNIT_ASSERT_EQUAL(contig1[i], index->display(1, i + 1)); 
    }
}

/**
 * Test searching an an FMD index.
 */
void FMDIndexTests::testSearch() {
    
    // Doesn't find spurious things
    CPPUNIT_ASSERT_EQUAL((size_t) 0, index->count("GATTACA").getLength());
    
    // Finds things which appear once
    CPPUNIT_ASSERT_EQUAL((size_t) 1, index->count("TCTTTT").getLength());
    
    // On both strands
    CPPUNIT_ASSERT_EQUAL((size_t) 1, index->count("AAAAGA").getLength());
    
    // Finds things which appear twice
    CPPUNIT_ASSERT_EQUAL((size_t) 2, index->count("TTCG").getLength());
    
    // Finds whole strands
    CPPUNIT_ASSERT_EQUAL((size_t) 1, 
        index->count("CGGGCGCATCGCTATTATTTCTTTCTCTTTTCACA").getLength());
    
}

/**
 * Test locating the positions of things you search up.
 */
void FMDIndexTests::testLocate() {

    // Find and locate a unique thing.
    TextPosition base = index->locate(index->count("TCTTTT").getForwardStart());
    
    // This is on the second sequence, forward strand.
    CPPUNIT_ASSERT(base.getText() == 2);
    // And it has 25 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 25);
    
    // Now let's try all of the first sequence, forwards
    base = index->locate(
        index->count("CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT").getForwardStart());
        
    // This is on the first sequence, forward strand.
    CPPUNIT_ASSERT(base.getText() == 0);
    // And it has 0 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 0);
    
    // And backwards
    base = index->locate(
        index->count("AGAGTCGCAGATGAGCGTCGAATCGCCGAAGCATG").getForwardStart());
        
    // This is on the first sequence, reverse strand.
    CPPUNIT_ASSERT(base.getText() == 1);
    // And it has 0 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 0);
}

/**
 * Test iterating over the suffix tree.
 */
void FMDIndexTests::testIterate() {

    for(int contextLength = 1; contextLength <= 25; contextLength++) {
        // Try all context lengths shorter than the contigs we put in.
        
        for(FMDIndex::iterator i = index->begin(contextLength); 
            i != index->end(contextLength); ++i) {
            // For each pair of suffix and position in the suffix tree
            
            // Unpack the iterator into pattern and FMDPosition at which it
            // happens.
            std::string pattern = (*i).first;
            // This is in SA coordinates.
            FMDPosition range = (*i).second;
        
            // Make sure we get the same range searching as iterating.
            CPPUNIT_ASSERT(index->count(pattern) == range);
        }
        
    }

}

/**
 * Make sure the LCP array is not lying.
 */
void FMDIndexTests::testLCP() {

    // Save the previous reconstructed string at every step.
    std::string lastReconstruction;

    for(size_t i = 0; i < index->getBWTLength(); i++) {
        // Look at each LCP entry
        size_t lcp = index->getLCP(i);
        
        // Reconstruct the string.
        std::string reconstruction;

        // Start at the last character of this row in the BWT.
        int64_t curIndex = i;
        
        do {
            // Put the current character on the reconstruction.
            reconstruction.push_back(index->display(curIndex));
            
            // Go to the previous character in the BWT row.
            curIndex = index->getLF(curIndex);
        } while(curIndex != i);      
    
        // Flip the string around forwards. See
        // <http://www.cplusplus.com/forum/beginner/11633/>
        std::string forwardReconstruction(reconstruction.rbegin(),
            reconstruction.rend());
            
        
        if(i == 0) {
            // First entry must always be 0
            CPPUNIT_ASSERT(lcp == 0);
        } else {
            // This reconstruction has to mismatch the last one at the
            // appropriate index.
            
            size_t measuredLCP = 0;
            
            while(lastReconstruction[measuredLCP] == 
                forwardReconstruction[measuredLCP] && 
                lastReconstruction[measuredLCP] != '$') {
                
                // While they are characters that can match ("$" never matches
                // another string's "$")...
                measuredLCP++;
            }
            
            // Make sure they got what we got.
            CPPUNIT_ASSERT(measuredLCP == lcp);
            
        }
        
        // Save the reconstruction for the next test
        lastReconstruction = forwardReconstruction;
        
    }
    
}

/**
 * Test left-extending and right-retracting.
 */
void FMDIndexTests::testRetract() {

    // Select the whole thing.
    FMDPosition search = index->getCoveringPosition();
    
    // Search some stuff
    index->extendLeftOnly(search, 'C');
    index->extendLeftOnly(search, 'T');
    index->extendLeftOnly(search, 'T');
    
    // Make an interval to compare against
    FMDPosition countInterval = index->count("TTC");
    
    Log::output() << "For string TTC, got " << search << " by left extend, " <<
        countInterval << " by count." << std::endl;
    
    // Make sure they match in the forward-strand interval, which is the one
    // used here.
    CPPUNIT_ASSERT(search.getForwardStart() == countInterval.getForwardStart());
    CPPUNIT_ASSERT(search.getEndOffset() == countInterval.getEndOffset());
 
    // Reverse back 1 character (to a search string of 2).
    index->retractRightOnly(search, 2);
    
    // Search that instead
    countInterval = index->count("TT");
    
    Log::output() << "For string TT, got " << search << " by retract, " <<
        countInterval << " by count." << std::endl;
    
    // Compare
    CPPUNIT_ASSERT(search.getForwardStart() == countInterval.getForwardStart());
    CPPUNIT_ASSERT(search.getEndOffset() == countInterval.getEndOffset());
    
    // Reverse back to an empty search string.
    index->retractRightOnly(search, 0);
    
    // Search that
    countInterval = index->count("");
    
    Log::output() << "For empty string, got " << search << 
        " by left extend, " << countInterval << " by count." << std::endl;
    
    // Compare
    CPPUNIT_ASSERT(search.getForwardStart() == countInterval.getForwardStart());
    CPPUNIT_ASSERT(search.getEndOffset() == countInterval.getEndOffset());
    
}












