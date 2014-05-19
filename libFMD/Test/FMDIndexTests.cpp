// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "FMDIndex.hpp"
#include "FMDIndexBuilder.hpp"
#include "util.hpp"

#include "FMDIndexTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( FMDIndexTests );

// Define constants
const std::string FMDIndexTests::filename = "test/haplotypes.fa";

void FMDIndexTests::setUp() {
    // We ened a built index as a fixture.

    // Set up a temporary directory to put the index in.
    tempDir = make_tempdir();
    
    // Make a new index builder.
    FMDIndexBuilder builder(tempDir + "/index.basename");
    
    // Add the haplotypes file
    builder.add(filename);
    
    // Finish the index.
    builder.close();
    
}


void FMDIndexTests::tearDown() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

/**
 * Test index metadata.
 */
void FMDIndexTests::testMetadata() {
    // Load the index up
    FMDIndex index(tempDir + "/index.basename");
    
    // Make sure it has the right number of characters.
    CPPUNIT_ASSERT(index.getTotalLength() == 35 * 2 * 2);
    
    // Make sure it has the right number of BWT positions (characters + texts).
    CPPUNIT_ASSERT(index.getBWTLength() == index.getTotalLength() + 4);
}

/**
 * Test the LF mapping.
 */
void FMDIndexTests::testLF() {

    // Load the index up
    FMDIndex index(tempDir + "/index.basename");

    // Make sure LF behaves sanely and maps the first $ to the first $.
    std::cout << "Index 21 maps LF to " << index.getLF(21) << std::endl;
    for(int i = 0; i < 21; i++) {
        CPPUNIT_ASSERT(index.display(i) != '$');
    }
    CPPUNIT_ASSERT(index.display(21) == '$');
    CPPUNIT_ASSERT(index.displayFirst(0) == '$');
    CPPUNIT_ASSERT(index.getLF(21) == 0);
    
    // Make sure it is consistent and always maps at least to an instance of the
    // correct character.
    for(int i = 0; i < index.getBWTLength(); i++) {
        // Get the first character of where we go, and compare it to the last
        // character of where we are.
        CPPUNIT_ASSERT(index.displayFirst(index.getLF(i)) == index.display(i));
    }

}

/**
 * Test dumping an FMD index's BWT.
 */
void FMDIndexTests::testDump() {
    
    // Load the index up
    FMDIndex index(tempDir + "/index.basename");
    
    // Dump the entire BWT
    for(int i = 0; i < index.getBWTLength(); i++) {
        // Reconstruct the string.
        std::string reconstruction;
        
        
        // Start at the last character of this row in the BWT.
        int64_t curIndex = i;
        
        do {
            // Put the current character on the reconstruction.
            reconstruction.push_back(index.display(curIndex));
            
            // Go to the previous character in the BWT row.
            curIndex = index.getLF(curIndex);
        } while(curIndex != i);      
    
        // Flip the string around forwards. See
        // <http://www.cplusplus.com/forum/beginner/11633/>
        std::string forwardReconstruction(reconstruction.rbegin(),
            reconstruction.rend());
    
        std::cout << i << ": " << index.displayFirst(i) << " " << 
            forwardReconstruction << " " << index.display(i) << std::endl;
    }
    
}

/**
 * Test searching an an FMD index.
 */
void FMDIndexTests::testSearch() {
    
    // Load the index up
    FMDIndex index(tempDir + "/index.basename");
    
    // Doesn't find spurious things
    CPPUNIT_ASSERT(index.count("GATTACA").getLength() == 0);
    
    // Finds things which appear once
    CPPUNIT_ASSERT(index.count("TCTTTT").getLength() == 1);
    
    // On both strands
    CPPUNIT_ASSERT(index.count("AAAAGA").getLength() == 1);
    
    // Finds things which appear twice
    CPPUNIT_ASSERT(index.count("TTCG").getLength() == 2);
    
    // Finds whole strands
    CPPUNIT_ASSERT(
        index.count("CGGGCGCATCGCTATTATTTCTTTCTCTTTTCACA").getLength() == 1);
    
}

/**
 * Test locating the positions of things you search up.
 */
void FMDIndexTests::testLocate() {

    // Load the index up
    FMDIndex index(tempDir + "/index.basename");
    
    // Find and locate a unique thing.
    TextPosition base = index.locate(index.count("TCTTTT").getForwardStart());
    
    // This is on the second sequence, forward strand.
    CPPUNIT_ASSERT(base.getText() == 2);
    // And it has 25 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 25);
    
    // Now let's try all of the first sequence, forwards
    base = index.locate(
        index.count("CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT").getForwardStart());
        
    // This is on the first sequence, forward strand.
    CPPUNIT_ASSERT(base.getText() == 0);
    // And it has 0 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 0);
    
    // And backwards
    base = index.locate(
        index.count("AGAGTCGCAGATGAGCGTCGAATCGCCGAAGCATG").getForwardStart());
        
    // This is on the first sequence, reverse strand.
    CPPUNIT_ASSERT(base.getText() == 1);
    // And it has 0 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 0);
}

/**
 * Test iterating over the suffix tree.
 */
void FMDIndexTests::testIterate() {

    // Load the index up
    FMDIndex index(tempDir + "/index.basename");

    for(int contextLength = 1; contextLength <= 25; contextLength++) {
        // Try all context lengths shorter than the contigs we put in.
        
        for(FMDIndex::iterator i = index.begin(contextLength); 
            i != index.end(contextLength); ++i) {
            // For each pair of suffix and position in the suffix tree
            
            // Unpack the iterator into pattern and FMDPosition at which it
            // happens.
            std::string pattern = (*i).first;
            // This is in SA coordinates.
            FMDPosition range = (*i).second;
        
            // Make sure we get the same range searching as iterating.
            CPPUNIT_ASSERT(index.count(pattern) == range);
        }
        
    }

}











