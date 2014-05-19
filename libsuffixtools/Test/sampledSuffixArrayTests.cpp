// Test the SampledSuffixArray's ability to build a Sampled Suffix Array.

#include "ReadTable.h"
#include "SuffixArray.h"
#include "SampledSuffixArray.h"

#include "sampledSuffixArrayTests.h"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( SampledSuffixArrayTests );

// Define constants
const std::string SampledSuffixArrayTests::filename = "Test/haplotypes.fa";

void SampledSuffixArrayTests::setUp() {
    // Set up the basic suffix array
    readTable = new ReadTable(filename);
    infoTable = new ReadInfoTable(filename);
    suffixArray = new SuffixArray(readTable, 1);
}


void SampledSuffixArrayTests::tearDown() {
    // Clean up
    delete suffixArray;
    delete infoTable;
    delete readTable;
}

/**
 * Test building a sampled suffix array.
 */
void SampledSuffixArrayTests::testConstruction() {
    
    // Make a BWT object from the normal suffix array
    BWT* bwt = new BWT(suffixArray, readTable);
    CPPUNIT_ASSERT(bwt != NULL);
    
    // Make a sampled suffix array
    SampledSuffixArray* sampled = new SampledSuffixArray();
    CPPUNIT_ASSERT(sampled != NULL);
    
    // Build it out from the BWT and the sequence info
    sampled->build(bwt, infoTable, 5);
    
    // Validate it
    sampled->validate(filename, bwt);
    
    // Practice saving it to files
    std::string devnull("/dev/null");
    sampled->writeSSA(devnull);
    sampled->writeLexicoIndex(devnull);
    
    delete sampled;
    delete bwt;
    
    
    
    

}
