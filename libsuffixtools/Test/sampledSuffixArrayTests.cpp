// Test the SampledSuffixArray's ability to build a Sampled Suffix Array.

#include "Util/ReadTable.h"
#include "SuffixArray.cpp"

#include "sampledSuffixArrayTests.h"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( SampledSuffixArrayTests );

// Define constants
const std::string SampledSuffixArrayTests::filename = "Test/haplotypes.fa";

void SampledSuffixArrayTests::setUp() {
}


void SampledSuffixArrayTests::tearDown() {
}

/**
 * Test building a BWT. TODO: Replace!
 */
void SampledSuffixArrayTests::testConstruction() {
    
    // Make a read table from the headers in the file
	ReadTable* readTable = new ReadTable(filename);
    CPPUNIT_ASSERT(readTable != NULL);
    
    // Make a suffix array in 1 thread.
    SuffixArray* suffixArray = new SuffixArray(readTable, 1);
    CPPUNIT_ASSERT(suffixArray != NULL);
    
    // Validate it
    suffixArray->validate(readTable);
    
    // Try writing it out to nowhere.
    std::string devnull("/dev/null");
    suffixArray->writeBWT(devnull, readTable);
    suffixArray->writeIndex(devnull);
    
    // Clean up
    delete readTable;
    delete suffixArray;

}
