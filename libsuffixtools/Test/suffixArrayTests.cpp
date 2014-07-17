// Test the SuffixArray.

#include "../ReadTable.h"
#include "../SuffixArray.h"

#include "suffixArrayTests.h"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( SuffixArrayTests );

// Define constants
const std::string SuffixArrayTests::filename = "Test/haplotypes.fa";

void SuffixArrayTests::setUp() {
}


void SuffixArrayTests::tearDown() {
}

/**
 * Test building a BWT.
 */
void SuffixArrayTests::testBWT() {
    
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
