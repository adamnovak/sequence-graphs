// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <iostream>

#include "Util/ReadTable.h"
#include "SuffixArray.h"

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
 * Test dumping an FMD index's BWT.
 */
void FMDIndexTests::testDump() {
    
    // Load the index up
    FMDIndex index(tempDir + "/index.basename");
    
    // Make sure it has the right number of characters.
    CPPUNIT_ASSERT(index.getTotalLength() == 35 * 2 * 2);
    
    // Make sure it has the right number of BWT positions (characters + texts).
    CPPUNIT_ASSERT(index.getBWTLength() == index.getTotalLength() + 4);
    
    // Dump it
    for(int i = 0; i < index.getBWTLength(); i++) {
        std::cout << i << ": " << index.displayFirst(i) << " ... " << index.display(i) << std::endl;
    }
    
}
