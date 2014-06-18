// Test the BWT generation.

#include <boost/filesystem.hpp>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "FMDIndexBuilderTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( FMDIndexBuilderTests );

// Define constants
const std::string FMDIndexBuilderTests::filename = "Test/haplotypes.fa";

void FMDIndexBuilderTests::setUp() {
    // Set up a temporary directory to put the index in.
    tempDir = make_tempdir();
}


void FMDIndexBuilderTests::tearDown() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

/**
 * Test building an FMD index.
 */
void FMDIndexBuilderTests::testBuild() {
    
    // Make a new index builder.
    FMDIndexBuilder builder(tempDir + "/index.basename");
    
    // Add the haplotypes file
    builder.add(filename);
    
    // Finish the index.
    FMDIndex* index = builder.build();
    
    // Don't leak 
    delete index;
}
