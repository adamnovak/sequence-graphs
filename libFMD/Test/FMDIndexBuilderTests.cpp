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

/**
 * Test the per-genome bitvectors. Belongs to the builder because it is
 * responsible for making these.
 */
void FMDIndexBuilderTests::testGenomes() {
    
    // Make a new index builder.
    FMDIndexBuilder builder(tempDir + "/index.basename");
    
    // Add the haplotypes file twice as two genomes
    builder.add(filename);
    builder.add(filename);
    
    // Finish the index.
    FMDIndex* index = builder.build();
    
    // Grab the mask for the first genome
    auto& genomeBitvector = index->getGenomeMask(0);
    
    // There were 35 bases in each contig, and we should count both the contigs
    // on both the strands, including their end of text characters.
    CPPUNIT_ASSERT(genomeBitvector.rank(genomeBitvector.getSize()) == 
        (35 + 1) * 2 * 2);
        
    // Grab the mask for the second genome
    auto& genomeBitvector2 = index->getGenomeMask(1);
    CPPUNIT_ASSERT(genomeBitvector.rank(genomeBitvector2.getSize()) == 
        (35 + 1) * 2 * 2);
        
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        // Make sure each BWT position is assigned to exactly one genome.
        CPPUNIT_ASSERT_MESSAGE(
            "One genome claims position " + std::to_string(i), 
            genomeBitvector.isSet(i) != genomeBitvector2.isSet(i));
    }
    
    // Don't leak the index
    delete index;
}
