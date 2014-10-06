// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "FMDIndexMismatchTest.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( FMDIndexMismatchTest );

// Define constants
const std::string FMDIndexMismatchTest::filename = "Test/haplotypes2.fa";

FMDIndexMismatchTest::FMDIndexMismatchTest() {

    // We need a built index as a fixture, and we don't want to rebuild it for
    // every test.

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

FMDIndexMismatchTest::~FMDIndexMismatchTest() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

void FMDIndexMismatchTest::setUp() {
    
}


void FMDIndexMismatchTest::tearDown() {
    
}


/**
 * Make sure we properly handle mismatches and let them produce ambiguity.
 */
void FMDIndexMismatchTest::testMismatch() {

    // We'll map this contig with only one distinctive base to itself.
    std::string query = "AAAAAAAAAAAAACAAAAAAAAAA";
    
    // Declare everything to be a range
    GenericBitVector bv;
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        bv.addBit(i);
    }
    bv.finish(index->getBWTLength());
    
    // Map it with right contexts, allowing for 1 mismatch.
    // TODO: can we specify things by name=value in c++11 to avoid breaking this
    // sneakily when we add new arguments?
    std::vector<Mapping> mappings = index->misMatchMap(bv,
        query, (GenericBitVector*)NULL, 0, 0, 0.0, 0.0, 1);
        
    for(size_t i = 0; i < query.size(); i++) {
        if(query[i] == 'A') {
            // None of the As should map, because under 1 mismatch the single C
            // can't place them.
            CPPUNIT_ASSERT(mappings[i].getRange() == -1);
        } else if(query[i] == 'C') {
            // The C should map.
            CPPUNIT_ASSERT(mappings[i].getRange() != -1);
        }
    }
}









