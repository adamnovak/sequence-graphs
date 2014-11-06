// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "FMDIndexMismatchTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( FMDIndexMismatchTests );

// Define constants
const std::string FMDIndexMismatchTests::filename = "Test/haplotypes2.fa";

FMDIndexMismatchTests::FMDIndexMismatchTests() {

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

FMDIndexMismatchTests::~FMDIndexMismatchTests() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

void FMDIndexMismatchTests::setUp() {
    
}


void FMDIndexMismatchTests::tearDown() {
    
}


/**
 * Make sure we properly handle mismatches and let them produce ambiguity.
 */
void FMDIndexMismatchTests::testMismatch() {

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

/**
 * Make sure we key in a mapping on the first base, even when using mismatches.
 */
void FMDIndexMismatchTests::testMapOnMismatch() {

    // This string appears once exactly, but differentiated only by its last
    // base.
    std::string query = "AAAAAAAAAAAC";
    
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
        
    // Make sure it doesn't map, since the C is the last base and not forced to
    // match exactly.
    CPPUNIT_ASSERT(mappings[11].getRange() == -1);
    
    
    // Do the same thing in a different orientation, showing that it maps due to
    // a mismatch
    std::vector<Mapping> rcMappings = index->misMatchMap(bv,
        reverseComplement(query), (GenericBitVector*)NULL, 0, 0, 0.0, 0.0, 1);
        
    // Make sure it maps, due to the C being forced to match exactly.
    CPPUNIT_ASSERT(rcMappings[0].getRange() != -1);
    
    // And try again with bases to the right of the C, which should not matter.
    query = "AAAAAAAAAAACAAA";
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 0, 0, 0.0,
        0.0, 1);
    rcMappings = index->misMatchMap(bv, reverseComplement(query),
        (GenericBitVector*)NULL, 0, 0, 0.0, 0.0, 1);
        
    CPPUNIT_ASSERT(mappings[11].getRange() == -1);
    CPPUNIT_ASSERT(rcMappings[3].getRange() != -1);
    
}










