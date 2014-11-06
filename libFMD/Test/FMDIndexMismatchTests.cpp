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
        
    // Make sure nothing maps, the C is the last base and has no context to map
    // on.
    for(auto mapping : mappings) {
        CPPUNIT_ASSERT(mapping.getRange() == -1);
    }
    
    
    // Do the same thing in a different orientation, showing that it maps due to
    // a mismatch
    std::vector<Mapping> rcMappings = index->misMatchMap(bv,
        reverseComplement(query), (GenericBitVector*)NULL, 0, 0, 0.0, 0.0, 1);
        
    // Make sure the C maps, due to the C being forced to match exactly and
    // having some context.
    CPPUNIT_ASSERT(rcMappings[0].getRange() != -1);
    
    // And try again with a few bases to the right of the C, which should
    // differentiate it from the reverse complement of the G and let it map.
    query = "AAAAAAAAAAACAA";
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 0, 0, 0.0,
        0.0, 1);
    rcMappings = index->misMatchMap(bv, reverseComplement(query),
        (GenericBitVector*)NULL, 0, 0, 0.0, 0.0, 1);
        
    // The C should map
    CPPUNIT_ASSERT(mappings[11].getRange() != -1);
    // But things on either side of it shouldn't
    CPPUNIT_ASSERT(mappings[10].getRange() == -1);
    CPPUNIT_ASSERT(mappings[12].getRange() == -1);
    
    // And similarly for left mapping
    CPPUNIT_ASSERT(rcMappings[2].getRange() != -1);
    CPPUNIT_ASSERT(rcMappings[1].getRange() == -1);
    CPPUNIT_ASSERT(rcMappings[3].getRange() == -1);
    
}

/**
 * Test mismatch count.
 */
void FMDIndexMismatchTests::testMismatchCount() {

    // Declare everything to be a range
    GenericBitVector bv;
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        bv.addBit(i);
    }
    bv.finish(index->getBWTLength());

    // This string appears once within 1 mismatch, 0 times within 0.
    std::string query = "CAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAG";
    
    // No results in 0 mismatches
    MisMatchAttemptResults results = index->misMatchCount(bv, query, 0);
    CPPUNIT_ASSERT_EQUAL(results.positions.size(), (size_t)0);
    
    // 1 result in 1 mismatches
    results = index->misMatchCount(bv, query, 1);
    CPPUNIT_ASSERT_EQUAL(results.positions.size(), (size_t)1);
    CPPUNIT_ASSERT_EQUAL(results.is_mapped, true);
    
    // Many results in 1 mismatch for a shorter string.
    query = "CAAAAAAAAA";
    results = index->misMatchCount(bv, query, 1);
    CPPUNIT_ASSERT(results.positions.size() > 1);
    CPPUNIT_ASSERT_EQUAL(results.is_mapped, false);
    
}










