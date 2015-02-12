// Test the LR MappingScheme against a different FASTA file.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../MismatchResultSet.hpp"
#include "../util.hpp"

#include "LRMappingSchemeTests2.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( LRMappingSchemeTests2 );

// Define constants
const std::string LRMappingSchemeTests2::filename = "Test/haplotypes2.fa";

LRMappingSchemeTests2::LRMappingSchemeTests2() {
}

LRMappingSchemeTests2::~LRMappingSchemeTests2() {
}

void LRMappingSchemeTests2::setUp() {
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
    
    // Declare everything to be a range
    ranges = new GenericBitVector();
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        ranges->addBit(i);
    }
    ranges->finish(index->getBWTLength());

    // Make the mapping scheme
    scheme = new LRMappingScheme(*index, ranges);
}


void LRMappingSchemeTests2::tearDown() {
    // Clean up the mapping scheme
    delete scheme;
    
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
    
    // Delete the index
    delete index;
    
    // And the ranges bit vector
    delete ranges;
}


/**
 * Make sure we properly handle mismatches and let them produce ambiguity.
 */
void LRMappingSchemeTests2::testMismatch() {

    // We'll map this contig with only one distinctive base to itself.
    std::string query = "AAAAAAAAAAAAACAAAAAAAAAA";
    
    // Set for mismatches
    scheme->z_max = 1;
    
    size_t mappedBases = 0;
    
    // Map everything and get callbacks
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Only the C should ever map
        CPPUNIT_ASSERT_EQUAL('C', query[i]);
        mappedBases++;
    });
    
    // Exactly one C should map.
    CPPUNIT_ASSERT_EQUAL((size_t) 1, mappedBases);
}

/**
 * Make sure we key on a mapping on the first base, even when using mismatches.
 */
void LRMappingSchemeTests2::testMapOnMismatch() {

    // This string appears once exactly, but differentiated only by its last
    // base.
    std::string query = "AAAAAAAAAAAC";
    
    // Set for mismatches
    scheme->z_max = 1;
    
    size_t mappedBases = 0;
    
    // Map everything and get callbacks
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Only the C should ever map, and it will map on its left context.
        // TODO: We can't see that because we don't get a Mapping back. Check
        // that.
        CPPUNIT_ASSERT_EQUAL('C', query[i]);
        mappedBases++;
    });
    
    // Exactly one C should map.
    CPPUNIT_ASSERT_EQUAL((size_t) 1, mappedBases);
}










