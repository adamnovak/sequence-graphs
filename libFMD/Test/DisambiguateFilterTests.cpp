// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "DisambiguateFilterTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( DisambiguateFilterTests );

// Define constants
const std::string DisambiguateFilterTests::filename = "Test/haplotypes.fa";

DisambiguateFilterTests::DisambiguateFilterTests() {

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

DisambiguateFilterTests::~DisambiguateFilterTests() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

void DisambiguateFilterTests::setUp() {
    
}


void DisambiguateFilterTests::tearDown() {
    
}

/**
 * Test applying the filter and disambiguating stuff.
 */
void DisambiguateFilterTests::testApply() {
    
    // Make some left mappings
    std::vector<Mapping> leftMappings {
        Mapping(TextPosition(0, 0)),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(0, 3)),
        Mapping(TextPosition(0, 4))
    };
    
    // Make some right mappings
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(1, 34)),
        Mapping(TextPosition(1, 33)),
        Mapping(),
        Mapping(TextPosition(1, 31)),
        Mapping(TextPosition(1, 29))
    };
    
    // Make the filter to test.
    DisambiguateFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
        
    // Check all the results
    // Both agree
    CPPUNIT_ASSERT_EQUAL(result[0], Mapping(TextPosition(1, 34)));
    // Left only
    CPPUNIT_ASSERT_EQUAL(result[1], Mapping(TextPosition(1, 33)));
    // Neither
    CPPUNIT_ASSERT_EQUAL(result[2], Mapping());
    // Right only
    CPPUNIT_ASSERT_EQUAL(result[3], Mapping(TextPosition(1, 31)));
    // Disagree
    CPPUNIT_ASSERT_EQUAL(result[4], Mapping());
    
}









