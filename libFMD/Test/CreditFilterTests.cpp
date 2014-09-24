// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "CreditFilterTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( CreditFilterTests );

// Define constants
const std::string CreditFilterTests::filename = "Test/haplotypes.fa";

CreditFilterTests::CreditFilterTests() {

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

CreditFilterTests::~CreditFilterTests() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

void CreditFilterTests::setUp() {
    
}


void CreditFilterTests::tearDown() {
    
}

/**
 * Make sure credit works in an easy case where it should happen.
 */
void CreditFilterTests::testApply() {
    
    // Make some left mappings
    std::vector<Mapping> leftMappings {
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(0, 4), 5)
    };
    
    // Make some right mappings
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(1, 34), 5),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping()
    };
    
    // Make the filter to test.
    CreditFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
        
    for(auto x : result) {
        Log::info() << x << std::endl;
    }
        
    // Check all the results
    CPPUNIT_ASSERT_EQUAL(TextPosition(0, 0), result[0].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(0, 1), result[1].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(0, 2), result[2].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(0, 3), result[3].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(0, 4), result[3].getLocation());
    
}

/**
 * Make sure disagreement at bases that wold provide credit prevents credit.
 */
void CreditFilterTests::testDisagreement() {
    
    // Make some left mappings
    std::vector<Mapping> leftMappings {
        Mapping(TextPosition(0, 0)),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(2, 4), 5)
    };
    
    // Make some right mappings
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(3, 34), 5),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(1, 30))
    };
    
    // Make the filter to test.
    CreditFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
        
    // Check all the results
    for(size_t i = 0; i < 5; i++) {
        // Nothing should map, because neither of the bases that could provide
        // credit mapped, because there was disagreement at each.
        CPPUNIT_ASSERT_EQUAL(false, result[i].isMapped());
    }
    
}









