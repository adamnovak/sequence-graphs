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
    
    // Make some left mappings. Make sure we have a sentinel first.
    std::vector<Mapping> leftMappings {
        Mapping(TextPosition(0, 0)),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(0, 4), 5)
    };
    
    // Make some right mappings. Make sure we have a sentinel last.
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(1, 34), 5),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(1, 30))
    };
    
    // Make the filter to test.
    CreditFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
    
    for(auto x : result) {
        Log::info() << x << std::endl;
    }
    
    // Make sure the result is the right length.
    CPPUNIT_ASSERT_EQUAL((size_t)5, result.size());
    
    for(size_t i = 0; i < 5; i++) {
        // Everythign should map
        CPPUNIT_ASSERT(result[i].isMapped());
    }
        
    // Check all the results
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 34), result[0].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 33), result[1].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 32), result[2].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 31), result[3].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 30), result[4].getLocation());
    
}

/**
 * Make sure disagreement at bases that wold provide credit prevents credit.
 */
void CreditFilterTests::testDisagreement() {
    
    // Make some left mappings
    std::vector<Mapping> leftMappings {
        Mapping(TextPosition(1, 34)),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(3, 30), 5)
    };
    
    // Make some right mappings
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(2, 0), 5),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(0, 4))
    };
    
    // Make the filter to test.
    CreditFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
    
    // Make sure the result is the right length.
    CPPUNIT_ASSERT_EQUAL((size_t)5, result.size());
        
    // Check all the results
    for(size_t i = 0; i < 5; i++) {
        // Nothing should map, because neither of the bases that could provide
        // credit mapped, because there was disagreement at each.
        CPPUNIT_ASSERT(!result[i].isMapped());
    }
    
}

/**
 * Make sure conflicting placement on a side interferes with mapping from that
 * side.
 */
void CreditFilterTests::testConflictingCredit() {
    
    // Make some left mappings
    std::vector<Mapping> leftMappings {
        Mapping(TextPosition(0, 0)),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(2, 3), 4),
        Mapping(TextPosition(0, 4), 5)
    };
    
    // Make some right mappings
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(1, 34), 5),
        Mapping(TextPosition(3, 33), 4),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(1, 30))
    };
    
    // Make the filter to test.
    CreditFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
    
    // Make sure the result is the right length.
    CPPUNIT_ASSERT_EQUAL((size_t)5, result.size());
        
    CPPUNIT_ASSERT(result[0].isMapped());
    CPPUNIT_ASSERT(result[1].isMapped());
    // Only this middle base should not map.
    CPPUNIT_ASSERT(!result[2].isMapped());
    CPPUNIT_ASSERT(result[3].isMapped());
    CPPUNIT_ASSERT(result[4].isMapped());
    
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 34), result[0].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(3, 33), result[1].getLocation());
    // Middle base didn't map
    CPPUNIT_ASSERT_EQUAL(TextPosition(3, 31), result[3].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 30), result[4].getLocation());
    
}

/**
 * Make sure conflicting placement on one side still allows credit mapping from
 * the other.
 */
void CreditFilterTests::testConflictingCreditOneSideOnly() {
    
    // Make some left mappings
    std::vector<Mapping> leftMappings {
        Mapping(TextPosition(0, 0)),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(0, 4), 5)
    };
    
    // Make some right mappings
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(1, 34), 5),
        Mapping(TextPosition(3, 33), 4),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(1, 30))
    };
    
    // Make the filter to test.
    CreditFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
    
    // Make sure the result is the right length.
    CPPUNIT_ASSERT_EQUAL((size_t)5, result.size());
        
    // Everyone should map.
    CPPUNIT_ASSERT(result[0].isMapped());
    CPPUNIT_ASSERT(result[1].isMapped());
    CPPUNIT_ASSERT(result[2].isMapped());
    CPPUNIT_ASSERT(result[3].isMapped());
    CPPUNIT_ASSERT(result[4].isMapped());
    
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 34), result[0].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(3, 33), result[1].getLocation());
    // These two bases with no mappings should map on credit fom left contexts,
    // even though right contexts disagree about them.
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 32), result[2].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 31), result[3].getLocation());
    
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 30), result[4].getLocation());
    
}

/**
 * Make sure context lengths are honored. Context lengths are inclusive of the
 * base, so 2 is the minimum to map anything on credit.
 */
void CreditFilterTests::testDistance() {
    
    // Make some left mappings
    std::vector<Mapping> leftMappings {
        Mapping(TextPosition(0, 0)),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(0, 4), 2)
    };
    
    // Make some right mappings
    std::vector<Mapping> rightMappings {
        Mapping(TextPosition(1, 34), 2),
        Mapping(),
        Mapping(),
        Mapping(),
        Mapping(TextPosition(1, 30))
    };
    
    // Make the filter to test.
    CreditFilter filter(*index);
    
    // Apply the filter
    std::vector<Mapping> result = filter.apply(leftMappings, rightMappings);
    
    // Make sure the result is the right length.
    CPPUNIT_ASSERT_EQUAL((size_t)5, result.size());
        
    // Check all the results
    
    CPPUNIT_ASSERT(result[0].isMapped());
    CPPUNIT_ASSERT(result[1].isMapped());
    // Only this middle base should not map.
    CPPUNIT_ASSERT(!result[2].isMapped());
    CPPUNIT_ASSERT(result[3].isMapped());
    CPPUNIT_ASSERT(result[4].isMapped());
    
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 34), result[0].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 33), result[1].getLocation());
    // Middle base didn't map
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 31), result[3].getLocation());
    CPPUNIT_ASSERT_EQUAL(TextPosition(1, 30), result[4].getLocation());
    
}









