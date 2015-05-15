// Test the credit strategy, which works on graphs.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "CreditStrategyTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( CreditStrategyTests );

// Define constants
const std::string CreditStrategyTests::filename = "Test/duplicated.fa";

CreditStrategyTests::CreditStrategyTests() {
}

CreditStrategyTests::~CreditStrategyTests() {
}

void CreditStrategyTests::setUp() {
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
    
    // Declare every 2 adjacent positions to be a range (since everything
    // appears twice in the input).
    ranges = new GenericBitVector();
    // We need to keep track of the range owners.
    std::map<size_t, TextPosition> rangeOwners;
    for(size_t i = 0; i < index->getBWTLength(); i += 2) {
        // Set every other bit
        ranges->addBit(i);
        
        Log::debug() << index->locate(i) << " should own " << i/2 << std::endl;
        
        if(i % 3 == 0) {
            // And record it owned by the appropriate base. Only do this for
            // some ranges, because it's the exact same thing as the implicit
            // ownership, and we want to test both.
            rangeOwners[i/2] = index->locate(i);
        }
    }
    ranges->finish(index->getBWTLength());

    // Make an FMDIndexView. The CreditStrategy doesn't take ownership of it,
    // since it needs to be able to share it with a MappingScheme.
    view = new FMDIndexView(*index, nullptr, ranges, rangeOwners);

    // Make the credit applier. Leave the mask empty so everything is masked in,
    // but use our ranges. Don't care at all about what positions are assigned
    // to ranges.
    credit = new CreditStrategy(*view);
}


void CreditStrategyTests::tearDown() {
    // Clean up the credit applier
    delete credit;
    
    // And the FMDIndexView
    delete view;
    
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
    
    // Delete the index
    delete index;
    
    // And the ranges bit vector
    delete ranges;
}

/**
 * Make sure credit works.
 */
void CreditStrategyTests::testCredit() {

    // Grab all of the duplicated contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";
    
    // Make a vector of mappings saying everything is unmapped.
    std::vector<Mapping> mappings(query.size(), Mapping());
    
    // Map the first base to contig 0, forward strand, base 0.
    mappings[0] = Mapping(TextPosition(0, 0));
    // map the last base to contig 0, forward strand, base last
    mappings[query.size() - 1] = Mapping(TextPosition(0, query.size() - 1));

    // Apply credit and update the mappings
    credit->applyCredit(query, mappings);

    Log::info() << "Mappings from credit:" << std::endl;
    for(auto mapping : mappings) {
        Log::info() << mapping << std::endl;
    }

    for(size_t i = 0; i < mappings.size(); i++) {
        // Make sure every base is mapped
        CPPUNIT_ASSERT(mappings[i].isMapped());
        
        // Get where each is mapped
        auto mappedTo = mappings[i].getLocation();
        
        // Make sure it is mapped to the right place on text 0 (contig 0 strand
        // 0).
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
    }
}

/**
 * Make sure credit works over mismatches, but doesn't align the mismatches.
 */
void CreditStrategyTests::testCreditOverMismatches() {

    // Grab all of the duplicated contig, with mismatches
    //                      5 v                   28 v
    std::string query = "CATGCCTCGGCGATTCGACGCTCATCTGAGACTCT";
    
    // Make sure mismatches are on.
    credit->maxMismatches = 1;
    
    // Make a vector of mappings saying everything is unmapped.
    std::vector<Mapping> mappings(query.size(), Mapping());
    
    // Map the first base to contig 0, forward strand, base 0.
    mappings[0] = Mapping(TextPosition(0, 0));
    // map the last base to contig 0, forward strand, base last
    mappings[query.size() - 1] = Mapping(TextPosition(0, query.size() - 1));

    // Apply credit and update the mappings
    credit->applyCredit(query, mappings);

    Log::info() << "Mappings from credit:" << std::endl;
    for(auto mapping : mappings) {
        Log::info() << mapping << std::endl;
    }

    for(size_t i = 0; i < mappings.size(); i++) {
        if(i == 5 || i == 28) {
            // These particular bases should be unmapped
            CPPUNIT_ASSERT(!mappings[i].isMapped());
        } else {
    
            // Make sure all other bases are mapped
            CPPUNIT_ASSERT(mappings[i].isMapped());
            
            // Get where each is mapped
            auto mappedTo = mappings[i].getLocation();
            
            // Make sure it is mapped to the right place on text 0 (contig 0 strand
            // 0).
            CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
            CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
    
        }
    
    }
    
}


