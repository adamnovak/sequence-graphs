// Test the Natural MappingScheme

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "NaturalMappingSchemeTests2.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( NaturalMappingSchemeTests2 );

// Define constants
const std::string NaturalMappingSchemeTests2::filename = "Test/duplicated.fa";

NaturalMappingSchemeTests2::NaturalMappingSchemeTests2() {
}

NaturalMappingSchemeTests2::~NaturalMappingSchemeTests2() {
}

void NaturalMappingSchemeTests2::setUp() {
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
    for(size_t i = 0; i < index->getBWTLength(); i += 2) {
        ranges->addBit(i);
    }
    ranges->finish(index->getBWTLength());

    // Make the mapping scheme. Leave the mask empty so everything is masked in,
    // but use our ranges. Don't care at all about what positions are assigned
    // to ranges.
    scheme = new NaturalMappingScheme(*index, nullptr, ranges);
}


void NaturalMappingSchemeTests2::tearDown() {
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
 * Make sure mapping works
 */
void NaturalMappingSchemeTests2::testMap() {
    // Grab all of the duplicated contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";
    
    size_t mappedBases = 0;
    
    // Map everything and get callbacks
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Everything that maps should be mapped to the correct position on text
        // 0 or 2, which are the ones we should have merged in this direction.
        CPPUNIT_ASSERT(mappedTo.getText() == 0 || mappedTo.getText() == 2);
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    // All the bases inside the protection of the outermost min matchings should
    // map. The first mappable position is 4, and the last is 31, giving 28
    // mappable positions.
    CPPUNIT_ASSERT_EQUAL(query.size() - 4 - 3, mappedBases);
    
    // Grab all of the first contig's reverse strand.
    std::string query2 = "AGAGTCGCAGATGAGCGTCGAATCGCCGAAGCATG";
    
    // Reset for the next mapping
    mappedBases = 0;
    
    // Map everything and get callbacks again
    scheme->map(query2, [&](size_t i, TextPosition mappedTo) {
        // Everything that maps should be mapped to the correct position on text
        // 1 or 3, which are the ones we should have merged in this direction.
        CPPUNIT_ASSERT(mappedTo.getText() == 1 || mappedTo.getText() == 3);
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    // All mappable bases should map; the uniqueness structure is the same as
    // the forward strand.
    CPPUNIT_ASSERT_EQUAL(query2.size() - 4 - 3, mappedBases);
}

