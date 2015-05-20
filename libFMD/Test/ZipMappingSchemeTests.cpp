// Test the Zip MappingScheme

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "ZipMappingSchemeTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( ZipMappingSchemeTests );

// Define constants
const std::string ZipMappingSchemeTests::filename = "Test/duplicated.fa";

ZipMappingSchemeTests::ZipMappingSchemeTests() {
}

ZipMappingSchemeTests::~ZipMappingSchemeTests() {
}

void ZipMappingSchemeTests::setUp() {
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
    scheme = new ZipMappingScheme<FMDPosition>(FMDIndexView(*index, nullptr,
        ranges));
}


void ZipMappingSchemeTests::tearDown() {
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
void ZipMappingSchemeTests::testMap() {
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
    
    // All the bases should map. This scheme is only weakly stable.
    CPPUNIT_ASSERT_EQUAL(query.size(), mappedBases);
    
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
    
    // All the abses should map. This scheme is only weakly stable.
    CPPUNIT_ASSERT_EQUAL(query2.size(), mappedBases);
}

/**
 * Make sure mapping works with a mask on.
 */
void ZipMappingSchemeTests::testMapWithMask() {

    // Grab all of the duplicated contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";

    // Make a mask for only one of the duplicates
    auto mask = new GenericBitVector();
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        if(index->locate(i).getContigNumber() == 0) {
            // Only include the first contig
            mask->addBit(i);
        }
    }
    mask->finish(index->getBWTLength());

    // Turn the genome restriction on, and don't use a merged ranges vector
    delete scheme;
    scheme = new ZipMappingScheme<FMDPosition>(FMDIndexView(*index, mask));
    
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps in order to the first text.
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    delete mask;
    
    // All of the bases should map
    CPPUNIT_ASSERT_EQUAL(query.size(), mappedBases);    
    
    
}

/**
 * Make sure mapping works with a mask on and ranges.
 */
void ZipMappingSchemeTests::testMapWithMaskAndRanges() {

    // Grab all of the duplicated contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";

    // Make a mask for the other of the duplicates
    auto mask = new GenericBitVector();
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        if(index->locate(i).getContigNumber() == 1) {
            // Only include the second contig
            mask->addBit(i);
        }
    }
    mask->finish(index->getBWTLength());

    // Turn the genome restriction on, and do use the merged ranges vector
    delete scheme;
    scheme = new ZipMappingScheme<FMDPosition>(FMDIndexView(*index, mask,
        ranges));
    
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps in order to text 2 (the one we masked in).
        CPPUNIT_ASSERT_EQUAL((size_t) 2, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    delete mask;
    
    // All of the bases should map
    CPPUNIT_ASSERT_EQUAL(query.size(), mappedBases);    
}

/**
 * Make sure mapping works with FMDPosiutionGroups, but no mismatches.
 */
void ZipMappingSchemeTests::testMapWithGroups() {

    // Grab all of the duplicated contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";

    // Make a mask for the other of the duplicates
    auto mask = new GenericBitVector();
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        if(index->locate(i).getContigNumber() == 1) {
            // Only include the second contig
            mask->addBit(i);
        }
    }
    mask->finish(index->getBWTLength());

    // Turn the genome restriction on, and do use the merged ranges vector
    delete scheme;
    scheme = new ZipMappingScheme<FMDPositionGroup>(FMDIndexView(*index, mask,
        ranges));
    
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps in order to text 2 (the one we masked in).
        CPPUNIT_ASSERT_EQUAL((size_t) 2, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    delete mask;
    
    // All of the bases should map
    CPPUNIT_ASSERT_EQUAL(query.size(), mappedBases);    
}

/**
 * Make sure mapping works with FMDPosiutionGroups and mismatches.
 */
void ZipMappingSchemeTests::testMapWithMismatches() {

    // Grab all of the duplicated contig, but introduce a couple mismatches at
    // 6 and 31.
    //                         v                        v
    std::string query = "CATGCTCCGGCGATTCGACGCTCATCTGCGAATCT";

    // Make a mask for the other of the duplicates
    auto mask = new GenericBitVector();
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        if(index->locate(i).getContigNumber() == 1) {
            // Only include the second contig
            mask->addBit(i);
        }
    }
    mask->finish(index->getBWTLength());

    // Turn the genome restriction on, and do use the merged ranges vector
    delete scheme;
    scheme = new ZipMappingScheme<FMDPositionGroup>(FMDIndexView(*index, mask,
        ranges));
    
    // Allow for mismatches.
    ((ZipMappingScheme<FMDPositionGroup>*) scheme)->mismatchTolerance = 1;
    
    // But require some nearby pins
    ((ZipMappingScheme<FMDPositionGroup>*) scheme)->minUniqueStrings = 2;
    
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure mismatched bases do not map.
        CPPUNIT_ASSERT(i != 6);
        CPPUNIT_ASSERT(i != 31);
        
        // Make sure each other base maps in order to text 2 (the one we masked
        // in).
        CPPUNIT_ASSERT_EQUAL((size_t) 2, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    delete mask;
    
    // All but two of the bases should map
    CPPUNIT_ASSERT_EQUAL(query.size() - 2, mappedBases);    
}

