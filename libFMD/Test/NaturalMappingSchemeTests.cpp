// Test the Natural MappingScheme

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../MismatchResultSet.hpp"
#include "../util.hpp"

#include "NaturalMappingSchemeTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( NaturalMappingSchemeTests );

// Define constants
const std::string NaturalMappingSchemeTests::filename = "Test/haplotypes.fa";

NaturalMappingSchemeTests::NaturalMappingSchemeTests() {
}

NaturalMappingSchemeTests::~NaturalMappingSchemeTests() {
}

void NaturalMappingSchemeTests::setUp() {
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
    scheme = new NaturalMappingScheme(*index, *ranges);
}


void NaturalMappingSchemeTests::tearDown() {
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
void NaturalMappingSchemeTests::testMap() {
    // Grab all of the first contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";
    
    size_t mappedBases = 0;
    
    // Map everything and get callbacks
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Everything that maps should be mapped to the correct position on text
        // 0, since we threw in all of text 0.
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    // All the bases inside the protection of the outermost min matchings should
    // map. The first mappable position is 4, and the last is 30, giving 27
    // mappable positions.
    CPPUNIT_ASSERT_EQUAL(query.size() - 4 - 4, mappedBases);
    
    // Grab all of the first contig's reverse strand.
    std::string query2 = "AGAGTCGCAGATGAGCGTCGAATCGCCGAAGCATG";
    
    // Reset for the next mapping
    mappedBases = 0;
    
    // Map everything and get callbacks again
    scheme->map(query2, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps in order, but now to the other strand.
        CPPUNIT_ASSERT_EQUAL((size_t) 1, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    // All mappable bases should map; the uniqueness structure is the same as
    // the forward strand.
    CPPUNIT_ASSERT_EQUAL(query2.size() - 4 - 4, mappedBases);
    
    // Concatenate them together
    std::string query3 = ("CATGCTTCGGCGATTCGACGCTCATCTGCGACTCTAGAGTCGCAGATGAGCG"
        "TCGAATCGCCGAAGCATG");
    
    mappedBases = 0;
    scheme->map(query3, [&](size_t i, TextPosition mappedTo) {
        // Make sure each part of the mapping is right
        if(i < 33) {
            // Should be text 0 (but last 2 are unmapped due to left vs right
            // mapping conflict)
            CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
            CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        } else if(i >= 37){
            // Should be text 1 (but first two are unmapped due to left vs right
            // mapping conflict)
            CPPUNIT_ASSERT_EQUAL((size_t) 1, mappedTo.getText());
            CPPUNIT_ASSERT_EQUAL(i - 35, mappedTo.getOffset());
        } else {
            // There's a bubble between them. The callback should never be
            // called for these bases.
            CPPUNIT_ASSERT(false);
        }
        mappedBases++;
    });
    
    // All the bases not in that gap and not too close to the ends should have
    // reported in.
    CPPUNIT_ASSERT_EQUAL(query3.size() - (37 - 33) - 4 - 4, mappedBases);
}

void NaturalMappingSchemeTests::testMapWithMask() {

    // Grab all of the first contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";

    // Turn the genome restriction on
    delete scheme;
    scheme = new NaturalMappingScheme(*index, *ranges,
        &index->getGenomeMask(0));
    
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps in order still.
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    // Same number of bases as in previous test should map since there is just
    // the one genome.
    CPPUNIT_ASSERT_EQUAL(query.size() - 4 - 4, mappedBases);    
}











