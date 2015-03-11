// Test the Natural MappingScheme

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
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
    
    // Make the mapping scheme, using a view of the whole index with no merging.
    scheme = new NaturalMappingScheme(FMDIndexView(*index));
}


void NaturalMappingSchemeTests::tearDown() {
    // Clean up the mapping scheme
    delete scheme;
    
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
    
    // Delete the index
    delete index;
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

/**
 * Make sure mapping works with a (trivial) mask on.
 */
void NaturalMappingSchemeTests::testMapWithMask() {

    // Grab all of the first contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";

    // Turn the genome restriction on
    delete scheme;
    scheme = new NaturalMappingScheme(FMDIndexView(*index,
        &index->getGenomeMask(0)));
    
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

/**
 * Make sure synteny chaining can cross a mismatch.
 */
void NaturalMappingSchemeTests::testSkipMismatches() {
    // Grab all of the first contig, with an error in the middle.
    //                                     v Here
    std::string query = "CATGCTTCGGCGATTCGAGGCTCATCTGCGACTCT";
    
    // Configure the scheme to be slightly more tolerant, but also so that it
    // needs to chain max matchings..
    scheme->minHammingBound = 6;
    scheme->maxHammingDistance = 1;
    
    // Map the query and see how much maps.
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps in order still.
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        CPPUNIT_ASSERT_EQUAL(i, mappedTo.getOffset());
        mappedBases++;
    });
    
    // We should not map that changed base, but the rest should map as before.
    CPPUNIT_ASSERT_EQUAL(query.size() - 4 - 4 - 1, mappedBases);   
}

/**
 * Make sure synteny chaining can cross an insertion.
 */
void NaturalMappingSchemeTests::testSkipInserts() {
    // Grab all of the first contig, with an insertion in the middle.
    //                               v Here      
    std::string query = "CATGCTTCGGCGTATTCGACGCTCATCTGCGACTCT";
    
    // Configure the scheme to be slightly more tolerant, but also so that it
    // needs to chain max matchings.
    scheme->minHammingBound = 6;
    scheme->maxHammingDistance = 1;
    
    // Map the query and see how much maps.
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps on the right text
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        // TODO: check the offsets.
        mappedBases++;
        
        Log::info() << i << ": " << query[i] << " mapped to " << mappedTo <<
        std::endl;
        
    });
    
    // We should not map that inserted base, but the rest should map as before.
    CPPUNIT_ASSERT_EQUAL(query.size() - 4 - 4 - 1, mappedBases);   
}

/**
 * Make sure synteny chaining can cross a deletion.
 */
void NaturalMappingSchemeTests::testSkipDeletes() {
    // Grab all of the first contig, with a deletion in the middle.
    //                                 v Was here  
    std::string query = "CATGCTTCGGCGATCGACGCTCATCTGCGACTCT";
    
    // Configure the scheme to be slightly more tolerant, but also so that it
    // needs to chain max matchings.
    scheme->minHammingBound = 6;
    scheme->maxHammingDistance = 1;
    
    // Map the query and see how much maps.
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps on the right text
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        // TODO: check the offsets.
        mappedBases++;
    });
    
    // We should not map the one remaining T of the double T, but the rest
    // should map as before.
    CPPUNIT_ASSERT_EQUAL(query.size() - 4 - 4 - 1, mappedBases);   
}

/**
 * Make sure synteny chaining can cross several errors.
 */
void NaturalMappingSchemeTests::testSkipAll() {
    // Grab all of the first contig, with an insertion, mismatch, and deletion.
    // Chosen so as not to create spurious matchings and thus blacklistings I
    // don't expect.
    //                                v     v Here    v
    std::string query = "CATGCTTCGGCGAGTTCGAGGCTCATCTGGACTCT";
    
    // Configure the scheme so that we can't map without chaining max matchings
    scheme->minHammingBound = 6;
    
    // Map the query and see how much maps, without mismatch tolerance
    size_t mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Count bases that map
        mappedBases++;
    });
    
    // Nothing should map.
    CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedBases);
    
    // Configure the scheme to be quite tolerant, so we can actually do the
    // chaining.
    scheme->maxHammingDistance = 4;
    
    // Map the query and see how much maps.
    mappedBases = 0;
    scheme->map(query, [&](size_t i, TextPosition mappedTo) {
        // Make sure each base maps to nthe right text.
        // TODO: also make sure theyt map to the right places.
        CPPUNIT_ASSERT_EQUAL((size_t) 0, mappedTo.getText());
        mappedBases++;
        
        Log::info() << i << " (" << query[i] << "): " << mappedTo << std::endl;
    });
    
    // We should not map the changed or inserted bases, or very near the ends,
    // but the rest should map
    CPPUNIT_ASSERT_EQUAL(query.size() - 4 - 4 - 2, mappedBases);   
}









