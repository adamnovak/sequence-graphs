// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include <ReadTable.h>
#include <SuffixArray.h>

#include "../FMDIndex.hpp"
#include "../FMDIndexBuilder.hpp"
#include "../util.hpp"

#include "FMDIndexTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( FMDIndexTests );

// Define constants
const std::string FMDIndexTests::filename = "Test/haplotypes.fa";

FMDIndexTests::FMDIndexTests() {

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

FMDIndexTests::~FMDIndexTests() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}

void FMDIndexTests::setUp() {
    
}


void FMDIndexTests::tearDown() {
    
}

/**
 * Test index metadata.
 */
void FMDIndexTests::testMetadata() {
    
    // Make sure it has the right number of characters.
    CPPUNIT_ASSERT(index->getTotalLength() == 35 * 2 * 2);
    
    // Make sure it has the right number of BWT positions (characters + texts).
    CPPUNIT_ASSERT(index->getBWTLength() == index->getTotalLength() + 4);
    
    // Make sure it has the right number of contigs
    CPPUNIT_ASSERT(index->getNumberOfContigs() == 2);
    
    // Make sure it has the right number of genomes.
    CPPUNIT_ASSERT(index->getNumberOfGenomes() == 1);
    
    // Make sure the genome contains the contigs
    CPPUNIT_ASSERT(index->getGenomeContigs(0).first == 0);
    CPPUNIT_ASSERT(index->getGenomeContigs(0).second == 2);
    
    // Make sure the contigs belong to the genome
    CPPUNIT_ASSERT(index->getContigGenome(0) == 0);
    CPPUNIT_ASSERT(index->getContigGenome(1) == 0);
    
    // Make sure all the positions belong to the one genome
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        CPPUNIT_ASSERT_MESSAGE("Genome claims position " + std::to_string(i), 
            index->getGenomeMask(0).isSet(i));
    }
}

/**
 * Test the LF mapping.
 */
void FMDIndexTests::testLF() {

    // Make sure LF behaves sanely and maps the first $ to the first $.
    std::cout << "Index 21 maps LF to " << index->getLF(21) << std::endl;
    for(int i = 0; i < 21; i++) {
        CPPUNIT_ASSERT(index->display(i) != '$');
    }
    CPPUNIT_ASSERT(index->display(21) == '$');
    CPPUNIT_ASSERT(index->displayFirst(0) == '$');
    CPPUNIT_ASSERT(index->getLF(21) == 0);
    
    // Make sure it is consistent and always maps at least to an instance of the
    // correct character.
    for(int i = 0; i < index->getBWTLength(); i++) {
        // Get the first character of where we go, and compare it to the last
        // character of where we are.
        CPPUNIT_ASSERT(index->displayFirst(index->getLF(i)) == 
            index->display(i));
    }

}

/**
 * Test dumping an FMD index's BWT.
 */
void FMDIndexTests::testDump() {
    
    // Dump the entire BWT
    for(int i = 0; i < index->getBWTLength(); i++) {
        // Reconstruct the string.
        std::string reconstruction;
        
        
        // Start at the last character of this row in the BWT.
        int64_t curIndex = i;
        
        do {
            // Put the current character on the reconstruction.
            reconstruction.push_back(index->display(curIndex));
            
            // Go to the previous character in the BWT row.
            curIndex = index->getLF(curIndex);
        } while(curIndex != i);      
    
        // Flip the string around forwards. See
        // <http://www.cplusplus.com/forum/beginner/11633/>
        std::string forwardReconstruction(reconstruction.rbegin(),
            reconstruction.rend());
    
        std::cout << i << ": " << index->displayFirst(i) << " " << 
            forwardReconstruction << " " << index->display(i) << " " << 
            index->isInGenome(i, 0) << " " << index->getLCP(i) << std::endl;
    }
    
}

/**
 * Test pulling out a contig.
 */
void FMDIndexTests::testDisplay() {
    
    // Make sure the first contig comes out right
    CPPUNIT_ASSERT(index->displayContig(0) == 
        "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT");
        
    // And the second
    CPPUNIT_ASSERT(index->displayContig(1) == 
        "CGGGCGCATCGCTATTATTTCTTTCTCTTTTCACA");
        
    std::string contig1("CGGGCGCATCGCTATTATTTCTTTCTCTTTTCACA");
        
    for(size_t i = 0; i < contig1.size(); i++) {
        CPPUNIT_ASSERT_EQUAL(contig1[i], index->display(1, i + 1)); 
    }
}

/**
 * Test searching an an FMD index.
 */
void FMDIndexTests::testSearch() {
    
    // Doesn't find spurious things
    CPPUNIT_ASSERT(index->count("GATTACA").getLength() == 0);
    
    // Finds things which appear once
    CPPUNIT_ASSERT(index->count("TCTTTT").getLength() == 1);
    
    // On both strands
    CPPUNIT_ASSERT(index->count("AAAAGA").getLength() == 1);
    
    // Finds things which appear twice
    CPPUNIT_ASSERT(index->count("TTCG").getLength() == 2);
    
    // Finds whole strands
    CPPUNIT_ASSERT(
        index->count("CGGGCGCATCGCTATTATTTCTTTCTCTTTTCACA").getLength() == 1);
    
}

/**
 * Test locating the positions of things you search up.
 */
void FMDIndexTests::testLocate() {

    // Find and locate a unique thing.
    TextPosition base = index->locate(index->count("TCTTTT").getForwardStart());
    
    // This is on the second sequence, forward strand.
    CPPUNIT_ASSERT(base.getText() == 2);
    // And it has 25 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 25);
    
    // Now let's try all of the first sequence, forwards
    base = index->locate(
        index->count("CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT").getForwardStart());
        
    // This is on the first sequence, forward strand.
    CPPUNIT_ASSERT(base.getText() == 0);
    // And it has 0 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 0);
    
    // And backwards
    base = index->locate(
        index->count("AGAGTCGCAGATGAGCGTCGAATCGCCGAAGCATG").getForwardStart());
        
    // This is on the first sequence, reverse strand.
    CPPUNIT_ASSERT(base.getText() == 1);
    // And it has 0 characters before it.
    CPPUNIT_ASSERT(base.getOffset() == 0);
}

/**
 * Test iterating over the suffix tree.
 */
void FMDIndexTests::testIterate() {

    for(int contextLength = 1; contextLength <= 25; contextLength++) {
        // Try all context lengths shorter than the contigs we put in.
        
        for(FMDIndex::iterator i = index->begin(contextLength); 
            i != index->end(contextLength); ++i) {
            // For each pair of suffix and position in the suffix tree
            
            // Unpack the iterator into pattern and FMDPosition at which it
            // happens.
            std::string pattern = (*i).first;
            // This is in SA coordinates.
            FMDPosition range = (*i).second;
        
            // Make sure we get the same range searching as iterating.
            CPPUNIT_ASSERT(index->count(pattern) == range);
        }
        
    }

}

/**
 * Make sure disambiguating of Mappings works.
 */
void FMDIndexTests::testDisambiguate() {

    // Make some positions    
    Mapping leftMapped(TextPosition(1, 1), 5, 0);
    Mapping rightMapped(TextPosition(1, 1), 0, 4);
    Mapping otherSide(TextPosition(0, 33), 0, 5); 
    Mapping unmapped;
    Mapping elsewhere(TextPosition(2, 10));
    
    // Make sure disambiguate does the right things
    
    CPPUNIT_ASSERT(index != NULL);
    
    // Map things that agree
    CPPUNIT_ASSERT(!index->disambiguate(leftMapped, otherSide).is_mapped);
    CPPUNIT_ASSERT(!index->disambiguate(otherSide, leftMapped).is_mapped);
    
    // Don't map things that disagree, especially if someone put in things with
    // differing orientations
    CPPUNIT_ASSERT_EQUAL(leftMapped, index->disambiguate(leftMapped, unmapped));
    CPPUNIT_ASSERT_EQUAL(rightMapped, index->disambiguate(unmapped, 
        rightMapped));
    
    // Properly integrate contexts when things map.
    CPPUNIT_ASSERT(index->disambiguate(leftMapped, rightMapped).isMapped());
    CPPUNIT_ASSERT_EQUAL(leftMapped.getLocation(), 
        index->disambiguate(leftMapped, rightMapped).getLocation());
    CPPUNIT_ASSERT_EQUAL((size_t)5, index->disambiguate(leftMapped,
        rightMapped).getLeftMaxContext());
    CPPUNIT_ASSERT_EQUAL((size_t)4, index->disambiguate(leftMapped,
        rightMapped).getRightMaxContext());
        
    // Don't map when the two inputs go to completely different places
    CPPUNIT_ASSERT(!index->disambiguate(unmapped, unmapped).is_mapped);
    CPPUNIT_ASSERT(!index->disambiguate(elsewhere, leftMapped).is_mapped);
}

/**
 * Make sure mapping works
 */
void FMDIndexTests::testMap() {
    
    // Grab all of the first contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";
    
    // Map all of the first contig left, right, and both
    std::vector<Mapping> leftMappings = index->mapLeft(query);
    std::vector<Mapping> rightMappings = index->mapRight(query);
    std::vector<Mapping> mappings = index->mapBoth(query);
    
    for(size_t i = 0; i < query.size(); i++ ) {
        Log::output() << leftMappings[i] << "\t" << mappings[i] << "\t" << 
            rightMappings[i] << std::endl;
    }
    
    for(size_t i = 0; i < query.size(); i++ ) {
        // Make sure each base maps in order
        if(leftMappings[i].is_mapped) {
            CPPUNIT_ASSERT(leftMappings[i].location.getText() == 0);
            CPPUNIT_ASSERT(leftMappings[i].location.getOffset() == i);
        }
        
        if(rightMappings[i].is_mapped) {
            CPPUNIT_ASSERT(rightMappings[i].location.getText() == 0);
            CPPUNIT_ASSERT(rightMappings[i].location.getOffset() == i);
        }
        
        CPPUNIT_ASSERT(mappings[i].is_mapped == true);
        CPPUNIT_ASSERT(mappings[i].location.getText() == 0);
        CPPUNIT_ASSERT(mappings[i].location.getOffset() == i);
    }
    
    // Grab all of the first contig's reverse strand.
    std::string query2 = "AGAGTCGCAGATGAGCGTCGAATCGCCGAAGCATG";
    
    // Map it
    mappings = index->mapBoth(query2);
    
    for(size_t i = 0; i < query2.size(); i++ ) {
        // Make sure each base maps in order, but now to the other strand.
        CPPUNIT_ASSERT(mappings[i].is_mapped == true);
        CPPUNIT_ASSERT(mappings[i].location.getText() == 1);
        CPPUNIT_ASSERT(mappings[i].location.getOffset() == i);
    }
    
    // Test again with the genome restriction on
    mappings = index->mapBoth(query, 0);
    
    for(size_t i = 0; i < query2.size(); i++ ) {
        // Make sure each base maps in order still.
        CPPUNIT_ASSERT(mappings[i].is_mapped == true);
        CPPUNIT_ASSERT(mappings[i].location.getText() == 0);
        CPPUNIT_ASSERT(mappings[i].location.getOffset() == i);
    }
    
    // Concatenate them together
    std::string query3 = ("CATGCTTCGGCGATTCGACGCTCATCTGCGACTCTAGAGTCGCAGATGAGCG"
        "TCGAATCGCCGAAGCATG");
    
    // Map all of the concatenated result
    leftMappings = index->mapLeft(query3);
    rightMappings = index->mapRight(query3);
    mappings = index->mapBoth(query3);
    
    Log::output() << "Concatenated:" << std::endl;
    
    for(size_t i = 0; i < query3.size(); i++ ) {
        Log::output() << query3[i] << "\t" << leftMappings[i] << "\t" << 
            mappings[i] << "\t" << rightMappings[i] << std::endl;
    }
    
    for(size_t i = 0; i < query3.size(); i++ ) {
        // Make sure each part of the mapping is right
        if(i < 33) {
            // Should be text 0 (but last 2 are unmapped)
            CPPUNIT_ASSERT(mappings[i].is_mapped == true);
            CPPUNIT_ASSERT(mappings[i].location.getText() == 0);
            CPPUNIT_ASSERT(mappings[i].location.getOffset() == i);
        } else if(i >= 37){
            // Should be text 1 (but first two are unmapped)
            CPPUNIT_ASSERT(mappings[i].is_mapped == true);
            CPPUNIT_ASSERT(mappings[i].location.getText() == 1);
            CPPUNIT_ASSERT(mappings[i].location.getOffset() == i - 35);
        } else {
            // There's a bubble between them
            CPPUNIT_ASSERT(mappings[i].is_mapped == false);
        }
        
    }
}

/**
 * Make sure we can properly map over a mismatch that isn't the first one
 */
void FMDIndexTests::testMapOverMismatch() {

    // Grab all of the first contig, with a mismatch in the middle.
    //                                     v here          v and here
    std::string query = "CATGCTTCGGCGATTCGAAGCTCATCTGCGACTCC";
    
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
        
    // Make sure it maps right before the mismatch, having read through it
    CPPUNIT_ASSERT(mappings[17].getRange() != -1);
    CPPUNIT_ASSERT_EQUAL((size_t)17, mappings[17].getRightMaxContext());
    // Not on the mismatch
    CPPUNIT_ASSERT(mappings[18].getRange() == -1);
    // It can't get further than 1 base max.
    CPPUNIT_ASSERT_EQUAL((size_t)1, mappings[18].getRightMaxContext());
    // But again right after the mismatch
    CPPUNIT_ASSERT(mappings[19].getRange() != -1);
    CPPUNIT_ASSERT_EQUAL((size_t)16, mappings[19].getRightMaxContext());

}

/**
 * Make sure minimum context length is respected.
 */
void FMDIndexTests::testContextLimit() {
    
    // Grab all of the first contig.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";
    
    // Map all of the first contig.
    std::vector<Mapping> mappings = index->mapBoth(query, -1, 100);
    
    for(size_t i = 0; i < query.size(); i++ ) {
        // Make sure no base maps
        CPPUNIT_ASSERT(mappings[i].is_mapped == false);
    }
    
    // Grab all of the first contig's reverse strand.
    std::string query2 = "AGAGTCGCAGATGAGCGTCGAATCGCCGAAGCATG";
    
    // Map it
    mappings = index->mapBoth(query2, -1, 100);
    
    for(size_t i = 0; i < query2.size(); i++ ) {
        // Make sure no base maps this way either
        CPPUNIT_ASSERT(mappings[i].is_mapped == false);
    }
    
    // Test again with the genome restriction on
    mappings = index->mapBoth(query, 0, 100);
    
    for(size_t i = 0; i < query2.size(); i++ ) {
        // Make sure no base maps this way either.
        CPPUNIT_ASSERT(mappings[i].is_mapped == false);
    }
    
    // And one where we actually get mappings.
    mappings = index->mapLeft(query, -1, 10);
    
    for(size_t i = 0; i < query2.size(); i++ ) {
        // Make sure only bases at positions 9+ map. Position 9 is at the end of
        // a length 10 string and so should map.
        if(i < 9) {
            CPPUNIT_ASSERT(mappings[i].is_mapped == false);
        } else {
            CPPUNIT_ASSERT(mappings[i].is_mapped == true);
        }
    }
}

/**
 * Make sure additional context requirement after becoming unique is working.
 */
void FMDIndexTests::testAddContext() {

    // Grab all of the first contig. It becomes unique at "ACTCT" from the right
    // end, since "ACT" is unique, and it then has 2 extra context.
    // "GCGACTCT" has 3 extra context on the right.
    // "CGACTCT" has only 2 extra context on the right.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";
    
    // Declare everything to be a range
    GenericBitVector bv;
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        bv.addBit(i);
    }
    bv.finish(index->getBWTLength());
    
    // Grab (range number, context length) pairs for mapping with 0 additional
    // context and 0 mismatches.
    std::vector<Mapping> mappings = index->misMatchMap(bv,
        query, (GenericBitVector*)NULL, 0, 0, 0.0, 0.0, 0);
        
    std::vector<std::pair<int64_t,size_t>> mappingsExact = index->map(bv, query,
        (GenericBitVector*)NULL, 0, 0);
        
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
        
    // The A in ACTCT should be the last mapped thing
    CPPUNIT_ASSERT(mappings[30].getRange() != -1);
    CPPUNIT_ASSERT(mappings[31].getRange() == -1);
    
    // Try again with 3 context required after uniqueness
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 0, 3, 0.0, 0.0, 0);
    
    mappingsExact = index->map(bv, query, (GenericBitVector*)NULL, 0, 3);
    
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
    
    // The G in GCGACTCT should be the last mapped thing
    CPPUNIT_ASSERT(mappings[27].getRange() != -1);
    CPPUNIT_ASSERT(mappings[28].getRange() == -1);
    
    // And a nonzero min context
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 10, 0, 0.0, 0.0, 0);
    
    mappingsExact = index->map(bv, query, (GenericBitVector*)NULL, 10, 0);
    
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
    
    // The first C in CTGCGACTCT should be the last mapped thing, because it's
    // at the start of a length-10 string.
    CPPUNIT_ASSERT(mappings[25].getRange() != -1);
    CPPUNIT_ASSERT(mappings[26].getRange() == -1);

    // And a nonzero min context with an addContext
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 10, 8, 0.0, 0.0, 0);
    
    mappingsExact = index->map(bv, query, (GenericBitVector*)NULL, 10, 8);
    
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
    
    // ATCT is unique, and then GCGACTCT after it is 8 extra characters, so
    // that first A should be the last thing to map.
    CPPUNIT_ASSERT(mappings[23].getRange() != -1);
    CPPUNIT_ASSERT(mappings[24].getRange() == -1);

}

/**
 * Make sure additional context multiplier is working.
 */
void FMDIndexTests::testMultContext() {

    // Grab all of the first contig. It becomes unique at "ACTCT" from the right
    // end, since "ACT" is unique, and it then has 2 extra context.
    // "GCGACTCT" has 3 extra context on the right.
    // "CGACTCT" has only 2 extra context on the right.
    std::string query = "CATGCTTCGGCGATTCGACGCTCATCTGCGACTCT";
    
    // Declare everything to be a range
    GenericBitVector bv;
    for(size_t i = 0; i < index->getBWTLength(); i++) {
        bv.addBit(i);
    }
    bv.finish(index->getBWTLength());
    
    // Grab (range number, context length) pairs for mapping with 0 additional
    // context and 0 mismatches.
    std::vector<Mapping> mappings = index->misMatchMap(bv,
        query, (GenericBitVector*)NULL, 0, 0, 0.0, 0.0, 0);
        
    std::vector<std::pair<int64_t,size_t>> mappingsExact = index->map(bv, query,
        (GenericBitVector*)NULL, 0, 0, 0.0);
        
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
        
    // The A in ACTCT should be the last mapped thing
    CPPUNIT_ASSERT(mappings[30].getRange() != -1);
    CPPUNIT_ASSERT(mappings[31].getRange() == -1);
    
    // Try again with 2x multiplier.
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 0, 0, 2.0,
        0.0, 0);
    
    mappingsExact = index->map(bv, query, (GenericBitVector*)NULL, 0, 0, 2.0);
    
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
    
    // The C in CTGCGACTCT should be the last mapped thing
    CPPUNIT_ASSERT(mappings[25].getRange() != -1);
    CPPUNIT_ASSERT(mappings[26].getRange() == -1);
    
    // And a nonzero min context
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 10, 0, 0.0, 0.0, 0);
    
    mappingsExact = index->map(bv, query, (GenericBitVector*)NULL, 10, 0);
    
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
    
    // The first C in CTGCGACTCT should be the last mapped thing, because it's
    // at the start of a length-10 string.
    CPPUNIT_ASSERT(mappings[25].getRange() != -1);
    CPPUNIT_ASSERT(mappings[26].getRange() == -1);

    // And a nonzero min context with a multContext
    mappings = index->misMatchMap(bv, query, (GenericBitVector*)NULL, 10, 0, 
        2.0, 0.0, 0);
    
    mappingsExact = index->map(bv, query, (GenericBitVector*)NULL, 10, 0,
        2.0);
    
    for(size_t i = 0; i < mappings.size(); i++) {
        Log::debug() << "Mapping " << i << ": " << mappings[i].getRange() <<
            "," << mappings[i].getRightMaxContext() << " vs " <<
            mappingsExact[i].first  << "," << mappingsExact[i].second <<
            std::endl;
        CPPUNIT_ASSERT(mappings[i].getRange() == mappingsExact[i].first);
        if(mappings[i].isMapped()) {
            CPPUNIT_ASSERT_EQUAL(mappings[i].getRightMaxContext(),
                mappingsExact[i].second);
        }
    }
    
    // CTG is unique, and then CGACTCT after it is 7 extra characters (which
    // more than exceeds 2 * 3 total), and we have 10 total characters, so that
    // first C should be the last thing to map.
    CPPUNIT_ASSERT(mappings[25].getRange() != -1);
    CPPUNIT_ASSERT(mappings[26].getRange() == -1);

}

/**
 * Make sure the LCP array is not lying.
 */
void FMDIndexTests::testLCP() {

    // Save the previous reconstructed string at every step.
    std::string lastReconstruction;

    for(size_t i = 0; i < index->getBWTLength(); i++) {
        // Look at each LCP entry
        size_t lcp = index->getLCP(i);
        
        // Reconstruct the string.
        std::string reconstruction;

        // Start at the last character of this row in the BWT.
        int64_t curIndex = i;
        
        do {
            // Put the current character on the reconstruction.
            reconstruction.push_back(index->display(curIndex));
            
            // Go to the previous character in the BWT row.
            curIndex = index->getLF(curIndex);
        } while(curIndex != i);      
    
        // Flip the string around forwards. See
        // <http://www.cplusplus.com/forum/beginner/11633/>
        std::string forwardReconstruction(reconstruction.rbegin(),
            reconstruction.rend());
            
        
        if(i == 0) {
            // First entry must always be 0
            CPPUNIT_ASSERT(lcp == 0);
        } else {
            // This reconstruction has to mismatch the last one at the
            // appropriate index.
            
            size_t measuredLCP = 0;
            
            while(lastReconstruction[measuredLCP] == 
                forwardReconstruction[measuredLCP] && 
                lastReconstruction[measuredLCP] != '$') {
                
                // While they are characters that can match ("$" never matches
                // another string's "$")...
                measuredLCP++;
            }
            
            // Make sure they got what we got.
            CPPUNIT_ASSERT(measuredLCP == lcp);
            
        }
        
        // Save the reconstruction for the next test
        lastReconstruction = forwardReconstruction;
        
    }
    
}

/**
 * Test left-extending and right-retracting.
 */
void FMDIndexTests::testRetract() {

    // Select the whole thing.
    FMDPosition search = index->getCoveringPosition();
    
    // Search some stuff
    index->extendLeftOnly(search, 'C');
    index->extendLeftOnly(search, 'T');
    index->extendLeftOnly(search, 'T');
    
    // Make an interval to compare against
    FMDPosition countInterval = index->count("TTC");
    
    Log::output() << "For string TTC, got " << search << " by left extend, " <<
        countInterval << " by count." << std::endl;
    
    // Make sure they match in the forward-strand interval, which is the one
    // used here.
    CPPUNIT_ASSERT(search.getForwardStart() == countInterval.getForwardStart());
    CPPUNIT_ASSERT(search.getEndOffset() == countInterval.getEndOffset());
 
    // Reverse back 1 character (to a search string of 2).
    index->retractRightOnly(search, 2);
    
    // Search that instead
    countInterval = index->count("TT");
    
    Log::output() << "For string TT, got " << search << " by retract, " <<
        countInterval << " by count." << std::endl;
    
    // Compare
    CPPUNIT_ASSERT(search.getForwardStart() == countInterval.getForwardStart());
    CPPUNIT_ASSERT(search.getEndOffset() == countInterval.getEndOffset());
    
    // Reverse back to an empty search string.
    index->retractRightOnly(search, 0);
    
    // Search that
    countInterval = index->count("");
    
    Log::output() << "For empty string, got " << search << 
        " by left extend, " << countInterval << " by count." << std::endl;
    
    // Compare
    CPPUNIT_ASSERT(search.getForwardStart() == countInterval.getForwardStart());
    CPPUNIT_ASSERT(search.getEndOffset() == countInterval.getEndOffset());
    
}












