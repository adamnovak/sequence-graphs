// Test the changes I made to quicksort.

#include <vector>

#include "../ReadTable.h"
#include "../SuffixArray.h"
#include "../SuffixCompare.h"
#include "../mkqs.h"

#include "mkqsTests.h"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( MKQSTests );

// Define constants
const std::string MKQSTests::filename = "Test/haplotypes.fa";

void MKQSTests::setUp() {

    std::cout << "Loading MKQS test data..." << std::endl;

    // Set up the basic suffix array
    readTable = new ReadTable(filename);
    suffixArray = new SuffixArray(readTable, 1);

}


void MKQSTests::tearDown() {
    delete suffixArray;
    delete readTable;
}

/**
 * Test building a sampled suffix array.
 */
void MKQSTests::testSort() {
    
    // Make the comparison objects to use in the sort.
    SuffixCompareRadix radix_compare(readTable, 6);
    SuffixCompareIndex index_compare;
    
    // Get our own copy of the SA data since it's private
    SAElemVector elems;
    for(size_t i = 0; i < suffixArray->getSize(); i++) {
        elems.push_back(suffixArray->get(i));
    }
    
    // Sort the whole suffix array
    std::cout << "Sorting " << elems.size() << " elements" << std::endl;
    // Do it in parallel.
    parallel_mkqs(&elems[0], elems.size(), 100, radix_compare,
        index_compare);
    std::cout << "Checking order" << std::endl;
    
    for(int i = 0; i < (int)elems.size() - 1; i++) {
        // Scan it and assert order.
        
        SAElem a = elems[i];
        SAElem b = elems[i + 1];
        
        // We need a < b

        // First compare the stings, and ensure a's string <= b's.
        int comparison = strcmp(radix_compare.getChrPtr(a),
            radix_compare.getChrPtr(b));
        
        CPPUNIT_ASSERT(comparison <= 0);
        
        if(comparison == 0) {
        
            // If they're the same, index_compare must sort them right.
            CPPUNIT_ASSERT(index_compare(a, b));
        
        }
        
        
        
    }

}
