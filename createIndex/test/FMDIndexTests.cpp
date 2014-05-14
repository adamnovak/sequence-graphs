// Test the BWT generation.

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

#include "Util/ReadTable.h"
#include "SuffixArray.h"

#include "FMDIndex.hpp"
#include "FMDIndexBuilder.hpp"
#include "util.hpp"

#include "FMDIndexTests.hpp"


// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( FMDIndexTests );

// Define constants
const std::string FMDIndexTests::filename = "test/haplotypes.fa";

void FMDIndexTests::setUp() {
    // We ened a built index as a fixture.

    // Set up a temporary directory to put the index in.
    tempDir = make_tempdir();
    
    // Make a new index builder.
    FMDIndexBuilder builder(tempDir + "/index.basename");
    
    // Add the haplotypes file
    builder.add(filename);
    
    // Finish the index.
    builder.close();
    
}


void FMDIndexTests::tearDown() {
    // Get rid of the temporary index directory
    boost::filesystem::remove_all(tempDir);
}


/**
 * Test dumping an FMD index's BWT.
 */
void FMDIndexTests::testDump() {
    
    // Load the index up
    FMDIndex index(tempDir + "/index.basename");
    
    // Make sure it has the right number of characters.
    CPPUNIT_ASSERT(index.getTotalLength() == 35 * 2 * 2);
    
    // Make sure it has the right number of BWT positions (characters + texts).
    CPPUNIT_ASSERT(index.getBWTLength() == index.getTotalLength() + 4);
    
    // make sure LF behaves sanely and maps the first $ to the first $.
    std::cout << "Index 21 maps LF to " << index.getLF(21) << std::endl;
    for(int i = 0; i < 21; i++) {
        CPPUNIT_ASSERT(index.display(i) != '$');
    }
    CPPUNIT_ASSERT(index.display(21) == '$');
    CPPUNIT_ASSERT(index.displayFirst(0) == '$');
    CPPUNIT_ASSERT(index.getLF(21) == 0);
    
    // Dump it
    for(int i = 0; i < index.getBWTLength(); i++) {
        // Reconstruct the string.
        std::string reconstruction;
        
        
        // Start at the last character of this row in the BWT.
        int64_t curIndex = i;
        
        do {
            // Put the current character on the reconstruction.
            reconstruction.push_back(index.display(curIndex));
            
            // Go to the previous character in the BWT row.
            curIndex = index.getLF(curIndex);
        } while(curIndex != i);      
    
        // Flip the string around forwards. See
        // <http://www.cplusplus.com/forum/beginner/11633/>
        std::string forwardReconstruction(reconstruction.rbegin(),
            reconstruction.rend());
    
        std::cout << i << ": " << index.displayFirst(i) << " " << 
            forwardReconstruction << " " << index.display(i) << std::endl;
    }
    
}











