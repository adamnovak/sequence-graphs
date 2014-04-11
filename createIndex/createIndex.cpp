#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "FMDIndexBuilder.hpp"

/**
 * createIndex: command-line tool to create a multi-level reference structure.
 */
int main(int argc, char** argv) {
    if(argc < 3) {
        // They forgot their arguments.
        std::cout << "Usage: " << argv[0] << " <index directory> <fasta> "
            << "[<fasta> [<fasta> ...]]" << std::endl;
        std::cout << "The index directory will be deleted and recreated, if it "
            << "exists." << std::endl;
        return 1;
    }
    
    // If we get here, we have the right arguments.
    
    // This holds the directory for the reference structure to build.
    std::string indexDirectory(argv[1]);
    
    // This holds a list of FASTA filenames to load and index.
    std::vector<std::string> fastas;
    for(int i = 2; i < argc; i++) {
        // Put each filename in the vector.
        fastas.push_back(std::string(argv[i]));
    }
    
    // Make sure an empty indexDirectory exists.
    if(boost::filesystem::exists(indexDirectory)) {
        // Get rid of it if it exists already.
        boost::filesystem::remove_all(indexDirectory);
    }
    // Create it again.
    boost::filesystem::create_directory(indexDirectory);
    
    // Build the bottom-level index.
    
    // Make a new builder
    FMDIndexBuilder builder(indexDirectory + "/index.basename");
    for(std::vector<std::string>::iterator i = fastas.begin(); i < fastas.end();
        ++i) {
        
        // Add each FASTA file to the index.
        builder.add(*i);
    }
    
    
    // Load it up with RLCSA.
    
    // Pull contig lengths from the index.
    // Make a ThreadSet with one thread per contig.
    // To construct the non-symmetric merged graph with p context:
        // Traverse the suffix tree down to depth p + 1
        // At each leaf, you have a range:
            // Range = all bases that match this character and share a downstream context of length p.
            // Locate each base in the range.
            // Pinch them all in the correct orientations in the ThreadSet.
                // This does the union-find that we were doing with the complementary connected components problem
        // Then for each block in the ThreadSet:
            // Zip all the segments to get sets of oriented bases.
                // To make the full graph:
                    // Produce two Sides and a Site for each base set.
                    // Connect to upstream/downstream stuff
                        // Accounting for the ends of blocks where you can go many possible places.
                // To allow mapping:
                    // Un-locate everything in the base set.
                    // Make two Sides for each base.
                    // Put them in under every 1-element located base interval in an IntervalTree.
                        // This will re-construct the intervals we pulled out in the suffix-tree traversal.
                    // Save the IntervalTree
                        // Which can make the RangeVector when needed
                        // And stores the bi-directional range<->Position mappings.
    
    return 0;
}
