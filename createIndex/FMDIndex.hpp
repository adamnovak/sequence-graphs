#include <string>
#include <vector>
#include <rlcsa/fmd.h>


/**
 * A class that encapsulates access to an FMDIndex, consisting of the underlying
 * FMD and some auxilliary structures (cointig names and lengths) that we store
 * along with it.
 */
class FMDIndex {
public:
    /**
     * Load an FMD and metadata from the given basename.
     */
    FMDIndex(std::string basename);
    
    /**
     * Get the contig number (not the text number) from a (text, offset) pair.
     */
    CSA::usint getContigNumber(CSA::pair_type base) const;
    
    /**
     * Get the strand from a (text, offset) pair: either 0 or 1.
     */
    bool getStrand(CSA::pair_type base) const;
    
    /**
     * Given a pair_type representing a base, and a vector of contig lengths by
     * number, determine the specified base's offset from the left of its
     * contig, 1-based.
     */
    CSA::usint getOffset(CSA::pair_type base) const;
    
    /**
     * Get a unique string name for a position.
     */
    std::string getName(CSA::pair_type base) const;
    
    /**
     * Get the total length of all contigs, on both strands.
     */
    CSA::usint getTotalLength() const;
    
    /**
     * Get the character for a base at the given offset on the given contig,
     * oriented as on the given strand. Offset is 1-based from contig start.
     */
    char display(CSA::usint contig, CSA::usint offset, bool strand) const;
    
    /**
     * Get the character for a base.
     */
    char display(CSA::pair_type base) const; 
    
    /**
     * Holds the actual underlying FMD.
     */
    CSA::FMD fmd;
    
    /**
     * Holds the names of all the contigs.
     */
    std::vector<std::string> names;
    
    /**
     * Holds the lengths of all the contigs, in the same order.
     */
    std::vector<long long int> lengths;

};
