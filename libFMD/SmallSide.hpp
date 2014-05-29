#ifndef SMALLSIDE_HPP
#define SMALLSIDE_HPP

#include <fstream>

/**
 * Represents a Side (i.e. a base ID and a left or right Face) in a succinct
 * way. Borrows one bit off the ID field to store the face/orientation, and thus
 * manages to cram the whole thing into 8 bytes.
 */
class SmallSide {

public:
    /**
     * Make a new SmallSide for the given coordinate (63 bits) and face (1 bit).
     */
    inline SmallSide(size_t coordinate, bool face): 
        bits((coordinate << 1) | face) {
        
        // Nothing to do. Already packed coordinate and face into the same
        // field.
        
    }
    
    /**
     * Get the coordinate from a SmallSide.
     */
    inline size_t getCoordinate() {
        // This is the high 63 bits, without sign extension.
        return bits >> 1;
    }
    
    /**
     * Get the face (0 for left, 1 for right) from a SmallSide.
     */
    inline bool getFace() {
        // This is the low 1 bit.
        return bits & 0x1;
    }
    
    /**
     * Write this SmallSide to a stream, in a binary format.
     */
    inline void write(std::ostream& stream) {
        // Write our bits to the stream.
        stream.write((char*)&bits, sizeof(bits));
    }
    
private:
    // Pack everything into 8 bytes. Assuming a size_t really is 8 bytes.
    size_t bits;    

};

#endif
