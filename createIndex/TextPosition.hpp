#ifndef TEXTPOSITION_HPP
#define TEXTPOSITION_HPP

#include <utility>

/**
 * Represents a (text, offset) pair.
 */
class TextPosition {
public:

    /**
     * Create a new TextPosition on the given text at the given 0-based offset.
     */    
    inline TextPosition(size_t text, size_t offset): text(text), 
        offset(offset) {
        
        // Nothing to do!
            
    }
    
    /**
     * Get the text number that this TextPosition is on.
     */
    inline size_t getText() {
        return text;
    }
    
    /**
     * Get the 0-based offset in that text that is being referred to.
     */
    inline size_t getOffset() {
        return offset;
    }
    
protected:
    // The text that the position is on.
    size_t text;
    // The (0-based) offset into that text.
    size_t offset;
} 
 
 
typedef std::pair<size_t, size_t> TextPosition;

#endif
