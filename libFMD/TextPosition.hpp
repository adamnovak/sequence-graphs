#ifndef TEXTPOSITION_HPP
#define TEXTPOSITION_HPP

#include <utility>

/**
 * Represents a (text, offset) pair. Every contig ends up as two sequential
 * texts, one for the forwards strand and one for the reverse.
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
    inline size_t getText() const {
        return text;
    }
    
    /**
     * Set the text number to the given value.
     */
    inline void setText(size_t value) {
        text = value;
    }
    
    /**
     * Get the 0-based offset in that text that is being referred to.
     */
    inline size_t getOffset() const {
        return offset;
    }
    
    /**
     * Set the 0-based offset in the text to the given value.
     */
    inline void setOffset(size_t value) {
        offset = value;        
    }
    
    /**
     * Provide equality comparison for testing.
     */
    inline bool operator==(const TextPosition& other) const {
        // Just compare our texts and offsets.
        return getText() == other.getText() && getOffset() == other.getOffset();
    }
    
    /**
     * Provide inequality comparison for testing.
     */
    inline bool operator!=(const TextPosition& other) const {
        return !(*this == other);
    }
    
protected:
    // The text that the position is on.
    size_t text;
    // The (0-based) offset into that text.
    size_t offset;
}; 

#endif
