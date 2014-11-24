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
     * Return which strand our text is: 0 (forward) or 1 (reverse)
     */
    inline bool getStrand() const {
        return getText() % 2 == 1;
    }
    
    /**
     * Return which contig our text is.
     */
    inline size_t getContigNumber() const {
        return getText() / 2;
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
     * Shift forwards along the forward strand, or backwards along the reverse
     * strand, by this amount.
     */
    inline void addOffset(int64_t offset) {
        // Offset in the correct direction depending on its strand. Go backwards
        // on reverse strands and forwards on forwards strands.
        setOffset(getOffset() + (getStrand() ? -1 : 1) * offset);
    }
    
    /**
     * Flip to the other text that belongs to the same contig, converting offset
     * with the given contig length.
     */
    inline void flip(size_t textLength) {
        // Flip the LSB of the text
        text ^= 1;
        // Compute the offset from the other end. -1 since end is exclusive.
        offset = (int64_t) textLength - (int64_t) offset - 1;
    }
    
    /**
     * Returns true if the other position is exactly offset after this position
     * on the same strand, and false otherwise.
     */
    inline bool isConsistent(TextPosition other, int64_t offset) {
        // They have to be on the same text, with the same offset between them.
        return getText() == other.getText() && 
            other.getOffset() - getOffset() == offset;
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
    
    /**
     * Provide less-than comparison for sets.
     */
    inline bool operator<(const TextPosition& other) const {
        // Just compare our texts and offsets, with text ranked as more
        // important.
        return getText() < other.getText() || (getText() == other.getText() &&
            getOffset() < other.getOffset());
    }
    
    /**
     * Provide greater-than comparison for sets.
     */
    inline bool operator>(const TextPosition& other) const {
        // Just compare our texts and offsets, with text ranked as more
        // important.
        return getText() > other.getText() || (getText() == other.getText() &&
            getOffset() > other.getOffset());
    }
    
protected:
    // The text that the position is on.
    size_t text;
    // The (0-based) offset into that text.
    size_t offset;
}; 

#endif
