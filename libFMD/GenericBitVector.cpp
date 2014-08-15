#include "GenericBitVector.hpp"
#include <stdexcept>

GenericBitVector::GenericBitVector(std::ifstream& stream): encoder(NULL), 
    bitvector(new BitVector(stream)), size(bitvector->getSize()) {
    
    // Nothing to do, loaded from the stream.
}

GenericBitVector::GenericBitVector(std::ifstream&& stream): encoder(NULL), 
    bitvector(new BitVector(stream)), size(bitvector->getSize()) {
    
    // Nothing to do, loaded from the stream.
}


GenericBitVector::GenericBitVector(const std::string& filename): 
    GenericBitVector(std::ifstream(filename.c_str(), std::ios::binary)) {
    
    // Nothing to do, delegated to the first constructor.
}

GenericBitVector::GenericBitVector(): encoder(new BitVectorEncoder(32)), 
    bitvector(NULL), size(0) {

    // Nothing to do, already made the encoder.
}

GenericBitVector::~GenericBitVector() {
    // We own all our fields, and we can't be assigned or anything silly. Clean
    // them up.
    
    if(encoder != NULL) {
        delete encoder;
    }
    
    if(bitvector != NULL) {
        delete bitvector;
    }
}

void GenericBitVector::addBit(size_t index) {
    if(encoder == NULL) {
        throw std::runtime_error("Can't add to a vector we didn't create!");
    }
    
    // Pass on valid requests.
    encoder->addBit(index);
}

void GenericBitVector::finish(size_t length) {
    if(encoder == NULL) {
        throw std::runtime_error("Can't finish a vector we didn't create!");
    }
    
    // Finish encoding
    encoder->flush();
    // Make the BitVector
    bitvector = new BitVector(*encoder, length);
    
    // Set our size
    size = length;    
}

void GenericBitVector::writeTo(std::ofstream& stream) const {
    bitvector->writeTo(stream);
}

size_t GenericBitVector::getSize() const {
    return size;
}

size_t GenericBitVector::rank(size_t index) const {
    // Take the rank of a position (not in at_least mode)
    return BitVectorIterator(*bitvector).rank(index, false);
}

size_t GenericBitVector::rank(size_t index, bool atLeast) const {
    if(atLeast) {
        // Implement this mode ourselves so we can be certain of exactly what it
        // does.
        if(index == 0) {
            // rank(-1) would be 0, but we can't express that with size_ts.
            // So we onlu check if we are in range.
            return (index < size);
        }
        // Get the rank of the previous position, then add 1 if this position
        // could contain a 1.
        return rank(index-1) + (index < size);
    } else {
        return rank(index);
    }
}

bool GenericBitVector::isSet(size_t index) const {
    return BitVectorIterator(*bitvector).isSet(index);
}

size_t GenericBitVector::select(size_t one) const {
    // Go select the right position
    return BitVectorIterator(*bitvector).select(one);
}

std::pair<size_t, size_t> GenericBitVector::valueAfter(size_t index) const {
    // Get the rank of the next 1 at or after the given position.
    size_t nextRank = rank(index, true) + 1;
    return std::make_pair(select(nextRank), nextRank);
}

std::pair<size_t, size_t> GenericBitVector::valueBefore(size_t index) const {
    // Get the rank of the 1 before this position.
    size_t lastRank = rank(index);
    return std::make_pair(select(lastRank), lastRank);
}

GenericBitVector* GenericBitVector::createUnion(GenericBitVector& other) {
    // We need to hack a BitVector into a GenericBitVector.
    
    // Make one to populate.
    GenericBitVector* toReturn = new GenericBitVector();
    
    // make the bitvector to put in
    BitVector* unionBitvector = bitvector->createUnion(*(other.bitvector));
    
    // Put it in
    toReturn->bitvector = unionBitvector;
    // Populate the size
    toReturn->size = unionBitvector->getSize();
    
    // Return the new GenericBitVector holding the BitVector we made.
    return toReturn;
    
}



















