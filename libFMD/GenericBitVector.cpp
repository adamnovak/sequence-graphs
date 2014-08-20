#include "GenericBitVector.hpp"

GenericBitVector::GenericBitVector(std::ifstream& stream): encoder(NULL), 
    bitvector(new BitVector(stream)), 
    iterator(new BitVectorIterator(*bitvector)), size(bitvector->getSize()) {
    
    // Nothing to do, loaded from the stream.
}

GenericBitVector::GenericBitVector(std::ifstream&& stream): encoder(NULL), 
    bitvector(new BitVector(stream)), 
    iterator(new BitVectorIterator(*bitvector)), size(bitvector->getSize()) {
    
    // Nothing to do, loaded from the stream.
}


GenericBitVector::GenericBitVector(const std::string& filename): 
    GenericBitVector(std::ifstream(filename.c_str(), std::ios::binary)) {
    
    // Nothing to do, delegated to the first constructor.
}

GenericBitVector::GenericBitVector(): encoder(new BitVectorEncoder(32)), 
    bitvector(NULL), iterator(NULL), size(0) {

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
    
    if(iterator != NULL) {
        delete iterator;
    }
}

void GenericBitVector::finish() {
    // Finish encoding
    encoder->flush();
    // Make the BitVector
    bitvector = new BitVector(*encoder, size);    
    // And keep a persistent iterator to look in it.
    iterator = new BitVectorIterator(*bitvector);
}

size_t GenericBitVector::rank(size_t index) {
    // Take the rank of a position (not in at_least mode)
    return iterator->rank(index, false);
}

size_t GenericBitVector::select(size_t one) {
    // Go select the right position
    return iterator->select(one);
}



















