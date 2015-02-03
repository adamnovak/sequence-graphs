#include "GenericBitVector.hpp"
#include <stdexcept>
#include <thread>

#ifdef BITVECTOR_CSA
GenericBitVector::GenericBitVector(std::ifstream& stream): encoder(NULL), 
    bitvector(new BitVector(stream)), size(bitvector->getSize()), 
    iterator(new BitVectorIterator(*bitvector)), iteratorMutex() {
    
    // Nothing to do, loaded from the stream.
}
#endif

#ifdef BITVECTOR_SDSL
GenericBitVector::GenericBitVector(std::ifstream& stream): bitvector(),
    rankSupport(NULL), selectSupport(NULL) {
    
    // Deserialize the bitvector
    sdsl::load(bitvector, stream);
    
    // Make the supports for rank and select
    rankSupport = new sdsl::rank_support_v5<>(&(bitvector));
    selectSupport = new sdsl::select_support_mcl<>(&(bitvector));
    
    
    
}
#endif

#ifdef BITVECTOR_CSA
GenericBitVector::GenericBitVector(std::ifstream&& stream): encoder(NULL), 
    bitvector(new BitVector(stream)), size(bitvector->getSize()), 
    iterator(new BitVectorIterator(*bitvector)), iteratorMutex() {
    
    // Nothing to do, loaded from the stream.
}
#endif

#ifdef BITVECTOR_SDSL
GenericBitVector::GenericBitVector(std::ifstream&& stream): bitvector(),
    rankSupport(NULL), selectSupport(NULL) {
    
    // Deserialize the bitvector
    sdsl::load(bitvector, stream);
    
    // Make the supports for rank and select
    rankSupport = new sdsl::rank_support_v5<>(&(bitvector));
    selectSupport = new sdsl::select_support_mcl<>(&(bitvector));
    
}
#endif


GenericBitVector::GenericBitVector(const std::string& filename): 
    GenericBitVector(std::ifstream(filename.c_str(), std::ios::binary)) {
    
    // Nothing to do, delegated to the first constructor.
}

#ifdef BITVECTOR_CSA
GenericBitVector::GenericBitVector(size_t sizeHint): encoder(new BitVectorEncoder(32)), 
    bitvector(NULL), size(0), iterator(NULL), iteratorMutex() {

    // Nothing to do, already made the encoder.
    // TODO: actually use the size hint.
}
#endif

#ifdef BITVECTOR_SDSL
GenericBitVector::GenericBitVector(size_t sizeHint): bitvector(size_hint),
    rankSupport(NULL), selectSupport(NULL) {

    // Nothing to do
}
#endif

#ifdef BITVECTOR_CSA
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
#endif

#ifdef BITVECTOR_SDSL
GenericBitVector::~GenericBitVector() {
    // Clean up any supports we made.
    
    if(rankSupport != NULL) {
        delete rankSupport;
    }
    
    if(selectSupport != NULL) {
        delete selectSupport;
    }
    
}
#endif

#ifdef BITVECTOR_CSA
void GenericBitVector::finish(size_t length) {
    if(encoder == NULL) {
        throw std::runtime_error("Can't finish a vector we didn't create!");
    }
    
    // Finish encoding
    encoder->flush();
    // Make the BitVector
    bitvector = new BitVector(*encoder, length);
    // And an iterator for it
    iterator = new BitVectorIterator(*bitvector);
    
    // Set our size
    size = length;    
}
#endif

#ifdef BITVECTOR_SDSL
void GenericBitVector::finish(size_t length) {
    if(rankSupport != NULL || selectSupport != NULL) {
        throw std::runtime_error("Can't finish a finished/loaded vector!");
    }

    // What's the size of the bitvector before we make it bigger?
    size_t oldSize = getSize();

    if(length > oldSize) {
        // Make sure we have room
        bitvector.resize(length);
        
        for(size_t i = oldSize; i < length; i++) {
            // Make sure all the bits we just added in are 0s.
            bitvector[i] = 0;
        }
        
    }
    
    // Make the supports for rank and select
    rankSupport = new sdsl::rank_support_v5<>(&(bitvector));
    selectSupport = new sdsl::select_support_mcl<>(&(bitvector));
}
#endif

#ifdef BITVECTOR_CSA
void GenericBitVector::writeTo(std::ofstream& stream) const {
    bitvector->writeTo(stream);
}
#endif

#ifdef BITVECTOR_SDSL
void GenericBitVector::writeTo(std::ofstream& stream) const {
    sdsl::serialize(bitvector, stream);
    
    // TODO: save supports.
}
#endif

#ifdef BITVECTOR_CSA
GenericBitVector* GenericBitVector::createUnion(
    const GenericBitVector& other) const {

    // We need to hack a BitVector into a GenericBitVector.
    
    // Make one to populate.
    GenericBitVector* toReturn = new GenericBitVector();
    
    // make the bitvector to put in
    BitVector* unionBitvector = bitvector->createUnion(*(other.bitvector));
    
    // Put it in
    toReturn->bitvector = unionBitvector;
    // Populate the size
    toReturn->size = unionBitvector->getSize();
    // And fill in the iterator
    toReturn->iterator = new BitVectorIterator(*unionBitvector);
    
    // Return the new GenericBitVector holding the BitVector we made.
    return toReturn;
    
}
#endif

#ifdef BITVECTOR_SDSL
GenericBitVector* GenericBitVector::createUnion(
    const GenericBitVector& other) const {

    if(getSize() > other.getSize()) {
        // Make sure we are the shorter bitvector.
        return other.createUnion(*this);
    }

    // We need to hack a BitVector into a GenericBitVector.
    
    // Make one to populate.
    GenericBitVector* toReturn = new GenericBitVector();
    
    for(size_t i = 0; i < getSize(); i++) {
        // OR all the bits we both have
        if(isSet(i) || other.isSet(i)) {
            toReturn->addBit(i);
        }
    }
    
    for(size_t i = getSize(); i < other.getSize(); i++) {
        // And add in extra bits that only the other has.
        if(other.isSet(i)) {
            toReturn->addBit(i);
        }
    }

    // Finish it to the longer length.
    toReturn->finish(other.getSize());
    
    // Return the new GenericBitVector we made.
    return toReturn;
    
}
#endif


















