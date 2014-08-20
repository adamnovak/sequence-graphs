#include "GenericBitVector.hpp"
#include <stdexcept>

#ifdef BITVECTOR_CSA
GenericBitVector::GenericBitVector(std::ifstream& stream): encoder(NULL), 
    bitvector(new BitVector(stream)), size(bitvector->getSize()) {
    
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
    bitvector(new BitVector(stream)), size(bitvector->getSize()) {
    
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
GenericBitVector::GenericBitVector(): encoder(new BitVectorEncoder(32)), 
    bitvector(NULL), size(0) {

    // Nothing to do, already made the encoder.
}
#endif

#ifdef BITVECTOR_SDSL
GenericBitVector::GenericBitVector(): bitvector(), rankSupport(NULL), 
    selectSupport(NULL) {

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
void GenericBitVector::addBit(size_t index) {
    if(encoder == NULL) {
        throw std::runtime_error("Can't add to a vector we didn't create!");
    }
    
    // Pass on valid requests.
    encoder->addBit(index);
}
#endif

#ifdef BITVECTOR_SDSL
void GenericBitVector::addBit(size_t index) {
    if(rankSupport != NULL || selectSupport != NULL) {
        throw std::runtime_error("Can't add to a finished/loaded vector!");
    }

    if(index >= getSize()) {
        // Make sure we have room
        bitvector.resize(index + 1);
    }
    
    // Set the bit
    bitvector[index] = 1;
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
    
    // Set our size
    size = length;    
}
#endif

#ifdef BITVECTOR_SDSL
void GenericBitVector::finish(size_t length) {
    if(rankSupport != NULL || selectSupport != NULL) {
        throw std::runtime_error("Can't finish a finished/loaded vector!");
    }

    if(length > getSize()) {
        // Make sure we are that long
        bitvector.resize(length);
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
size_t GenericBitVector::getSize() const {
    return size;
}
#endif

#ifdef BITVECTOR_SDSL
size_t GenericBitVector::getSize() const {
    return bitvector.size();
}
#endif

#ifdef BITVECTOR_CSA
size_t GenericBitVector::rank(size_t index) const {
    // Take the rank of a position (not in at_least mode)
    return BitVectorIterator(*bitvector).rank(index, false);
}
#endif

#ifdef BITVECTOR_SDSL
size_t GenericBitVector::rank(size_t index) const {
    if(index >= getSize()) {
        // Get the total number of 1s in the bitvector
        return rank(getSize() - 1) + isSet(getSize() - 1);
    }
    
    // Take the rank of a position (not in at_least mode). Correct for
    // disagreement over what rank is.
    return (*rankSupport)(index + 1);
}
#endif

size_t GenericBitVector::rank(size_t index, bool atLeast) const {
    if(atLeast) {
        // Implement this mode ourselves so we can be certain of exactly what it
        // does.
        if(index == 0) {
            // rank(-1) would be 0, but we can't express that with size_ts.
            // So we onlu check if we are in range.
            return (index < getSize());
        }
        // Get the rank of the previous position, then add 1 if this position
        // could contain a 1.
        return rank(index-1) + (index < getSize());
    } else {
        return rank(index);
    }
}

#ifdef BITVECTOR_CSA
bool GenericBitVector::isSet(size_t index) const {
    return BitVectorIterator(*bitvector).isSet(index);
}
#endif

#ifdef BITVECTOR_SDSL
bool GenericBitVector::isSet(size_t index) const {
    return bitvector[index];
}
#endif

#ifdef BITVECTOR_CSA
size_t GenericBitVector::select(size_t one) const {
    // Go select the right position.
    return BitVectorIterator(*bitvector).select(one);
}
#endif

#ifdef BITVECTOR_SDSL
size_t GenericBitVector::select(size_t one) const {
    // Go select the right position. We need to compensate since 0 means the
    // first one to us, but 1 means the first one to SDSL.
    return (*selectSupport)(one + 1);
}
#endif

#ifdef BITVECTOR_CSA
std::pair<size_t, size_t> GenericBitVector::valueAfter(size_t index) const {
    // Get the rank of the next 1 at or after the given position.
    size_t nextRank = rank(index, true) - 1;
    if(index >= getSize()) {
        // We ought to be giving the past-the-end value even though rank dips by
        // 1 at the end of the bitvector.
        nextRank += 1;
    }
    
    // Look up where it is and return that and the rank.
    return std::make_pair(select(nextRank), nextRank);
}
#endif

#ifdef BITVECTOR_SDSL
std::pair<size_t, size_t> GenericBitVector::valueAfter(size_t index) const {
    if(index >= getSize()) {
        // We ought to be giving the past-the-end value
        return std::make_pair(getSize(), rank(getSize() - 1));
    }
    // Assume we're in bounds
    
    // Get the rank of the next 1 at or after the given position.
    size_t nextRank = rank(index, true) - 1;
    
    if(nextRank >= rank(getSize() - 1)) {
        // Too close to the end to do a select. Return past the end value again.
        return std::make_pair(getSize(), rank(getSize() - 1));
    }
    
    // Look up where it is and return that and the rank.
    return std::make_pair(select(nextRank), nextRank);
}
#endif

#ifdef BITVECTOR_CSA
std::pair<size_t, size_t> GenericBitVector::valueBefore(size_t index) const {
    // Get the rank of the 1 before this position, 1-based apparently.
    size_t lastRank = rank(index);
    if(lastRank == 0 || index >= getSize()) {
        // We ought to be giving the past-the-end value.
        lastRank = rank(getSize() + 1);
    } else {
        // Budge left to get the actual 0-based rank.
        lastRank--;
    }
    
    // Look up where it is and return that and the rank.
    return std::make_pair(select(lastRank), lastRank);
}
#endif

#ifdef BITVECTOR_SDSL
std::pair<size_t, size_t> GenericBitVector::valueBefore(size_t index) const {
    // Get the rank of the 1 before this position, 1-based apparently.
    size_t lastRank = rank(index);
    
    if(lastRank == 0) {
        // There is no value before, do the past-the-end result.
        return std::make_pair(getSize(), rank(getSize() - 1));
    }
    
    // Assume we are in bounds.
    
    // Budge left to get the actual 0-based rank.
    lastRank--;
    
    // Look up where it is and return that and the rank.
    return std::make_pair(select(lastRank), lastRank);
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


















