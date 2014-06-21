#ifndef CSA_BITVECTOR_HPP
#define CSA_BITVECTOR_HPP

#include <fstream>

#include "BitVectorBase.hpp"

namespace CSA {

/*
  This class is used to construct a BitVector.
*/

class BitVectorEncoder : public VectorEncoder
{
  public:
    BitVectorEncoder(size_t block_bytes, size_t superblock_size = VectorEncoder::SUPERBLOCK_SIZE);
    ~BitVectorEncoder();

    void setBit(size_t value);
    void setRun(size_t start, size_t len);

    void addBit(size_t value);
    void addRun(size_t start, size_t len);
    void flush();

    // FIXME for gap encoding
    inline void nibbleEncode(size_t diff, size_t len)
    {
      this->size += diff + len - 1;
      this->items += len;
      this->buffer->writeNibbleCode(diff);
      this->buffer->writeNibbleCode(len);
    }

  protected:
    pair_type run;

    // These are not allowed.
    BitVectorEncoder();
    BitVectorEncoder(const BitVectorEncoder&);
    BitVectorEncoder& operator = (const BitVectorEncoder&);
};


/*
  This bit vector uses nibble coding. Each block is either run-length encoded or
  gap encoded, depending on the first nibble.

  // FIXME reverting to gap encoding not implemented yet
*/

class BitVector : public BitVectorBase
{
  public:
    typedef BitVectorEncoder Encoder;

    explicit BitVector(std::ifstream& file);
    explicit BitVector(FILE* file);
    BitVector(Encoder& encoder, size_t universe_size);
    ~BitVector();

//--------------------------------------------------------------------------

    size_t reportSize() const;
    
//--------------------------------------------------------------------------

    class Iterator : public BitVectorBase::Iterator
    {
      public:
        explicit Iterator(const BitVector& par);
        ~Iterator();

        size_t rank(size_t value, bool at_least = false);

        size_t select(size_t index);
        size_t selectNext();

        pair_type valueBefore(size_t value);
        pair_type valueAfter(size_t value);
        pair_type nextValue();

        pair_type selectRun(size_t index, size_t max_length);
        pair_type selectNextRun(size_t max_length);

        bool isSet(size_t value);

        size_t countRuns();

        static const size_t GAP_ENCODING = 0;
        static const size_t RUN_LENGTH_ENCODING = 1;

      protected:

        void valueLoop(size_t value);

        inline void getSample(size_t sample_number)
        {
          BitVectorBase::Iterator::getSample(sample_number);
          this->run = 0;
//           this->use_rle = this->buffer.readNibble();
        }

        bool use_rle;

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

//--------------------------------------------------------------------------
  
  protected:

    // These are not allowed.
    BitVector();
    BitVector(const BitVector&);
    BitVector& operator = (const BitVector&);
};

// Expose the iterator type.
typedef BitVector::Iterator BitVectorIterator;

}

#endif
