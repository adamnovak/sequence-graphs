#ifndef CSA_NIBBLEVECTOR_H
#define CSA_NIBBLEVECTOR_H

#include <fstream>

#include "BitVectorBase.hpp"

namespace CSA {

/*
  This class is used to construct a NibbleVector.
*/

class NibbleEncoder : public VectorEncoder
{
  public:
    NibbleEncoder(size_t block_bytes, size_t superblock_size = VectorEncoder::SUPERBLOCK_SIZE);
    ~NibbleEncoder();

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
    NibbleEncoder();
    NibbleEncoder(const NibbleEncoder&);
    NibbleEncoder& operator = (const NibbleEncoder&);
};


/*
  This bit vector uses nibble coding. Each block is either run-length encoded or
  gap encoded, depending on the first nibble.

  // FIXME reverting to gap encoding not implemented yet
*/

class NibbleVector : public BitVectorBase
{
  public:
    typedef NibbleEncoder Encoder;

    explicit NibbleVector(std::ifstream& file);
    explicit NibbleVector(FILE* file);
    NibbleVector(Encoder& encoder, size_t universe_size);
    ~NibbleVector();

//--------------------------------------------------------------------------

    size_t reportSize() const;

//--------------------------------------------------------------------------

    class Iterator : public BitVectorBase::Iterator
    {
      public:
        explicit Iterator(const NibbleVector& par);
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
    NibbleVector();
    NibbleVector(const NibbleVector&);
    NibbleVector& operator = (const NibbleVector&);
};

}

#endif