#ifndef CSA_BITVECTOR_H
#define CSA_BITVECTOR_H

#include <cstdio>
#include <fstream>
#include <list>

#include "BitBuffer.hpp"

namespace CSA {

/**
 * This class provides the core functionality for encoding a bit vector. It is 
 * abstract and cannot be directly instantiated.
*/

class VectorEncoder
{
  public:
    static const size_t SUPERBLOCK_SIZE = 1024;

    // We assume superblock size is divisible by block and sample size.
    VectorEncoder(size_t block_bytes, size_t superblock_size = SUPERBLOCK_SIZE,
        bool _use_small_blocks = true);
    
    ~VectorEncoder();

    /*
     * These must be implemented in any inherited class.
     */

    // These functions are assumed to be greedy, encoding the 1-bits
    // immediately.
    
    // Values must be in increasing order.
    virtual void setBit(size_t value) = 0;  
    virtual void setRun(size_t start, size_t len) = 0;

    // These versions may combine 1-bits into maximal runs.
    // Use flush() to finish the encoding when using these functions.
    // Do not mix with the greedy versions.
    virtual void addBit(size_t value) = 0;
    virtual void addRun(size_t start, size_t len) = 0;
    virtual void flush() = 0;


    void addNewBlock();
    void setFirstBit(size_t value);

    size_t size, items, blocks;
    size_t block_size, superblock_bytes;
    bool  use_small_blocks;

    WriteBuffer*      buffer;

    std::list<size_t*> array_blocks;
    size_t*            array;
    size_t             blocks_in_superblock, current_blocks;

    std::list<size_t*> sample_blocks;
    size_t*            samples;
    size_t             samples_in_superblock, current_samples;

  protected:

    // These are not allowed.
    VectorEncoder();
    VectorEncoder(const VectorEncoder&);
    VectorEncoder& operator = (const VectorEncoder&);
};


/*
  This class provides the core functionality for a bit vector.
  A bit vector must have at least one 1-bit. This class is abstract.
*/

class BitVectorBase
{
  public:
    static const size_t INDEX_RATE = 5;

    explicit BitVectorBase(std::ifstream& file);
    explicit BitVectorBase(FILE* file);
    BitVectorBase(VectorEncoder& encoder, size_t universe_size);
    explicit BitVectorBase(WriteBuffer& vector);
    ~BitVectorBase();

//--------------------------------------------------------------------------

    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;

    inline size_t getSize() const { return this->size; }
    inline size_t getNumberOfItems() const { return this->items; }
    inline size_t getBlockSize() const { return this->block_size; }

    // This returns only the sizes of the dynamically allocated structures.
    size_t reportSize() const;

    size_t getCompressedSize() const;

    // Removes structures not necessary for merging.
    void strip();

//--------------------------------------------------------------------------

    class Iterator
    {
      public:
        explicit Iterator(const BitVectorBase& par);
        ~Iterator();

        inline bool hasNext() const
        {
          return (this->sample.first + this->cur < this->parent.items - 1);
        }
        
        /*
         * These functions must be implemented by derived classes.
         */
        
        // rank invalidates the "next" functionality
        // regular:   \sum_{i = 0}^{value} V[i]
        // at_least:  \sum_{i = 0}^{value - 1} V[i] + 1
        virtual size_t rank(size_t value, bool at_least = false) = 0;

        // \min value: \sum_{i = 0}^{value} V[i] = index + 1
        virtual size_t select(size_t index) = 0;      
        virtual size_t selectNext() = 0;

        // (\max i <= value: V[i] = 1, rank(i) - 1)
        // Returns (size, items) if not found.
        virtual pair_type valueBefore(size_t value) = 0;

        // (\min i >= value: V[i] = 1, rank(i) - 1)
        virtual pair_type valueAfter(size_t value) = 0; 
        virtual pair_type nextValue() = 0;

        // These versions of select return (value, length_of_run).
        // max_length is an upper bound for the length of the run returned.
        // V[value] is not included in the length of the run
        // These functions are not greedy: the actual length of the run can be
        // more than reported.
        // This can happen even if max_length was not reached.
        // length_of_run is actually the number of extra items returned past
        // value

        virtual pair_type selectRun(size_t index, size_t max_length) = 0;
        virtual pair_type selectNextRun(size_t max_length) = 0;

        // isSet invalidates the "next" functionality
        virtual bool isSet(size_t value) = 0; // V[value]

        // Counts the number of 1-bit runs.
        virtual size_t countRuns() = 0;

      protected:
        const BitVectorBase& parent;

        size_t      block;
        pair_type  sample;
        size_t      cur, val, run; // cur == 0 is the sample
        size_t      block_items;

        ReadBuffer buffer, samples;

        /*
          These functions return the sample corresponding to the block the given
          index/value might be found in. Parameters are assumed to be valid.
        */
        size_t sampleForIndex(size_t index);
        size_t sampleForValue(size_t value);

        inline size_t getSampledIndex(size_t sample_number)
        {
          return this->samples.readItem(2 * sample_number);
        }

        inline size_t getSampledValue(size_t sample_number)
        {
          return this->samples.readItem(2 * sample_number + 1);
        }

        inline void getSample(size_t sample_number)
        {
          this->block = sample_number;
          this->samples.goToItem(2 * sample_number);
          this->sample.first = this->samples.readItem();
          this->sample.second = this->samples.readItem();
          this->cur = 0;
          this->val = this->sample.second;
          this->block_items = this->samples.readItem() - this->sample.first - 1;
          this->buffer.moveBuffer(this->parent.array + (this->block * this->parent.block_size));
        }

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

//--------------------------------------------------------------------------

  protected:
    size_t size, items;

    const size_t* array;
    size_t        block_size;
    size_t        number_of_blocks;

    /*
      Each sample is (rank(i) - 1, i) where V[i] = 1.
      Number of samples is number_of_blocks + 1.
    */
    ReadBuffer*  samples;
    size_t        integer_bits;

    ReadBuffer*  rank_index;
    size_t        rank_rate;

    ReadBuffer*  select_index;
    size_t        select_rate;

    /*
       These functions build a higher level index for faster rank/select
       queries. The index consists of about (number of samples) / INDEX_RATE 
       pointers. The bit vector cannot be used without the index.
    */
    void indexForRank();
    void indexForSelect();

    // These are used in disk storage.
    void writeHeader(std::ofstream& file) const;
    void writeHeader(FILE* file) const;
    void writeArray(std::ofstream& file) const;
    void writeArray(FILE* file) const;
    void readHeader(std::ifstream& file);
    void readHeader(FILE* file);
    void readArray(std::ifstream& file);
    void readArray(FILE* file);

    void copyArray(VectorEncoder& encoder, bool use_directly = false);
private:
    // These are not allowed.
    BitVectorBase();
    BitVectorBase(const BitVectorBase&);
    BitVectorBase& operator = (const BitVectorBase&);
};

}

#endif 
