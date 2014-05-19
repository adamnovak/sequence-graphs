#ifndef BITBUFFER_H
#define BITBUFFER_H

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <climits>

// Stuff that used to be in RLCSA's definitions.h

typedef std::pair<size_t, size_t> pair_type;

const pair_type EMPTY_PAIR = pair_type(1, 0);

inline size_t popcount(size_t field)
{
  return __builtin_popcountl(field);
}

inline size_t length(size_t n)
{
  size_t b = 0;
  while(n > 0) { b++; n >>= 1; }
  return b;
}

inline size_t nextMultipleOf(size_t multiplier, size_t value)
{
  return multiplier * ((value / multiplier) + 1);
}


const size_t CHARS = ((size_t)1 << CHAR_BIT);
const size_t MEGABYTE = 1048576;
const size_t MILLION  = 1000000;
const size_t WORD_BITS = CHAR_BIT * sizeof(size_t);
const size_t WORD_MAX = ~((size_t)0);

// Previous GET was broken when BITS == WORD_BITS
// Current version works for size_ts and less
//#define GET(FIELD, BITS) ((FIELD) & ((1 << (BITS)) - 1))
#define GET(FIELD, BITS) ((FIELD) & (WORD_MAX >> (WORD_BITS - (BITS))))
#define LOWER(FIELD, N)  ((FIELD) >> (N))
#define HIGHER(FIELD, N) ((FIELD) << (N))

#define BITS_TO_BYTES(BITS) (((BITS) + CHAR_BIT - 1) / CHAR_BIT)
#define BYTES_TO_WORDS(BYTES) (((BYTES) + sizeof(size_t) - 1) / sizeof(size_t))
#define BITS_TO_WORDS(BITS) (((BITS) + WORD_BITS - 1) / WORD_BITS)

class ReadBuffer
{
  public:
    ReadBuffer(std::ifstream& file, size_t words);
    ReadBuffer(FILE* file, size_t words);
    ReadBuffer(std::ifstream& file, size_t _items, size_t item_size);
    ReadBuffer(FILE* file, size_t _items, size_t item_size);

    // These versions do not delete the data when deleted.
    ReadBuffer(const size_t* buffer, size_t words);
    ReadBuffer(const size_t* buffer, size_t _items, size_t item_size);
    ReadBuffer(const ReadBuffer& original);

    ~ReadBuffer();

//--------------------------------------------------------------------------

    void claimData();

    void writeTo(std::ofstream& file) const;
    inline void writeBuffer(std::ofstream& file) const { this->writeTo(file); }
    void writeTo(FILE* file) const;
    inline void writeBuffer(FILE* file) const { this->writeTo(file); }

    // The buffer will no longer own the data.
    void moveBuffer(const size_t* buffer);

    size_t reportSize() const;

//--------------------------------------------------------------------------

    inline void reset()
    {
      this->pos = 0;
      this->bits = WORD_BITS;
      this->current = 0;
    }

    inline void skipBits(size_t count)
    {
      if(count < this->bits)
      {
        this->bits -= count;
        return;
      }

      count -= this->bits;
      this->pos += 1 + count / WORD_BITS;
      this->bits = WORD_BITS - count % WORD_BITS;
    }

//--------------------------------------------------------------------------

    inline size_t bitsLeft() const
    {
      return this->bits + WORD_BITS * (this->size - this->pos - 1);
    }

    // Returns nonzero if bit is 1
    inline size_t isSet(size_t index) const
    {
      return this->data[index / WORD_BITS] & ((size_t)1 << (WORD_BITS - index % WORD_BITS - 1));
    }

    // Returns nonzero if bit is 1
    inline size_t readBit()
    {
      this->bits--;
      size_t bit = this->data[this->pos] & ((size_t)1 << this->bits);

      if(this->bits == 0) { this->pos++; this->bits = WORD_BITS; }

      return bit;
    }

    inline size_t readBits(size_t count)
    {
      size_t value = 0;

      if(count >= this->bits)
      {
        count -= this->bits;
        value |= HIGHER(GET(this->data[this->pos], this->bits), count);
        this->pos++; this->bits = WORD_BITS;
      }
      if(count > 0)
      {
        this->bits -= count;
        value |= GET(LOWER(this->data[this->pos], this->bits), count);
      }

      return value;
    }

//--------------------------------------------------------------------------

    /*
      These operations work on 4-bit nibbles.
      Do not use these with the bit-level operations.
    */

    inline size_t readNibble()
    {
      this->bits -= 4;
      size_t value = GET(LOWER(this->data[this->pos], this->bits), 4);

      if(this->bits == 0) { this->pos++; this->bits = WORD_BITS; }

      return value;
    }

    // Nibble code for positive integers.
    inline size_t readNibbleCode()
    {
      size_t temp, value = 0, shift = 0;
      do
      {
        temp = this->readNibble();
        value |= (temp & 0x7) << shift;
        shift += 3;
      }
      while((temp & 0x8) == 0);

      return value + 1;
    }

    // This version reads the code only if value <= limit.
    inline size_t readNibbleCode(size_t limit)
    {
       size_t _pos = this->pos, _bits = this->bits;
       size_t value = this->readNibbleCode();
       if(value > limit)
       {
         this->pos = _pos; this->bits = _bits;
         return 0;
       }
       return value;
    }

//--------------------------------------------------------------------------

    /*
      These operations work on fixed-size items. No sanity checks are made
      for parameter values.
    */

    inline size_t getItemSize() const
    {
      return this->item_bits;
    }

    inline size_t getNumberOfItems() const
    {
      return this->items;
    }

    inline void goToItem(size_t item)
    {
      size_t b = item * this->item_bits;
      this->pos = b / WORD_BITS;
      this->bits = WORD_BITS - b % WORD_BITS;
      this->current = item;
    }

    inline size_t nextItem()
    {
      this->current++;
      return this->readBits(this->item_bits);
    }

    inline size_t readItem() { return this->nextItem(); }

    inline size_t readItem(size_t item)
    {
      this->goToItem(item);
      return this->nextItem();
    }

    inline size_t readFirstItem()
    {
      return this->readItem(0);
    }

    inline size_t readItemConst(size_t item) const
    {
      size_t b = item * this->item_bits;
      size_t p = b / WORD_BITS;
      b = WORD_BITS - b % WORD_BITS;

      size_t c = this->item_bits;
      size_t value = 0;

      if(c >= b)
      {
        c -= b;
        value |= HIGHER(GET(this->data[p], b), c);
        p++; b = WORD_BITS;
      }
      if(c > 0)
      {
        b -= c;
        value |= GET(LOWER(this->data[p], b), c);
      }

      return value;
    }

    inline bool hasNextItem() const
    {
      return (this->current < this->items);
    }

    inline void skipItem()
    {
      this->skipBits(this->item_bits);
      this->current++;
    }

//--------------------------------------------------------------------------

    /*
      Delta coding for positive integers
    */

    inline size_t readDeltaCode()
    {
      size_t len = 0;
      while(this->readBit() == 0) { len++; }

      size_t temp = (((size_t)1 << len) | this->readBits(len)) - 1;
      temp = ((size_t)1 << temp) | this->readBits(temp);
      return temp;
    }

    // This version reads the code only if value <= limit.
    inline size_t readDeltaCode(size_t limit)
    {
       size_t _pos = this->pos, _bits = this->bits;
       size_t value = this->readDeltaCode();
       if(value > limit)
       {
         this->pos = _pos; this->bits = _bits;
         return 0;
       }
       return value;
    }

//--------------------------------------------------------------------------

  private:
    const size_t* data;
    size_t size, item_bits, items;
    bool  free_buffer;

    // Iterator data
    size_t pos, bits, current;

    inline static size_t bitsToWords(size_t _bits) { return (_bits + WORD_BITS - 1) / WORD_BITS; }

    // These are not allowed.
    ReadBuffer();
    ReadBuffer& operator = (const ReadBuffer&);
};


//--------------------------------------------------------------------------


class WriteBuffer
{
  public:
    explicit WriteBuffer(size_t words);
    WriteBuffer(size_t _items, size_t item_size);

    // These versions do not delete the data when deleted.
    WriteBuffer(size_t* buffer, size_t words);
    WriteBuffer(size_t* buffer, size_t _items, size_t item_size);

    ~WriteBuffer();

//--------------------------------------------------------------------------

    // This transfers the ownership of the data to the read buffer.
    ReadBuffer* getReadBuffer();

    void writeTo(std::ofstream& file) const;
    inline void writeBuffer(std::ofstream& file) const { this->writeTo(file); }
    void writeTo(FILE* file) const;
    inline void writeBuffer(FILE* file) const { this->writeTo(file); }

    // The buffer will no longer own the data.
    void moveBuffer(size_t* buffer);

    size_t reportSize() const;

//--------------------------------------------------------------------------

    inline void reset()
    {
      this->pos = 0;
      this->bits = WORD_BITS;
      this->current = 0;
    }

    inline void skipBits(size_t count)
    {
      if(count < this->bits)
      {
        this->bits -= count;
        return;
      }

      count -= this->bits;
      this->pos += 1 + count / WORD_BITS;
      this->bits = WORD_BITS - count % WORD_BITS;
    }

//--------------------------------------------------------------------------

    inline size_t bitsLeft() const
    {
      return this->bits + WORD_BITS * (this->size - this->pos - 1);
    }

    inline void writeBits(size_t value, size_t count)
    {
      if(count >= this->bits)
      {
        count -= this->bits;
        this->data[this->pos] |= GET(LOWER(value, count), this->bits);
        this->pos++; this->bits = WORD_BITS;
      }
      if(count > 0)
      {
        this->bits -= count;
        this->data[this->pos] |= HIGHER(GET(value, count), this->bits);
      }
    }

//--------------------------------------------------------------------------

    // Returns nonzero if bit is 1
    inline size_t isSet(size_t index) const
    {
      return this->data[index / WORD_BITS] & ((size_t)1 << (WORD_BITS - index % WORD_BITS - 1));
    }

    inline void setBit(size_t index)
    {
      this->data[index / WORD_BITS] |= (size_t)1 << (WORD_BITS - index % WORD_BITS - 1);
    }

    inline void unsetBit(size_t index)
    {
      this->data[index / WORD_BITS] &= ~((size_t)1 << (WORD_BITS - index % WORD_BITS - 1));
    }

//--------------------------------------------------------------------------

    /*
      These operations work on fixed-size items.
    */

    inline size_t getItemSize() const
    {
      return this->item_bits;
    }

    inline void goToItem(size_t item)
    {
      size_t b = item * this->item_bits;
      this->pos = b / WORD_BITS;
      this->bits = WORD_BITS - b % WORD_BITS;
      this->current = item;
    }

    inline bool hasNextItem() const
    {
      return (this->current < this->items);
    }

    inline void writeItem(size_t item)
    {
      this->writeBits(item, this->item_bits);
      this->current++;
    }

    inline void skipItem()
    {
      this->skipBits(this->item_bits);
      this->current++;
    }

//--------------------------------------------------------------------------

    /*
      Nibble coding for positive integers.
    */

    inline size_t nibbleCodeLength(size_t value) const
    {
      size_t b = 0;
      value--;

      do
      {
        b += 4;
        value >>= 3;
      }
      while(value > 0);

      return b;
    }

    // Something breaks very badly if value > 15.
    inline void writeNibble(size_t value)
    {
      this->bits -= 4;
      this->data[this->pos] |= HIGHER(value, this->bits);
      if(this->bits == 0) { this->pos++; this->bits = WORD_BITS; }
    }

    // It is assumed that there is enough space for the code.
    inline void writeNibbleCode(size_t value)
    {
      value--;
      while(value > 0x7)
      {
        this->writeNibble(value & 0x7);
        value >>= 3;
      }
      this->writeNibble(value | 0x8);
    }

//--------------------------------------------------------------------------

    /*
      Delta coding for positive integers
    */

    inline bool canDeltaCode(size_t value) const
    {
      return this->deltaCodeLength(value) <= this->bitsLeft();
    }

    inline size_t deltaCodeLength(size_t value) const
    {
      size_t len = length(value);
      size_t llen = length(len);
      return (len + llen + llen - 2);
    }

    // This version returns false if there is no space left for the encoding.
    inline bool writeDeltaCode(size_t value)
    {
      size_t len = length(value);
      size_t llen = length(len);

      if(len + llen + llen - 2 > this->bitsLeft()) { return false; }

      // this->writeBits(0, llen - 1); // Now included in the next writeBits()
      this->writeBits(len, llen + llen - 1);
      this->writeBits(value, len - 1);
      return true;
    }

    // This version assumes the code fits into the buffer.
    inline void writeDeltaCodeDirect(size_t value)
    {
      size_t len = length(value);
      size_t llen = length(len);

      // this->writeBits(0, llen - 1); // Now included in the next writeBits()
      this->writeBits(len, llen + llen - 1);
      this->writeBits(value, len - 1);
    }

    // We assume the code fits into size_t:
    //  32-bit:  value < 2^24
    //  64-bit:  value < 2^54
    inline void writeDeltaCodeFast(size_t value)
    {
      size_t len = length(value);

      value ^= ((size_t)1 << (len - 1));
      this->writeBits((len << (len - 1)) | value, len + 2 * length(len) - 2);
    }

//--------------------------------------------------------------------------

  private:
    size_t* data;
    size_t size, item_bits, items;
    bool free_buffer;

    // Iterator data
    size_t pos, bits, current;

    inline static size_t bitsToWords(size_t _bits) { return (_bits + WORD_BITS - 1) / WORD_BITS; }

    // These are not allowed.
    WriteBuffer();
    WriteBuffer(const WriteBuffer&);
    WriteBuffer& operator = (const WriteBuffer&);
};


//--------------------------------------------------------------------------


#endif // BITBUFFER_H
