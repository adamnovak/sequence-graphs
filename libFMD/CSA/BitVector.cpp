#include <cstdlib>

#include <stdexcept>

#include "BitVector.hpp"

namespace CSA {

BitVector::BitVector(std::ifstream& file) :
  BitVectorBase(file)
{
}

BitVector::BitVector(FILE* file) :
  BitVectorBase(file)
{
}

BitVector::BitVector(Encoder& encoder, size_t universe_size) :
  BitVectorBase(encoder, universe_size)
{
}

BitVector::~BitVector()
{
}

//--------------------------------------------------------------------------

size_t
BitVector::reportSize() const
{
  size_t bytes = sizeof(*this);
  bytes += BitVectorBase::reportSize();
  return bytes;
}

BitVector* 
BitVector::createUnion(const BitVector& other) const
{
  // Make iterators to actually look at the vectors
  BitVectorIterator us(*this);
  BitVectorIterator them(other);

  // How many bits are in the longer BitVector? That's how long we will need to
  // make the union.
  size_t newSize = std::max(getSize(), other.getSize());
  
  // We need an encoder to encode it.
  BitVectorEncoder encoder(block_size);
  
  // How many ones are we putting in?
  size_t ones = 0;
  
  for(size_t i = 0; i < newSize; i++)
  {
    // OR each bit
    if(us.isSet(i) || them.isSet(i)) {
        // If the bit is true, add a bit at this position.
        encoder.addBit(i);
        ones++;
    }
  }
  
  // Finish encoding.
  encoder.flush();
  
  // Build the BitVector
  BitVector* toReturn = new BitVector(encoder, newSize);
  
  // Check its bits
  BitVectorIterator iterator(*toReturn);
  
  if(iterator.rank(newSize) != ones) {
    // Check the number of 1s to make sure we did it right.
    throw std::runtime_error("Expected " + std::to_string(ones) + 
        " ones in union but found " + std::to_string(iterator.rank(newSize)));
  }
  
  // Return the BitVector. Caller is responsible for it.
  return toReturn;
  
}

//--------------------------------------------------------------------------

BitVector::Iterator::Iterator(const BitVector& par) :
  BitVectorBase::Iterator(par),
  use_rle(false)
{
}

BitVector::Iterator::~Iterator()
{
}

// FIXME gap encoding for all operations

size_t
BitVector::Iterator::rank(size_t value, bool at_least)
{
  const BitVector& par = (const BitVector&)(this->parent);

  if(value >= par.size) { return par.items; }

  this->valueLoop(value);

  size_t idx = this->sample.first + this->cur + 1;
  if(!at_least && this->val > value)
  {
    idx--;
  }
  if(at_least && this->val < value)
  {
    this->getSample(this->block + 1);
    idx = this->sample.first + this->cur + 1;
  }
  return idx;
}

size_t
BitVector::Iterator::select(size_t index)
{
  const BitVector& par = (const BitVector&)(this->parent);

  if(index >= par.items) { return par.size; }
  this->getSample(this->sampleForIndex(index));

  size_t lim = index - this->sample.first;
  while(this->cur < lim)
  {
    this->val += this->buffer.readNibbleCode();
    size_t temp = this->buffer.readNibbleCode();
    this->val += temp - 1;
    this->cur += temp;
  }
  if(this->cur > lim)
  {
    this->run = this->cur - lim;
    this->cur -= this->run;
    this->val -= this->run;
  }

  return this->val;
}

size_t
BitVector::Iterator::selectNext()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return this->val;
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer.readNibbleCode();
    this->run = this->buffer.readNibbleCode() - 1;
  }

  return this->val;
}

pair_type
BitVector::Iterator::valueBefore(size_t value)
{
  const BitVector& par = (const BitVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }

  this->getSample(this->sampleForValue(value));
  if(this->val > value) { return pair_type(par.size, par.items); }
  this->run = 0;

  while(this->cur < this->block_items && this->val < value)
  {
    size_t temp = this->buffer.readNibbleCode(value - this->val);
    if(temp == 0) { break; }
    this->val += temp;

    temp = this->buffer.readNibbleCode();
    this->cur += temp;
    this->val += temp - 1;
  }
  if(this->val > value)
  {
    this->run = this->val - value;
    this->val = value;
    this->cur -= this->run;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
BitVector::Iterator::valueAfter(size_t value)
{
  const BitVector& par = (const BitVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }

  this->valueLoop(value);

  if(this->val < value)
  {
    this->getSample(this->block + 1);
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
BitVector::Iterator::nextValue()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return pair_type(this->val, this->sample.first);
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer.readNibbleCode();
    this->run = this->buffer.readNibbleCode() - 1;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
BitVector::Iterator::selectRun(size_t index, size_t max_length)
{
  size_t value = this->select(index);

  size_t len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

pair_type
BitVector::Iterator::selectNextRun(size_t max_length)
{
  size_t value = this->selectNext();

  size_t len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

bool
BitVector::Iterator::isSet(size_t value)
{
  const BitVector& par = (const BitVector&)(this->parent);

  if(value >= par.size) { return false; }

  this->valueLoop(value);

  return (this->val == value);
}

size_t
BitVector::Iterator::countRuns()
{
  const BitVector& par = (const BitVector&)(this->parent);

  if(par.items == 0) { return 0; }

  size_t runs = 1;
  pair_type res = this->selectRun(0, par.items);
  size_t last = res.first + res.second;

  while(last < par.size)
  {
    res = this->selectNextRun(par.items);
    if(res.first < par.size && res.first > last + 1) { runs++; }
    last = res.first + res.second;
  }

  return runs;
}

// FIXME for gap encoding
void
BitVector::Iterator::valueLoop(size_t value)
{
  this->getSample(this->sampleForValue(value));

  if(this->val >= value) { return; }
  while(this->cur < this->block_items)
  {
    this->val += this->buffer.readNibbleCode();
    this->cur++;
    this->run = this->buffer.readNibbleCode() - 1;
    if(this->val >= value) { break; }

    this->cur += this->run;
    this->val += this->run;
    if(this->val >= value)
    {
      this->run = this->val - value;
      this->val = value;
      this->cur -= this->run;
      break;
    }
    this->run = 0;
  }
}

//--------------------------------------------------------------------------

BitVectorEncoder::BitVectorEncoder(size_t block_bytes, size_t superblock_size) :
  VectorEncoder(block_bytes, superblock_size),
  run(EMPTY_PAIR)
{
}

BitVectorEncoder::~BitVectorEncoder()
{
}

void
BitVectorEncoder::setBit(size_t value)
{
  this->setRun(value, 1);
}

// FIXME for gap encoding
void
BitVectorEncoder::setRun(size_t start, size_t len)
{
  if(this->items == 0)
  {
    this->setFirstBit(start);
    if(len > 1)
    {
      this->nibbleEncode(1, len - 1);
    }
    return;
  }
  if(start < this->size || len == 0) { return; }

  // Write as much into the buffer as possible.
  size_t diff = start + 1 - this->size;
  size_t free_bits = this->buffer->bitsLeft();
  size_t code_bits = this->buffer->nibbleCodeLength(diff);
  if(free_bits >= code_bits + 4) // At least a part of the run fits into the block.
  {
    free_bits -= code_bits;
    size_t run_bits = this->buffer->nibbleCodeLength(len);
    if(run_bits <= free_bits)
    {
      this->nibbleEncode(diff, len);
      return;
    }

    // Encode as much as possible and let the rest spill.
    size_t llen = (size_t)1 << (3 * (free_bits / 4));
    this->nibbleEncode(diff, llen);
    len -= llen;

    // A new sample will be added.
    this->size++;
    this->items++;
  }
  else
  {
    this->size = start + 1;
    this->items++;
  }

  // Didn't fit into the block. A new sample & block required.
  this->addNewBlock();
  if(len > 1)
  {
    this->nibbleEncode(1, len - 1);
  }
}

void
BitVectorEncoder::addBit(size_t value)
{
  this->addRun(value, 1);
}

void
BitVectorEncoder::addRun(size_t start, size_t len)
{
  if(this->run.second == 0)
  {
    this->run = pair_type(start, len);
  }
  else if(start == this->run.first + this->run.second)
  {
    this->run.second += len;
  }
  else
  {
    this->setRun(this->run.first, this->run.second);
    this->run = pair_type(start, len);
  }
}

void
BitVectorEncoder::flush()
{
  this->setRun(this->run.first, this->run.second);
  this->run.second = 0;
}

}
