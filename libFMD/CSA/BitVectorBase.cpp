#include <cstring>
#include <cstdlib>

#include "BitVector.hpp"

namespace CSA {

BitVectorBase::BitVectorBase(std::ifstream& file) :
  rank_index(0), select_index(0)
{
  this->readHeader(file);
  this->readArray(file);

  this->integer_bits = length(this->size);
  this->samples = new ReadBuffer(file, 2 * (this->number_of_blocks + 1), this->integer_bits);

  this->indexForRank();
  this->indexForSelect();
}

BitVectorBase::BitVectorBase(FILE* file) :
  rank_index(0), select_index(0)
{
  this->readHeader(file);
  this->readArray(file);

  this->integer_bits = length(this->size);
  this->samples = new ReadBuffer(file, 2 * (this->number_of_blocks + 1), this->integer_bits);

  this->indexForRank();
  this->indexForSelect();
}

BitVectorBase::BitVectorBase(VectorEncoder& encoder, size_t universe_size) :
  size(universe_size), items(encoder.items),
  block_size(encoder.block_size),
  number_of_blocks(encoder.blocks),
  rank_index(0), select_index(0)
{
  if(this->items == 0)
  {
    std::cerr << "BitVector: Cannot create a bit vector with no 1-bits!" << std::endl;
    return;
  }
  this->copyArray(encoder);

  this->integer_bits = length(this->size);
  WriteBuffer sample_buffer(2 * (this->number_of_blocks + 1), this->integer_bits);

  // Compress the samples.
  for(std::list<size_t*>::iterator iter = encoder.sample_blocks.begin(); iter != encoder.sample_blocks.end(); iter++)
  {
    size_t* buf = *iter;
    for(size_t i = 0; i < 2 * encoder.samples_in_superblock; i++)
    {
      sample_buffer.writeItem(buf[i]);
    }
  }
  for(size_t i = 0; i < 2 * encoder.current_samples; i++)
  {
    sample_buffer.writeItem(encoder.samples[i]);
  }
  sample_buffer.writeItem(this->items);
  sample_buffer.writeItem(this->size);

  this->samples = sample_buffer.getReadBuffer();

  this->indexForRank();
  this->indexForSelect();
}

BitVectorBase::BitVectorBase() :
  array(0), samples(0), rank_index(0), select_index(0)
{
}

BitVectorBase::~BitVectorBase()
{
  delete[] this->array;
  delete this->samples;
  delete this->rank_index;
  delete this->select_index;
}

//--------------------------------------------------------------------------

void
BitVectorBase::writeTo(std::ofstream& file) const
{
  this->writeHeader(file);
  this->writeArray(file);
  this->samples->writeBuffer(file);
}

void
BitVectorBase::writeTo(FILE* file) const
{
  this->writeHeader(file);
  this->writeArray(file);
  this->samples->writeBuffer(file);
}

void
BitVectorBase::writeHeader(std::ofstream& file) const
{
  file.write((char*)&(this->size), sizeof(this->size));
  file.write((char*)&(this->items), sizeof(this->items));
  file.write((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.write((char*)&(this->block_size), sizeof(this->block_size));
}

void
BitVectorBase::writeHeader(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(&(this->size), sizeof(this->size), 1, file);
  std::fwrite(&(this->items), sizeof(this->items), 1, file);
  std::fwrite(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file);
  std::fwrite(&(this->block_size), sizeof(this->block_size), 1, file);
}

void
BitVectorBase::writeArray(std::ofstream& file) const
{
  file.write((char*)(this->array), this->block_size * this->number_of_blocks * sizeof(size_t));
}

void
BitVectorBase::writeArray(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(this->array, this->block_size * sizeof(size_t), this->number_of_blocks, file);
}

void
BitVectorBase::readHeader(std::ifstream& file)
{
  file.read((char*)&(this->size), sizeof(this->size));
  file.read((char*)&(this->items), sizeof(this->items));
  file.read((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.read((char*)&(this->block_size), sizeof(this->block_size));
}

void
BitVectorBase::readHeader(FILE* file)
{
  if(file == 0) { return; }
  if(!std::fread(&(this->size), sizeof(this->size), 1, file)) { return; }
  if(!std::fread(&(this->items), sizeof(this->items), 1, file)) { return; }
  if(!std::fread(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file)) { return; }
  if(!std::fread(&(this->block_size), sizeof(this->block_size), 1, file)) { return; }
}

void
BitVectorBase::readArray(std::ifstream& file)
{
  size_t* array_buffer = new size_t[this->block_size * this->number_of_blocks];
  file.read((char*)(array_buffer), this->block_size * this->number_of_blocks * sizeof(size_t));
  this->array = array_buffer;
}

void
BitVectorBase::readArray(FILE* file)
{
  if(file == 0) { return; }
  size_t* array_buffer = new size_t[this->block_size * this->number_of_blocks];
  if(!std::fread(array_buffer, this->block_size * sizeof(size_t), this->number_of_blocks, file)) { return; }
  this->array = array_buffer;
}

//--------------------------------------------------------------------------

void
BitVectorBase::copyArray(VectorEncoder& encoder, bool use_directly)
{
  if(use_directly && encoder.array_blocks.size() == 0)
  {
    this->array = encoder.array; encoder.array = 0;
    return;
  }

  size_t* array_buffer = new size_t[this->block_size * this->number_of_blocks];
  size_t pos = 0, total_words = this->block_size * this->number_of_blocks;

  for(std::list<size_t*>::iterator iter = encoder.array_blocks.begin(); iter != encoder.array_blocks.end(); iter++)
  {
    memcpy(array_buffer + pos, *iter, encoder.superblock_bytes);
    pos += encoder.superblock_bytes / sizeof(size_t);
  }

  size_t remainder = std::min(total_words - pos, (size_t)(encoder.superblock_bytes / sizeof(size_t)));
  memcpy(array_buffer + pos, encoder.array, remainder * sizeof(size_t));
  pos += remainder;
  if(pos < total_words) { memset(array_buffer + pos, 0, (total_words - pos) * sizeof(size_t)); }

  this->array = array_buffer;
}

//--------------------------------------------------------------------------

size_t
BitVectorBase::reportSize() const
{
  // We assume the reportSize() of derived classes includes any class variables of BitVector.
  size_t bytes = this->block_size * this->number_of_blocks * sizeof(size_t);
  if(this->samples != 0) { bytes += this->samples->reportSize(); }
  if(this->rank_index != 0) { bytes += this->rank_index->reportSize(); }
  if(this->select_index != 0) { bytes += this->select_index->reportSize(); }
  return bytes;
}

size_t
BitVectorBase::getCompressedSize() const
{
  return this->block_size * this->number_of_blocks * sizeof(size_t);
}

//--------------------------------------------------------------------------

void
BitVectorBase::strip()
{
  delete this->rank_index; this->rank_index = 0;
}

//--------------------------------------------------------------------------

void
BitVectorBase::indexForRank()
{
  delete this->rank_index;

  size_t value_samples = (this->number_of_blocks + BitVectorBase::INDEX_RATE - 1) / BitVectorBase::INDEX_RATE;
  this->rank_rate = (this->size + value_samples - 1) / value_samples;
  value_samples = (this->size + this->rank_rate - 1) / this->rank_rate + 1;
  WriteBuffer index_buffer(value_samples, length(this->number_of_blocks - 1));

  // current is value, pointer is sample number.
  size_t current = 0, pointer = 0;
  this->samples->goToItem(2);
  while(this->samples->hasNextItem())
  {
    this->samples->skipItem();
    size_t limit = this->samples->readItem();  // Next sampled value.
    while(current < limit)
    {
      index_buffer.writeItem(pointer);
      current += this->rank_rate;
    }
    pointer++;
  }
  index_buffer.writeItem(this->number_of_blocks - 1);

  this->rank_index = index_buffer.getReadBuffer();
}

void
BitVectorBase::indexForSelect()
{
  delete this->select_index;

  size_t index_samples = (this->number_of_blocks + BitVectorBase::INDEX_RATE - 1) / BitVectorBase::INDEX_RATE;
  this->select_rate = (this->items + index_samples - 1) / index_samples;
  index_samples = (this->items + this->select_rate - 1) / this->select_rate + 1;
  WriteBuffer index_buffer(index_samples, length(this->number_of_blocks - 1));

  // current is index, pointer is sample number.
  size_t current = 0, pointer = 0;
  this->samples->goToItem(2);
  while(this->samples->hasNextItem())
  {
    size_t limit = this->samples->readItem();  // Next sampled index.
    this->samples->skipItem();
    while(current < limit)
    {
      index_buffer.writeItem(pointer);
      current += this->select_rate;
    }
    pointer++;
  }
  index_buffer.writeItem(this->number_of_blocks - 1);

  this->select_index = index_buffer.getReadBuffer();
}

//--------------------------------------------------------------------------

BitVectorBase::Iterator::Iterator(const BitVectorBase& par) :
  parent(par),
  buffer(par.array, par.block_size),
  samples(*(par.samples))
{
}

BitVectorBase::Iterator::~Iterator()
{
}

size_t
BitVectorBase::Iterator::sampleForIndex(size_t index)
{
  size_t low = this->parent.select_index->readItemConst(index / this->parent.select_rate);
  size_t high = this->parent.number_of_blocks - 1;

  this->samples.goToItem(2 * low + 2);
  for(; low < high; low++)
  {
    if(this->samples.readItem() > index) { return low; }
    this->samples.skipItem();
  }

  return low;
}

size_t
BitVectorBase::Iterator::sampleForValue(size_t value)
{
  size_t low = this->parent.rank_index->readItemConst(value / this->parent.rank_rate);
  size_t high = this->parent.number_of_blocks - 1;

  this->samples.goToItem(2 * low + 3);
  for(; low < high; low++)
  {
    if(this->samples.readItem() > value) { return low; }
    this->samples.skipItem();
  }

  return low;
}

//--------------------------------------------------------------------------

VectorEncoder::VectorEncoder(size_t block_bytes, size_t superblock_size, bool _use_small_blocks) :
  size(0), items(0), blocks(0),
  block_size(BYTES_TO_WORDS(block_bytes)),
  superblock_bytes(superblock_size),
  use_small_blocks(_use_small_blocks),
  blocks_in_superblock(1), current_blocks(0),
  samples(0), samples_in_superblock(0), current_samples(0)
{
  this->array = new size_t[this->superblock_bytes / sizeof(size_t)];
  memset(this->array, 0, this->superblock_bytes);
  if(use_small_blocks)
  {
    this->blocks_in_superblock = this->superblock_bytes / (sizeof(size_t) * this->block_size);
    this->buffer = new WriteBuffer(this->array, this->block_size);
    this->samples = new size_t[this->superblock_bytes / sizeof(size_t)];
    this->samples_in_superblock = this->superblock_bytes / (2 * sizeof(size_t));
  }
  else
  {
    this->buffer = new WriteBuffer(this->array, BYTES_TO_WORDS(this->superblock_bytes));
    this->blocks = this->current_blocks = 1;
  }
}

VectorEncoder::~VectorEncoder()
{
  delete[] this->array;

  delete this->buffer;
  for(std::list<size_t*>::iterator iter = this->array_blocks.begin(); iter != this->array_blocks.end(); iter++)
  {
    delete[] *iter;
  }

  delete[] this->samples;
  for(std::list<size_t*>::iterator iter = this->sample_blocks.begin(); iter != this->sample_blocks.end(); iter++)
  {
    delete[] *iter;
  }
}

void
VectorEncoder::addNewBlock()
{
  // Do we need a new superblock for the block?
  this->blocks++; this->current_blocks++;
  if(this->current_blocks > this->blocks_in_superblock)
  {
    this->array_blocks.push_back(this->array);
    this->array = new size_t[this->superblock_bytes / sizeof(size_t)];
    memset(this->array, 0, this->superblock_bytes);
    this->current_blocks = 1;
  }
  this->buffer->moveBuffer(this->array + (this->block_size * (this->current_blocks - 1)));

  // Do we need a new superblock for the sample?
  if(this->use_small_blocks)
  {
    this->current_samples++;
    if(this->current_samples > this->samples_in_superblock)
    {
      this->sample_blocks.push_back(this->samples);
      this->samples = new size_t[this->superblock_bytes / sizeof(size_t)];
      this->current_samples = 1;
    }
    this->samples[2 * this->current_samples - 2] = this->items - 1;
    this->samples[2 * this->current_samples - 1] = this->size - 1;
  }
}

void
VectorEncoder::setFirstBit(size_t value)
{
  this->samples[0] = 0;
  this->samples[1] = value;

  this->size = value + 1;
  this->items = 1;
  this->blocks = 1;

  this->current_blocks = 1;
  this->current_samples = 1;
}

}
