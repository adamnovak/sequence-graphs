#include <cstdlib>
#include <cstring>

#include "BitBuffer.hpp"

namespace CSA {

//--------------------------------------------------------------------------

ReadBuffer::ReadBuffer(std::ifstream& file, size_t words):
size(words),
item_bits(1),
items(0),
free_buffer(true) {
    size_t* buffer = new size_t[this->size];
    memset(buffer, 0, this->size * sizeof (size_t));
    file.read((char*) buffer, this->size * sizeof (size_t));
    this->data = buffer;
    this->reset();
}

ReadBuffer::ReadBuffer(FILE* file, size_t words):
size(words),
item_bits(1),
items(0),
free_buffer(true) {
    size_t* buffer = new size_t[this->size];
    memset(buffer, 0, this->size * sizeof (size_t));
    if(file != 0) {
        if(!std::fread(buffer, this->size * sizeof (size_t), 1, file)) {
            return;
        }
    }
    this->data = buffer;
    this->reset();
}

ReadBuffer::ReadBuffer(std::ifstream& file, size_t _items, size_t item_size):
item_bits(item_size),
items(_items),
free_buffer(true) {
    this->size = bitsToWords(this->items * this->item_bits);
    size_t* buffer = new size_t[this->size];
    memset(buffer, 0, this->size * sizeof (size_t));
    file.read((char*) buffer, this->size * sizeof (size_t));
    this->data = buffer;
    this->reset();
}

ReadBuffer::ReadBuffer(FILE* file, size_t _items, size_t item_size):
item_bits(item_size),
items(_items),
free_buffer(true) {
    this->size = bitsToWords(this->items * this->item_bits);
    size_t* buffer = new size_t[this->size];
    memset(buffer, 0, this->size * sizeof (size_t));
    if(file != 0) {
        if(!std::fread(buffer, this->size * sizeof (size_t), 1, file)) {
            return;
        }
    }
    this->data = buffer;
    this->reset();
}

ReadBuffer::ReadBuffer(const size_t* buffer, size_t words):
size(words),
item_bits(1),
items(0),
free_buffer(false) {
    this->data = buffer;
    this->reset();
}

ReadBuffer::ReadBuffer(const size_t* buffer, size_t _items, size_t item_size):
item_bits(item_size),
items(_items),
free_buffer(false) {
    this->size = bitsToWords(this->items * this->item_bits);
    this->data = buffer;
    this->reset();
}

ReadBuffer::ReadBuffer(const ReadBuffer& original):
data(original.data),
size(original.size),
item_bits(original.item_bits),
items(original.items),
free_buffer(false) {
    this->reset();
}

ReadBuffer::~ReadBuffer() {
    if(this->free_buffer) {
        delete[] this->data;
    }
}

//--------------------------------------------------------------------------

void
ReadBuffer::claimData() {
    this->free_buffer = true;
}

void
ReadBuffer::writeTo(std::ofstream& file) const {
    file.write((const char*) this->data, this->size * sizeof (size_t));
}

void
ReadBuffer::writeTo(FILE* file) const {
    if(file == 0) {
        return;
    }
    std::fwrite(this->data, this->size * sizeof (size_t), 1, file);
}

void
ReadBuffer::moveBuffer(const size_t* buffer) {
    if(this->free_buffer) {
        delete[] this->data;
    }
    this->free_buffer = false;

    this->data = buffer;
    this->reset();
}

size_t
ReadBuffer::reportSize() const {
    size_t bytes = sizeof (*this);
    if(this->free_buffer) {
        bytes += this->size * sizeof (size_t);
    }
    return bytes;
}

//--------------------------------------------------------------------------

WriteBuffer::WriteBuffer(size_t words):
size(words),
item_bits(1),
items(0),
free_buffer(true) {
    this->data = new size_t[words];
    memset(this->data, 0, this->size * sizeof (size_t));
    this->reset();
}

WriteBuffer::WriteBuffer(size_t _items, size_t item_size):
item_bits(item_size),
items(_items),
free_buffer(true) {
    this->size = bitsToWords(this->items * this->item_bits);
    this->data = new size_t[this->size];
    memset(this->data, 0, this->size * sizeof (size_t));
    this->reset();
}

WriteBuffer::WriteBuffer(size_t* buffer, size_t words):
size(words),
item_bits(1),
items(0),
free_buffer(false) {
    this->data = buffer;
    this->reset();
}

WriteBuffer::WriteBuffer(size_t* buffer, size_t _items, size_t item_size):
item_bits(item_size),
items(_items),
free_buffer(false) {
    this->size = bitsToWords(this->items * this->item_bits);
    this->data = buffer;
    this->reset();
}

WriteBuffer::~WriteBuffer() {
    if(this->free_buffer) {
        delete[] this->data;
    }
}

//--------------------------------------------------------------------------

ReadBuffer*
WriteBuffer::getReadBuffer() {
    ReadBuffer* buffer;
    if(this->items > 0) {
        buffer = new ReadBuffer(this->data, this->items, this->item_bits);
    } else {
        buffer = new ReadBuffer(this->data, this->size);
    }

    if(this->free_buffer) {
        buffer->claimData();
        this->free_buffer = false;
    }

    return buffer;
}

void
WriteBuffer::writeTo(std::ofstream& file) const {
    file.write((char*) this->data, this->size * sizeof (size_t));
}

void
WriteBuffer::writeTo(FILE* file) const {
    if(file == 0) {
        return;
    }
    std::fwrite(this->data, this->size * sizeof (size_t), 1, file);
}

void
WriteBuffer::moveBuffer(size_t* buffer) {
    if(this->free_buffer) {
        delete[] this->data;
    }
    this->free_buffer = false;

    this->data = buffer;
    this->reset();
}

size_t
WriteBuffer::reportSize() const {
    size_t bytes = sizeof (*this);
    if(this->free_buffer) {
        bytes += this->size * sizeof (size_t);
    }
    return bytes;
}

}
//--------------------------------------------------------------------------
