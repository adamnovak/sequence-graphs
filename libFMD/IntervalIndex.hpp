#ifndef INTERVALINDEX_HPP
#define INTERVALINDEX_HPP

#include "GenericBitVector.hpp"
#include "Log.hpp"

#include <vector>
#include <algorithm>
#include <utility>
#include <memory>

/**
 * Represents a container for (start, length) intervals of size_ts, with an
 * Annotation associated with each interval.
 *
 * Supports finding the latest-starting interval that starts at or before a
 * position, and the earliest-ending interval that ends at or after a position.
 *
 * The Allocator type parameter is just to let us accept vectors of any
 * allocator in the constructor. TODO: Is there a way to get it to infer
 * properly from the allocator of whatever vector is passed to the constructor?
 */
template<typename Annotation, typename Allocator = 
    std::allocator<std::pair<std::pair<size_t, size_t>, Annotation>>>
class IntervalIndex {
public:
    // We define some typedefs like std::map does.

    /**
     * Type of the interval keys. Pair of start, length.
     */
    typedef std::pair<size_t, size_t> key_type; 
    
    /**
     * Type of the values we store for the keys.
     */
    typedef Annotation mapped_type;
    
    /**
     * Type of the conceptual records from intervals to annotations.
     */
    typedef std::pair<key_type, mapped_type> value_type;
    
    /**
     * Type of the const iterator over the contained mapping pairs.
     */
    typedef typename std::vector<value_type>::const_iterator const_iterator;
    
    
    /**
     * Create a new empty IntervalIndex.
     */
    IntervalIndex(): IntervalIndex(std::vector<value_type, Allocator>()) {
        // Nothing to do!
    }
    
    /**
     * Create a new interval index given the possibly unsorted vector of
     * intervals and their associated values.
     *
     * Takes O(n) time if intervals are already sorted by both start and end
     * coordinates, and O(n log n) time otherwise.
     */
    IntervalIndex(const std::vector<value_type, Allocator>& intervals): 
        records(intervals) {
    
        if(!std::is_sorted(records.begin(), records.end())) {
            // If not already sorted (check is O(n), sort is O(n log n))...
            // Sort records by start, length (and then value)
            std::sort(records.begin(), records.end());
        }
        
        // What's the past-the-end position for the last interval?
        size_t totalLength = 0;
        if(records.size() > 0) {
            // If we have a last interval, work that out. Grab the key and find
            // the after-the-end position.
            const key_type& lastKey = records[records.size() - 1].first;
            totalLength =  lastKey.first + lastKey.second;
        }
        
        // We're going to make a 1 in this bit vector at every place an interval
        // starts.
        startBits = new GenericBitVector();
        
        for(size_t i = 0; i < records.size(); i++) {
            // For each record
            
            if(i > 0 && records[i].first.first == records[i - 1].first.first) {
                // This interval starts at the same place as the other interval,
                // so we can ignore it in our index.
                continue;
            }
            
            // Note that this is an interval starting at this position.
            startBits->addBit(records[i].first.first);
            startRecords.push_back(i);
            Log::info() << "Start bit at " << records[i].first.first << std::endl;
        }
        
        // We're done marking start positions.
        startBits->finish(totalLength);
        
        // We're going to make a similar index of end positions facing the other
        // way, so  we need this list of end positions and interval numbers,
        // sorted by end.
        std::vector<std::pair<size_t, size_t>> ends;
        
        for(size_t i = 0; i < records.size(); i++) {
            // Make an entry with each interval's number under its end position.
            ends.push_back(std::make_pair(
                records[i].first.first + records[i].first.second - 1, i));
        }
        
        if(!std::is_sorted(ends.begin(), ends.end())) {
            // Sort the interval indices by end position, if they aren't in
            // order already.
            std::sort(ends.begin(), ends.end());
        }
        
        // We keep them in ascending order because it's not all that hard to get
        // the rank in either direction.
        
        // Make the bit vector for indexing the end positions.
        endBits = new GenericBitVector();
        
        for(size_t i = 0; i < ends.size(); i++) {
            // For each end, record number pair
            
            if(i > 0 && ends[i].first == ends[i - 1].first) {
                // This interval ends at the same place as the other interval,
                // so we can ignore it in our index.
                continue;
            }
            
            // Note that this points to an interval ending at this position.
            endBits->addBit(ends[i].first);
            endRecords.push_back(ends[i].second);
            Log::info() << "End bit at " << ends[i].first << std::endl;
        }
        
        // Finish off the bit vector with the total length of the region we care
        // about.
        endBits->finish(totalLength);
        
    }
    
    /**
     * Clean up an IntervalIndex.
     */
    ~IntervalIndex() {
        // Delete our dynamically-allocated memory.
        delete startBits;
        delete endBits;
    }

    /**
     * Copy an IntervalIndex.
     */
    IntervalIndex(const IntervalIndex& other): records(other.records),
        startBits(new GenericBitVector(*(other.startBits))),
        startRecords(other.startRecords), 
        endBits(new GenericBitVector(*(other.endBits))), 
        endRecords(other.endRecords) {
        
        // Nothing to do!
    
    }
    
    /**
     * Replace an IntervalIndex with a copy of another.
     */
    IntervalIndex& operator=(const IntervalIndex& other) {
        // Delete our dynamically-allocated memory.
        delete startBits;
        delete endBits;
        
        // Copy the data
        records = other.records;
        startBits = new GenericBitVector(*(other.startBits));
        startRecords = other.startRecords;
        endBits = new GenericBitVector(*(other.endBits)); 
        endRecords = other.endRecords;
    }
    
    /**
     * How many intervals are in this IntervalIndex?
     */
    size_t size() const {
        return records.size();
    }
    
    /**
     * Get the index'th interval in the index. Index must be less than size().
     * TODO: overload for key_type.
     */
    const value_type& operator[](size_t index) const {
        return records[index];
    }
    
    /**
     * Provide a begin iterator to iterate through the key, value pairs in this
     * IntervalIndex.
     */
    const_iterator begin() const {
        return records.begin();
    }
    
    /**
     * Provide an end iterator to iterate through the key, value pairs in this
     * IntervalIndex.
     */
    const_iterator end() const {
        return records.begin();
    }
    
    /**
     * Return true if an interval exists starting at or before the given
     * position, and false otherwise.
     */
    bool hasStartingBefore(size_t index) const {
        if(index >= startBits->getSize()) {
            // Going too far off the end.
            return hasStartingBefore(startBits->getSize() - 1);
        } else {
            return startBits->rank(index, false);
        }
    }
    
    /**
     * Get the latest starting interval that starts at or before the given
     * index, and its associated data value.
     */
    const value_type& getStartingBefore(size_t index) const {
    
        if(index >= startBits->getSize()) {
            // Going too far off the end.
            return getStartingBefore(startBits->getSize() - 1);
        }
    
        // How many positions where intervals start are before or at that
        // position?
        size_t rank = startBits->rank(index, false);
        
        if(rank == 0) {
            // No interval starts before or at the given position.
            throw std::runtime_error("No interval starting at or before " +
                std::to_string(index));
        }
        
        // We know we have a result, so get the record for that rank.
        return records[startRecords[rank - 1]];
    }
    
    /**
     * Return true if an interval exists ending at or before the given position,
     * and false otherwise.
     */
    bool hasEndingBefore(size_t index) const {
        if(index >= endBits->getSize()) {
            // Going too far off the end.
            return hasEndingBefore(endBits->getSize() - 1);
        } else {
            return endBits->rank(index, false);
        }
    }
    
    /**
     * Get the latest ending interval that ends at or before the given index,
     * and its associated data value.
     */
    const value_type& getEndingBefore(size_t index) const {
    
        if(index >= endBits->getSize()) {
            // Going too far off the end.
            return getEndingBefore(endBits->getSize() - 1);
        }
    
        // How many positions where intervals end are before or at that
        // position?
        size_t rank = endBits->rank(index, false);
        
        if(rank == 0) {
            // No interval ends before or at the given position.
            throw std::runtime_error("No interval ending at or before " +
                std::to_string(index));
        }
        
        // We know we have a result, so get the record for that rank.
        return records[endRecords[rank - 1]];
    }
    
    /**
     * Returns true if an interval exists ending at or after the given position,
     * and false otherwise.
     */
    bool hasEndingAfter(size_t index) const {
        // There is an interval ending at or after the given index if all of the
        // interval endpoints aren't already before the position.
        if(index >= endBits->getSize()) {
            // Nothing ends at or after the past-the-end position.
            return false;
        } else {
            
            return endBits->rank(index, true) - 1 < endRecords.size();
        }
    }
    
    /**
     * Get the earliest ending interval that ends at or after the given index,
     * and its associated data value.
     */
    const value_type& getEndingAfter(size_t index) const {
        // How many interval ending positions are before this index? If this is
        // 0, the soonest-ending interval ending here or later will be the
        // first-ending interval, and we count up from there.
        size_t rank = endBits->rank(index, true) - 1;
        
        // Go get and return that interval.
        return records[endRecords[rank]];
    }
    
    /**
     * Returns true if an interval exists starting at or after the given
     * position, and false otherwise.
     */
    bool hasStartingAfter(size_t index) const {
        // There is an interval starting at or after the given index if all of
        // the interval start points aren't already before the position.
        if(index >= startBits->getSize()) {
            // Nothing starts at or after the past-the-end position.
            return false;
        } else {
            
            return startBits->rank(index, true) - 1 < startRecords.size();
        }
    }
    
    /**
     * Get the earliest starting interval that starts at or after the given
     * index, and its associated data value.
     */
    const value_type& getStartingAfter(size_t index) const {
        // How many interval starting positions are before this index? If this
        // is 0, the soonest-starting interval starting here or later will be
        // the first-starting interval, and we count up from there.
        size_t rank = startBits->rank(index, true) - 1;
        
        // Go get and return that interval.
        return records[startRecords[rank]];
    }
     
     
        


private:
    /**
     * Holds all the intervals and their annotations, defining a size_t index
     * for each.
     */
    std::vector<value_type> records;
    
    /**
     * Holds a 1 at each position at which an interval starts.
     */
    GenericBitVector* startBits;
    
    /**
     * Holds the index of some interval that starts at a position, by the
     * position's bit rank in startBits.
     */
    std::vector<size_t> startRecords;
    
    /**
     * Holds a 1 at each position at shich an interval ends.
     */
    GenericBitVector* endBits;
    
    /**
     * Holds the index of some interval that ends at a position, by the
     * position's bit rank in endBits.
     */
    std::vector<size_t> endRecords;
    
};

#endif
