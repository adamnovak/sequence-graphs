#include "FMDPosition.hpp"

FMDPosition::FMDPosition(usint forward_start, usint reverse_start,
  usint end_offset): forward_start(forward_start), reverse_start(reverse_start),
  end_offset(end_offset)
{
}

FMDPosition::FMDPosition(): forward_start(0), reverse_start(0), end_offset(-1)
{
}

FMDPosition
FMDPosition::flip() const
{
  // Swap the two intervals of the bi-interval
  return FMDPosition(reverse_start, forward_start, end_offset);
}

bool FMDPosition::operator==(const FMDPosition& other) const
{
  // Compare all the fields.
  return
    forward_start == other.forward_start &&
    reverse_start == other.reverse_start &&
    end_offset == other.end_offset;
}

sint
FMDPosition::range(const RangeVector& ranges) const
{
  // Get an iterator for making rank queries.
  // TODO: Is it efficient to do this a lot?
  RangeVector::Iterator iter(ranges);

  // Look up the range that the forward starting position is in
  sint start_range = iter.rank(forward_start);
  
  // And the range the forward end is in
  sint end_range = iter.rank(forward_start + end_offset);
  
  if(start_range == end_range)
  {
    // Both ends of the interval are in the same range.
    return start_range;
  }
  else
  {
    // There is no range spanning our forward-strand interval.
    return -1;
  }
}

sint
FMDPosition::ranges(const RangeVector& ranges) const
{
  // Get an iterator for making rank queries.
  // TODO: Is it efficient to do this a lot?
  RangeVector::Iterator iter(ranges);

  // Look up the range that the starting position is in
  sint start_range = iter.rank(forward_start);
  
  // And the range the end is in
  sint end_range = iter.rank(forward_start + end_offset);
  
  // Return the number of ranges we intersect (1s hit plus 1)
  return(end_range - start_range + 1);
}

std::ostream& operator<< (std::ostream& o, FMDPosition const& position)
{
  // Report both the ranges that we represent.
  return o << position.forward_start << "-" << 
    (position.forward_start + position.end_offset) << "|" << 
    position.reverse_start << "-" << (position.reverse_start + 
    position.end_offset);
}
