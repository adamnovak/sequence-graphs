#include <stdexcept>

#include "FMDPosition.hpp"

FMDPosition::FMDPosition(int64_t forward_start, int64_t reverse_start,
  int64_t end_offset): forward_start(forward_start), 
  reverse_start(reverse_start), end_offset(end_offset) {
}

FMDPosition::FMDPosition(): forward_start(0), reverse_start(0), end_offset(-1) {
}

FMDPosition FMDPosition::flip() const {
    // Swap the two intervals of the bi-interval
    return FMDPosition(reverse_start, forward_start, end_offset);
}

bool FMDPosition::operator==(const FMDPosition& other) const {
    // Compare all the fields.
    return
        forward_start == other.forward_start &&
        reverse_start == other.reverse_start &&
        end_offset == other.end_offset;
}

bool FMDPosition::operator<(const FMDPosition& other) const {
    // Compare all the fields like normal didgit places.
    // TODO: Is there a way to generate this code?
    return
        forward_start < other.forward_start ||
        (forward_start == other.forward_start &&
        (reverse_start < other.reverse_start ||
        (reverse_start == other.reverse_start &&
        end_offset <= other.end_offset)));
}

std::ostream& operator<< (std::ostream& o, FMDPosition const& position) {
    // Report both the ranges that we represent.
    return o << position.forward_start << "-" << 
        (position.forward_start + position.end_offset) << "|" << 
        position.reverse_start << "-" << (position.reverse_start + 
        position.end_offset);
}
