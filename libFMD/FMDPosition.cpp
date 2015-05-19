#include <stdexcept>

#include "FMDPosition.hpp"
#include "FMDIndexView.hpp"

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

bool FMDPosition::isEmpty(const FMDIndexView& view) const {
    // Get our forward interval and see if it's empty.
    return view.isEmpty(getForwardStart(), getLength());
}

bool FMDPosition::isUnique(const FMDIndexView& view) const {
    // Get our forward interval and see if it's unique.
    return view.isUnique(getForwardStart(), getLength());
}

bool FMDPosition::isAmbiguous(const FMDIndexView& view) const {
    // Get our forward interval and see if it's ambiguous.
    return view.isAmbiguous(getForwardStart(), getLength());
}

TextPosition FMDPosition::getTextPosition(const FMDIndexView& view) const {
    // Get our text position from our interval
    return view.getTextPosition(getForwardStart(), getLength());
}

size_t FMDPosition::getApproximateNumberOfRanges(
    const FMDIndexView& view) const {
    
    // Get our range count from our interval
    return view.getApproximateNumberOfRanges(getForwardStart(), getLength());
}

size_t FMDPosition::getApproximateNumberOfNewRanges(const FMDIndexView& view,
    const FMDPosition& old) {
    
    // Get our range count from our interval and the old position's interval.
    return view.getApproximateNumberOfNewRanges(old.getForwardStart(),
        old.getLength(), getForwardStart(), getLength());
}
    
std::set<TextPosition> FMDPosition::getTextPositions(
    const FMDIndexView& view) const {
    
    // Get our text positions from our interval
    return view.getTextPositions(getForwardStart(), getLength());
}

std::set<TextPosition> FMDPosition::getNewTextPositions(
    const FMDIndexView& view, const FMDPosition& old) const {
    
    // Get our new text positions from our interval and the old position's
    // interval.
    return view.getNewTextPositions(old.getForwardStart(),
        old.getLength(), getForwardStart(), getLength());
}
