#include <iostream>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "FMDIndex.hpp"
#include "util.hpp"
#include "debug.hpp"

FMDIndex::FMDIndex(std::string basename, SuffixArray* fullSuffixArray): 
    names(), lengths(), cumulativeLengths(), bwt(basename + ".bwt"),
    suffixArray(basename + ".ssa"), fullSuffixArray(fullSuffixArray) {

    // We already loaded the index itself in the initializer. Go load the
    // length/order metadata.
    
    // Open the contig name/length file for reading.
    std::ifstream contigFile((basename + ".chrom.sizes").c_str());
    
    // Keep a cumulative length sum
    size_t lengthSum = 0;
    
    // Have a string to hold each line in turn.
    std::string line;
    while(std::getline(contigFile, line)) {
        // For each <contig>\t<length> pair...
        
        // Find the \t
        size_t separator = line.find('\t');
        
        // Split into name and length
        std::string contigName = line.substr(0, separator - 1);
        std::string contigLength = line.substr(separator + 1, line.size() - 1);
        
        // Parse the length
        long long int lengthNumber = atoll(contigLength.c_str());
        
        // Add it to the vector of names in number order.
        names.push_back(contigName);
        
        // And the vector of sizes in number order
        lengths.push_back(lengthNumber);
        
        // And the vector of cumulative lengths
        cumulativeLengths.push_back(lengthSum);
        lengthSum += lengthNumber;
    }
    // Close up the contig file. We read our contig metadata.
    contigFile.close();
    
}

FMDIndex::~FMDIndex() {
    if(fullSuffixArray != NULL) {
        // If we were holding a full SuffixArray, throw it out.
        delete fullSuffixArray;
    }
}

size_t FMDIndex::getContigNumber(TextPosition base) const {
    // What contig corresponds to that text? Contigs all have both strands.
    return base.getText() / 2;
}

bool FMDIndex::getStrand(TextPosition base) const {
    // What strand corresponds to that text? Strands are forward, then reverse.
    return base.getText() % 2 == 1;
}

size_t FMDIndex::getOffset(TextPosition base) const {
    // What base offset, 1-based, from the left, corresponds to this pair_type,
    // which may be on either strand.
    if(getStrand(base) == 0) {
        // We're on the forward strand. Make offset 1-based.
        return base.getOffset() + 1;
    } else {
        // We're on the reverse strand, so we measured from the end. Make it
        // 1-based.
        return getContigLength(getContigNumber(base)) - base.getOffset();
    }
}

std::string FMDIndex::getName(TextPosition base) const {

    // Unpack the coordinate parts.
    size_t contig = getContigNumber(base);
    size_t offset = getOffset(base);
    
    // Work out what to name the position.
    std::stringstream nameStream;
    // Leave off the strand.
    nameStream << "N" << contig << "B" << offset;
    return nameStream.str(); 
    
}

size_t FMDIndex::getBaseID(TextPosition base) const {
    // Get the cumulative total of bases by the start of the given contig
    size_t total = cumulativeLengths[getContigNumber(base)];
    
    // Add in the offset of this base from the start of its contig, convert back
    // to 0-based, and return.
    return total + getOffset(base) - 1;
}

size_t FMDIndex::getContigs() const {
    // How many contigs do we know about?
    return names.size();
}
    
const std::string& FMDIndex::getContigName(size_t index) const {
    // Get the name of that contig.
    return names[index];
}

size_t FMDIndex::getContigLength(size_t index) const {
    // Get the length of that contig.
    return lengths[index];
}

int64_t FMDIndex::getTotalLength() const {
    // Sum all the contig lengths and double (to make it be for both strands).
    // See <http://stackoverflow.com/a/3221813/402891>
    return std::accumulate(lengths.begin(), lengths.end(), 0) * 2;
}

int64_t FMDIndex::getBWTLength() const {
    return bwt.getBWLen();
}

FMDPosition FMDIndex::getCoveringPosition() const {
    // We want an FMDPosition that covers the entire BWT.
    
    // Just construct one. TODO: Will this extend properly since it goes to 0?
    return FMDPosition(0, 0, getBWTLength() - 1);
}
    
   
FMDPosition FMDIndex::getCharPosition(char c) const {
    // Starting a search with this character
    
    // See BWTAlgorithms::initInterval
    
    // Start the forward string with this character.
    int64_t forwardStart = bwt.getPC(c);
    
    // Start the reverse string with its complement.
    int64_t reverseStart = bwt.getPC(complement(c));
    
    // Get the offset to the end of the first interval (as well as the second).
    int64_t offset = bwt.getOcc(c, bwt.getBWLen() - 1) - 1;

    // Make the FMDPosition.
    return FMDPosition(forwardStart, reverseStart, offset);
    
}
     
   
FMDPosition FMDIndex::extend(FMDPosition range, char c, bool backward) const {
    // Extend the search with this character.
    
    // More or less directly implemented off of algorithms 2 and 3 in "Exploring
    // single-sample SNP and INDEL calling with whole-genome de novo assembly"
    // (Li, 2012). However, our character indices are one less, since we don't
    // allow search patterns to include the end-of-text symbol. We also use
    // alphabetical ordering instead of the paper's N-last ordering in the FM-
    // index, and consequently need to assign reverse ranges in alphabetical
    // order by reverse complement.
  
    if(!backward) {
        // We only really want to implement backwards search. Flip the interval,
        // do backwards search with the complement of the base, and then flip
        // back.
        return extend(range.flip(), complement(c), true).flip();        
    }

    if(c == '\0') {
        throw std::runtime_error("Can't extend with null byte!");
    }

    if(!isBase(c)) {
        std::string errorMessage = std::string("Character #");
        errorMessage.push_back(c);
        errorMessage += std::string(" is not a DNA base.");
        throw std::runtime_error(errorMessage);
    }

    DEBUG(std::cout << "Extending " << range << " backwards with " << c <<
        std::endl;)

    // We have an array of FMDPositions, one per base, that we will fill in by a
    // tiny dynamic programming.
    FMDPosition answers[NUM_BASES];

    for(size_t base = 0; base < NUM_BASES; base++) {
        // Go through the bases in arbitrary order.

        DEBUG(std::cout << "\tThinking about base " << base << "(" << 
            BASES[base] << ")" << std::endl;)

        // Count up the number of characters < this base.
        int64_t start = bwt.getPC(c);

        DEBUG(std::cout << "\t\tstart = " << start << std::endl;)

        // Get the rank among occurrences of the first instance of this base in
        // this slice.
        int64_t forwardStartRank = bwt.getOcc(BASES[base], 
            range.getForwardStart() - 1);
        
        // Get the same rank for the last instance. TODO: Is the -1 right here?
        int64_t forwardEndRank = bwt.getOcc(BASES[base], 
            range.getForwardStart() + range.getEndOffset()) - 1;

        // Fill in the forward-strand start position and range end offset for
        // this base's answer.
        answers[base].setForwardStart(start + forwardStartRank);
        answers[base].setEndOffset(forwardEndRank - forwardStartRank);

        DEBUG(std::cout << "\t\tWould go to: " << answers[base] << std::endl;)
    }

    // Since we don't keep an FMDPosition for the non-base end-of-text
    // character, we need to track its length separately in order for the DP
    // algorithm given in the paper to be implementable. We calculate
    // occurrences of the text end character (i.e. how much of the current range
    // is devoted to things where an endOfText comes next) implicitly: it's
    // whatever part of the length of the range is unaccounted-for by the other
    // characters. We need to use the length accessor because ranges with one
    // thing have the .end_offset set to 0.
    int64_t endOfTextLength = range.getLength();

    for(size_t base = 0; base < NUM_BASES; base++)
    {
        // Go through the bases in order and account for their lengths.
        endOfTextLength -= answers[base].getLength();
    }


    DEBUG(std::cout << "\tendOfTextLength = " << endOfTextLength << std::endl;)

    // The endOfText character is the very first character we need to account
    // for when subdividing the reverse range and picking which subdivision to
    // take.
    DEBUG(std::cout << "\tendOfText reverse_start would be " << 
        range.getReverseStart() << std::endl;)

    // Next, allocate the range for the base that comes first in alphabetical
    // order by reverse complement.
    answers[0].setReverseStart(range.getReverseStart() + endOfTextLength);
    DEBUG(std::cout << "\t" << BASES[0] << " reverse_start is " << 
    answers[0].getReverseStart() << std::endl;)

    for(size_t base = 1; base < NUM_BASES; base++)
    {
        // For each subsequent base in alphabetical order by reverse complement
        // (as stored in BASES), allocate it the next part of the reverse range.

        answers[base].setReverseStart(answers[base - 1].getReverseStart() + 
            answers[base - 1].getLength());
        DEBUG(std::cout << "\t" << BASES[base] << " reverse_start is " << 
        answers[base].getReverseStart() << std::endl;)
    }

    // Now all the per-base answers are filled in.

    for(size_t base = 0; base < NUM_BASES; base++)
    {
        // For each base in arbitrary order
        if(BASES[base] == c)
        {
            // This is the base we're actually supposed to be extending with. Return
            // its answer.
            
            DEBUG(std::cout << "Moving " << range << " to " << answers[base] << 
                " on " << BASES[base] << std::endl;)
            
            return answers[base];
        }
    }

    // If we get here, they gave us something not in BASES somehow, even though
    // we checked already.
    throw std::runtime_error("Unrecognized base");
}

FMDPosition FMDIndex::count(std::string pattern) const {
    if(pattern.size() == 0) {
        // We match everything! Say the whole range of the BWT.
        return getCoveringPosition();
    }
    

    // Start at the end and select the first character.
    FMDPosition position = getCharPosition(pattern[pattern.size() - 1]);
    
    for(int i = pattern.size() - 2; !position.isEmpty() && i >= 0; i--) {
        // Extend backwards with each character
        position = extend(position, pattern[i], true);
    }
    
    // We either ran out of matching locations or finished the pattern.
    return position;

}

TextPosition FMDIndex::locate(int64_t index) const {
    // Wrap up locate functionality.
    
    // We're going to fill in an SAElem: a composite thing where the high bits
    // encode the text number, and the low bits encode the offset.
    SAElem bitfield;

    if(fullSuffixArray != NULL) {
        // We can just look at the full suffix array cheat sheet.
        
        bitfield = fullSuffixArray->get(index);
        
    } else {
        // We need to use the sampled suffix array.
        
        // Run the libsuffixtools locate. 
        bitfield = suffixArray.calcSA(index, &bwt); 
        
    }
    
    // Unpack it and convert to our own format.
    return TextPosition(bitfield.getID(), bitfield.getPos());
}

char FMDIndex::display(int64_t index) const {
    // Just pull straight from the BWT string.
    return bwt.getChar(index);
}

char FMDIndex::displayFirst(int64_t index) const {
    // Our BWT supports this natively.
    return bwt.getF(index);
}

int64_t FMDIndex::getLF(int64_t index) const {
    // Find the character we're looking at
    char toFind = display(index);
    
    // Find the start of that character in the first column. It's just the
    // number of characters less than it, counting text stops.
    int64_t charBlockStart = bwt.getPC(toFind);
    
    // Find the rank of that instance of that character among instances of the
    // same character in the last column. Subtract 1 from occurrences since the
    // first copy should be rank 0.
    int64_t instanceRank = bwt.getOcc(toFind, index) - 1;
    
    // Add that to the start position to produce the LF mapping.
    return charBlockStart + instanceRank;
}

std::vector<Mapping> FMDIndex::map(const std::string& query, int start,
    int length) const {

    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }

    // We need a vector to return.
    std::vector<Mapping> mappings;

    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to track our FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.
    MapAttemptResult location;
    // Make sure the scratch position is empty so we re-start on the first base.
    // Other fields get overwritten.
    location.position = EMPTY_FMD_POSITION;

    for(size_t i = start; i < start + length; i++)
    {
        if(location.position.isEmpty())
        {
            INFO(std::cout << "Starting over by mapping position " << i <<
            std::endl;)
            // We do not currently have a non-empty FMDPosition to extend. Start
            // over by mapping this character by itself.
            location = this->mapPosition(query, i);
        } else {
            INFO(std::cout << "Extending with position " << i << std::endl;)
            // The last base either mapped successfully or failed due to multi-
            // mapping. Try to extend the FMDPosition we have to the right (not
            // backwards) with the next base.
            location.position = this->extend(location.position, query[i],
                false);
            location.characters++;
        }

        if(location.is_mapped && location.position.getLength() == 1) {
            // It mapped. We didn't do a re-start and fail, and there's exactly
            // one thing in our interval.

            // Take the first (only) thing in the bi-interval's forward strand
            // side.
            int64_t start = location.position.getForwardStart();

            // Locate it, and then report position as a (text, offset) pair.
            // This will give us the position of the first base in the pattern,
            // which lets us infer the position of the last base in the pattern.
            TextPosition textPosition = locate(start);

            INFO(std::cout << "Mapped " << location.characters << 
            " context to text " << textPosition.getText() << " position " << 
            textPosition.getOffset() << std::endl;)

            // Correct to the position of the last base in the pattern, by
            // offsetting by the length of the pattern that was used. A
            // 2-character pattern means we need to go 1 further right in the
            // string it maps to to find where its rightmost character maps.
            textPosition.setOffset(textPosition.getOffset() + 
                (location.characters - 1));

            // Add a Mapping for this mapped base.
            mappings.push_back(Mapping(textPosition));

            // We definitely have a non-empty FMDPosition to continue from

        } else {

            INFO(std::cout << "Failed (" << location.position.getLength() << 
            " options for " << location.characters << " context)." << std::endl;)

            if(location.is_mapped && location.position.isEmpty()) {
                // We extended right until we got no results. We need to try
                // this base again, in case we tried with a too-long left
                // context.

                INFO(std::cout << "Restarting from here..." << std::endl;)

                // Move the loop index back
                i--;

                // Since the FMDPosition is empty, on the next iteration we will
                // retry this base.

            } else {
                // It didn't map for some other reason:
                // - It was an initial mapping with too little left context to 
                //   be unique
                // - It was an initial mapping with a nonexistent left context
                // - It was an extension that was multimapped and still is

                // In none of these cases will re-starting from this base help
                // at all. If we just restarted here, we don't want to do it
                // again. If it was multimapped before, it had as much left
                // context as it could take without running out of string or
                // getting no results.

                // It didn't map. Add an empty/unmapped Mapping.
                mappings.push_back(Mapping());

                // Mark that the next iteration will be an extension (if we had
                // any results this iteration; if not it will just restart)
                location.is_mapped = true;

            }
        }
    }

    // We've gone through and attempted the whole string. Give back our answers.
    return mappings;

}

std::vector<int64_t> FMDIndex::map(const RangeVector& ranges,
    const std::string& query, int start, int length) const {
    
    // RIGHT-map to a range.

    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }

    // We need a vector to return.
    std::vector<int64_t> mappings;

    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to trackour FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.
    MapAttemptResult location;
    // Make sure the scratch position is empty so we re-start on the first base
    location.position = EMPTY_FMD_POSITION;

    for(int i = start + length - 1; i >= start; i--) {
        // Go from the end of our selected region to the beginning.

        DEBUG(std::cout << "On position " << i << " from " <<
            start + length - 1 << " to " << start << std::endl;)

        if(location.position.isEmpty()) {
            INFO(std::cout << "Starting over by mapping position " << i <<
            std::endl;)
            // We do not currently have a non-empty FMDPosition to extend. Start
            // over by mapping this character by itself.
            location = this->mapPosition(ranges, query, i);
        } else {
            INFO(std::cout << "Extending with position " << i << std::endl;)
            // The last base either mapped successfully or failed due to multi-
            // mapping. Try to extend the FMDPosition we have to the left
            // (backwards) with the next base.
            location.position = this->extend(location.position, query[i], true);
            location.characters++;
        }

        // What range index does our current left-side position (the one we just
        // moved) correspond to, if any?
        int64_t range = location.position.range(ranges);

        if(location.is_mapped && !location.position.isEmpty() && range != -1) {
            // It mapped. We didn't do a re-start and fail, and our interval is
            // nonempty and subsumed by a range.

            INFO(std::cout << "Mapped " << location.characters << 
            " context to range #" << range << " in range vector." << std::endl;)

            // Remember that this base mapped to this range
            mappings.push_back(range);

            // We definitely have a non-empty FMDPosition to continue from

        } else {

            INFO(std::cout << "Failed (" << location.position.ranges(ranges) << 
                " options for " << location.characters << " context)." << 
                std::endl;)

            if(location.is_mapped && location.position.isEmpty()) {
                // We extended right until we got no results. We need to try
                // this base again, in case we tried with a too-long left
                // context.

                INFO(std::cout << "Restarting from here..." << std::endl;)

                // Move the loop index towards the end we started from (right)
                i++;

                // Since the FMDPosition is empty, on the next iteration we will
                // retry this base.

            } else {
                // It didn't map for some other reason:
                // - It was an initial mapping with too little left context to 
                //   be unique to a range.
                // - It was an initial mapping with a nonexistent left context
                // - It was an extension that was multimapped and still is

                // In none of these cases will re-starting from this base help
                // at all. If we just restarted here, we don't want to do it
                // again. If it was multimapped before, it had as much left
                // context as it could take without running out of string or
                // getting no results.

                // It didn't map. Say it corresponds to no range.
                mappings.push_back(-1);

                // Mark that the next iteration will be an extension (if we had
                // any results this iteration; if not it will just restart)
                location.is_mapped = true;

            }
        }
    }

    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());

    // Give back our answers.
    return mappings;
    
    
}

FMDIndex::iterator FMDIndex::begin(size_t depth, bool reportDeadEnds) const {
    // Make a new suffix tree iterator that automatically searches out the first
    // suffix of the right length.
    return FMDIndex::iterator(*this, depth, false, reportDeadEnds);
}
     
FMDIndex::iterator FMDIndex::end(size_t depth, bool reportDeadEnds) const {
    // Make a new suffix tree iterator that is just a 1-past-the-end sentinel.
    return FMDIndex::iterator(*this, depth, true, reportDeadEnds);
}

MapAttemptResult FMDIndex::mapPosition(const std::string& pattern,
    size_t index) const {

    DEBUG(std::cout << "Mapping " << index << " in " << pattern << std::endl;)
  
    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MapAttemptResult result;

    // Do a backward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.position = this->getCharPosition(pattern[index]);
    result.characters = 1;
    if(result.position.isEmpty()) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if(result.position.getLength() == 1) {
        // We've already mapped.
        result.is_mapped = true;
        return result;
    }

    if(index == 0) {
        // The rest of the function deals with characters to the left of the one
        // we start at. If we start at position 0 there can be none.
        return result;
    }

    DEBUG(std::cout << "Starting with " << result.position << std::endl;)

    do {
        // Now consider the next character to the left.
        index--;

        // Grab the character to extend with.
        char character = pattern[index];

        DEBUG(std::cout << "Index " << index << " in " << pattern << " is " << 
            character << "(" << character << ")" << std::endl;)

        // Backwards extend with subsequent characters.
        FMDPosition next_position = this->extend(result.position, character,
            true);

        DEBUG(std::cout << "Now at " << next_position << " after " << 
            pattern[index] << std::endl;)
        if(next_position.isEmpty()) {
            // The next place we would go is empty, so return the result holding
            // the last position.
            return result;
        } else if(next_position.getLength() == 1) {
            // We have successfully mapped to exactly one place. Update our
            // result to reflect the additional extension and our success, and
            // return it.
            result.position = next_position;
            result.characters++;
            result.is_mapped = true;
            return result;      
        }

        // Otherwise, we still map to a plurality of places. Record the
        // extension and loop again.
        result.position = next_position;
        result.characters++;

    } while(index > 0);
    // Continue considering characters to the left until we hit the start of the
    // string.

    // If we get here, we ran out of upstream context and still map to multiple
    // places. Just give our multi-mapping FMDPosition and unmapped result.
    return result;
}

MapAttemptResult FMDIndex::mapPosition(const RangeVector& ranges, 
    const std::string& pattern, size_t index) const {
    
    
    // We're going to right-map so ranges match up with the things we can map to
    // (downstream contexts)

    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MapAttemptResult result;

    // Do a forward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.position = this->getCharPosition(pattern[index]);
    result.characters = 1;
    if(result.position.isEmpty()) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if (result.position.range(ranges) != -1) {
        // We've already mapped.
        result.is_mapped = true;
        return result;
    }

    DEBUG(std::cout << "Starting with " << result.position << std::endl;)

    for(index++; index < pattern.size(); index++) {
        // Forwards extend with subsequent characters.
        FMDPosition next_position = this->extend(result.position,
            pattern[index], false);

        DEBUG(std::cout << "Now at " << next_position << " after " << 
            pattern[index] << std::endl;)
        if(next_position.isEmpty()) {
            // The next place we would go is empty, so return the result holding
            // the last position.
            return result;
        }

        if(next_position.range(ranges) != -1) {
            // We have successfully mapped to exactly one range. Update our
            // result to reflect the additional extension and our success, and
            // return it.
            result.position = next_position;
            result.characters++;
            result.is_mapped = true;
            return result;      
        }

        // Otherwise, we still map to a plurality of ranges. Record the
        // extension and loop again.
        result.position = next_position;
        result.characters++;
    }

    // If we get here, we ran out of downstream context and still map to
    // multiple ranges. Just give our multi-mapping FMDPosition and unmapped
    // result.
    return result;


}
















