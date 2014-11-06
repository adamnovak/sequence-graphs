#include <iostream>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <Util.h>

#include "FMDIndex.hpp"
#include "util.hpp"
#include "Log.hpp"


FMDIndex::FMDIndex(std::string basename, SuffixArray* fullSuffixArray, 
    MarkovModel* markovModel): names(), starts(), lengths(), 
    cumulativeLengths(), genomeAssignments(), endIndices(), genomeRanges(), 
    genomeMasks(), bwt(basename + ".bwt"), suffixArray(basename + ".ssa"), 
    fullSuffixArray(fullSuffixArray), markovModel(markovModel),
    lcpArray(basename + ".lcp") {
    
    // TODO: Too many initializers

    Log::info() << "Loading " << basename << std::endl;

    // We already loaded the index itself in the initializer. Go load the
    // length/order metadata.
    
    // Open the contig name/start/length/genome file for reading.
    std::ifstream contigFile((basename + ".contigs").c_str());
    
    // Keep a cumulative length sum
    size_t lengthSum = 0;
    
    // Have a string to hold each line in turn.
    std::string line;
    while(std::getline(contigFile, line)) {
        // For each <contig>\t<start>\t<length>\t<genome> line...
        
        // Make a stringstream we can read out of.
        std::stringstream lineData(line);
        
        // Read in the name of the contig
        std::string contigName;
        lineData >> contigName;
        
        // Add it to the vector of names in number order.
        names.push_back(contigName);
        
        // Read in the contig start position on its scaffold
        size_t startNumber;
        lineData >> startNumber;
        
        // Add it to the vector of starts in number order.
        starts.push_back(startNumber);
        
        // Read in the contig's length
        size_t lengthNumber;
        lineData >> lengthNumber;
        
        // Add it to the vector of sizes in number order
        lengths.push_back(lengthNumber);
        
        // And to the vector of cumulative lengths
        cumulativeLengths.push_back(lengthSum);
        lengthSum += lengthNumber;
        
        // Read in the number of the genome that the contig belongs to.
        size_t genomeNumber;
        lineData >> genomeNumber;
        
        // Add it to the vector of genome assignments in number order
        genomeAssignments.push_back(genomeNumber);
    }
        
        
    // Close up the contig file. We read our contig metadata.
    contigFile.close();
    
    // Now read the genome bit masks.
    
    // What file are they in? Make sure to hold onto it while we construct the
    // stream with its c_str pointer.
    std::string genomeMaskFile = basename + ".msk";
    
    // Open the file where they live.
    std::ifstream genomeMaskStream(genomeMaskFile.c_str(), std::ios::binary);
    
    while(genomeMaskStream.peek() != EOF && !genomeMaskStream.eof()) {
        // As long as there is data left to read
        
        // Read a new GenericBitVector from the stream and put it in our list.
        genomeMasks.push_back(new GenericBitVector(genomeMaskStream));
        
        // This lets us autodetect how many genomes there are.
    }
    
    // Now invert the contig-to-genome index to make the genome-to-contig-range
    // index.
    
    // How many genomes are there?
    size_t numGenomes = genomeMasks.size();
    
    // Make the genome range vector big enough. Fill it with empty ranges for
    // genomes that somehow have no contigs.
    genomeRanges.resize(numGenomes, std::make_pair(0, 0));
    
    // Start out the first range.
    std::pair<size_t, size_t> currentRange = std::make_pair(0, 0);
    
    // Which genome are we on?
    size_t currentGenome;
    
    for(std::vector<size_t>::iterator i = genomeAssignments.begin(); 
        i != genomeAssignments.end(); ++i) {
    
        if(i == genomeAssignments.begin()) {
            // We're starting. Grab this genome as our current genome (since the
            // file may not start with genome 0).
            // TODO: Maybe just require that?
            currentGenome = *i;
        }
        
        if(*i == currentGenome) {
            // This is the same genome as the last one. Extend the range.
            currentRange.second++;
        } else {
            // This is now a different genome. Save the range.
            genomeRanges[currentGenome] = currentRange;
            
            // Make a new range that comes after it, with length 1, since we
            // include this contig.
            currentRange = std::make_pair(currentRange.second,
                currentRange.second + 1);
            
            // Remember what new genome we're looking at.
            currentGenome = *i;
            
            if(currentGenome >= numGenomes) {
                // Complain if we have a genome number higher than the number of
                // masks we loaded.
                throw std::runtime_error(
                    "Got a contig for a genome with no mask!");
            }
        }
    
    }
    
    // Save the last range
    genomeRanges[currentGenome] = currentRange;
    
    // First make sure the vector is big enough for them.
    endIndices.resize(getNumberOfContigs());
    
    for(int64_t i = 0; i < getNumberOfContigs() * 2; i++) {
        // The first #-of-texts rows in the BWT table have a '$' in the F
        // column, so the L column (what our BWT string actually is) will have
        // the last real character in some text.
        
        // Locate it to a text and offset
        TextPosition position = locate(i);
        
        if(position.getText() % 2 == 0) {
            // This is a forward strand. Save the index of the last real
            // character in the forward strand of the contig.
            endIndices[position.getText() / 2] = i;
        }
        
    }
    
    Log::info() << "Loaded " << names.size() << " contigs in " << numGenomes <<
        " genomes" << std::endl;
}

FMDIndex::~FMDIndex() {
    if(fullSuffixArray != NULL) {
        // If we were holding a full SuffixArray, throw it out.
        delete fullSuffixArray;
    }
    
    if(markovModel != NULL) {
        // If we were holding a markovModel, throw it out.
        delete markovModel;
    }
    
    for(std::vector<GenericBitVector*>::iterator i = genomeMasks.begin(); 
        i != genomeMasks.end(); ++i) {
        
        // Also delete all the genome masks we loaded.
        delete (*i);
    }
}

void FMDIndex::setMarkovModel(MarkovModel* model) {
    if(markovModel != NULL) {
        // Get rid of any Markov model we already had
        delete markovModel;
    }
    
    // Take ownership of this one
    markovModel = model;
}

size_t FMDIndex::getContigNumber(TextPosition base) const {
    // What contig corresponds to that text? Contigs all have both strands.
    return base.getContigNumber();
}

bool FMDIndex::getStrand(TextPosition base) const {
    // What strand corresponds to that text?
    return base.getStrand();
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

size_t FMDIndex::getNumberOfContigs() const {
    // How many contigs do we know about?
    return names.size();
}
    
const std::string& FMDIndex::getContigName(size_t index) const {
    // Get the name of that contig.
    return names[index];
}

size_t FMDIndex::getContigStart(size_t index) const {
    // Get the start of that contig.
    return starts[index];
}

size_t FMDIndex::getContigLength(size_t index) const {
    // Get the length of that contig.
    return lengths[index];
}

size_t FMDIndex::getContigGenome(size_t index) const {
    // Get the genome that that contig belongs to.
    return genomeAssignments[index];
}
    
size_t FMDIndex::getNumberOfGenomes() const {
    // Get the number of genomes that are in the index.
    // TODO: Several things are this length. There should only be one.
    return genomeMasks.size();
}
    
std::pair<size_t, size_t> FMDIndex::getGenomeContigs(size_t genome) const {
    // Get the range of contigs belonging to the given genome.
    return genomeRanges[genome];
}

size_t FMDIndex::getGenomeLength(size_t genome) {
    // What's the last contig in the range for the genome?
    size_t lastContig = getGenomeContigs(genome).second - 1;
    
    // Where does it end?
    return getContigStart(lastContig) + getContigLength(lastContig);
}

bool FMDIndex::isInGenome(int64_t bwtIndex, size_t genome) const {
    return genomeMasks[genome]->isSet(bwtIndex);
}

const GenericBitVector& FMDIndex::getGenomeMask(size_t genome) const {
    return *genomeMasks[genome];
}

TextPosition FMDIndex::getTextPosition(
    std::pair<std::pair<size_t, size_t>, bool> base) const {
    
    // Convert contig and face to text.
    size_t text = base.first.first * 2 + base.second;
    
    // Work out the offset.
    size_t offset;
    if(base.second) {
        // Convert to an offset from the start of the reverse strand, converting
        // to 0-based in the process.
        offset = getContigLength(base.first.first) - base.first.second;
    } else {
        // Just make 0-based.
        offset = base.first.second - 1;
    }
    
    // Return the conversion result.
    return TextPosition(text, offset);
        
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
     
void FMDIndex::extendFast(FMDPosition& range, char c, bool backward) const {
    // Extend the search with this character in an optimized way. We work on our
    // argument so we don't need to bother with copies.
    
    // Skip any sort of argument validation.
    
    if(!backward) {
        // Flip the arguments around so we can work on the reverse strand.
        c = complement(c);
        range.flipInPlace();
    }
    
    // Read occurrences of everything from the BWT
    
    // What rank among occurrences is the first instance of every character in
    // the BWT range?
    AlphaCount64 startRanks = bwt.getFullOcc(range.getForwardStart() - 1);
    
    // And the last? If endOffset() is 0, this will be 1 character later than
    // the call for startRanks, which is what we want.
    AlphaCount64 endRanks = bwt.getFullOcc(range.getForwardStart() + 
        range.getEndOffset());
        
    // Get the number of suffixes that had '$' (end of text) next. TODO: should
    // this be '\0' instead?
    range.setReverseStart(range.getReverseStart() + 
        (endRanks.get('$') - startRanks.get('$')));
        
    for(size_t base = 0; base < NUM_BASES; base++) {
        // For each base in alphabetical order by reverse complement
        
        // Work out the length of the interval this base gets.
        size_t intervalLength = endRanks.get(BASES[base]) - 
            startRanks.get(BASES[base]);
        
        if(BASES[base] == c) {
            // This is the base we're looking for. Finish up and break out of
            // the loop.
            
            // Range reverse start is already set.
            
            // Set the range forward start.
            range.setForwardStart(bwt.getPC(c) + startRanks.get(c));
            
            // Set the range length.
            range.setEndOffset((int64_t)intervalLength - 1);
            
            // Now we've put together the range.
            break;
            
        } else {
            // This is not the base we're looking for. Budge the reverse strand
            // interval over by the length of the interval, to account for the
            // bit this base took up.
            range.setReverseStart(range.getReverseStart() + intervalLength);
        }
    }
    
    if(!backward) {
        // Flip the result since we were working on the opposite strand.
        range.flipInPlace();
    }
    
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

    Log::trace() << "Extending " << range << " backwards with " << c <<
        std::endl;

    // We have an array of FMDPositions, one per base, that we will fill in by a
    // tiny dynamic programming.
    FMDPosition answers[NUM_BASES];

    for(size_t base = 0; base < NUM_BASES; base++) {
        // Go through the bases in arbitrary order.

        Log::trace() << "\tThinking about base " << base << "(" << 
            BASES[base] << ")" << std::endl;

        // Count up the number of characters < this base.
        int64_t start = bwt.getPC(c);

        Log::trace() << "\t\tstart = " << start << std::endl;

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

        Log::trace() << "\t\tWould go to: " << answers[base] << std::endl;
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


    Log::trace() << "\tendOfTextLength = " << endOfTextLength << std::endl;

    // The endOfText character is the very first character we need to account
    // for when subdividing the reverse range and picking which subdivision to
    // take.
    Log::trace() << "\tendOfText reverse_start would be " << 
        range.getReverseStart() << std::endl;

    // Next, allocate the range for the base that comes first in alphabetical
    // order by reverse complement.
    answers[0].setReverseStart(range.getReverseStart() + endOfTextLength);
    Log::trace() << "\t" << BASES[0] << " reverse_start is " << 
        answers[0].getReverseStart() << std::endl;

    for(size_t base = 1; base < NUM_BASES; base++)
    {
        // For each subsequent base in alphabetical order by reverse complement
        // (as stored in BASES), allocate it the next part of the reverse range.

        answers[base].setReverseStart(answers[base - 1].getReverseStart() + 
            answers[base - 1].getLength());
        Log::trace() << "\t" << BASES[base] << " reverse_start is " << 
        answers[base].getReverseStart() << std::endl;
    }

    // Now all the per-base answers are filled in.

    for(size_t base = 0; base < NUM_BASES; base++)
    {
        // For each base in arbitrary order
        if(BASES[base] == c)
        {
            // This is the base we're actually supposed to be extending with. Return
            // its answer.
            
            Log::trace() << "Moving " << range << " to " << answers[base] << 
                " on " << BASES[base] << std::endl;
            
            return answers[base];
        }
    }

    // If we get here, they gave us something not in BASES somehow, even though
    // we checked already.
    throw std::runtime_error("Unrecognized base");
}

void FMDIndex::extendLeftOnly(FMDPosition& range, char c) const {

    // Extend the search with this character in an optimized way. We work on our
    // argument so we don't need to bother with copies.
    
    // Skip any sort of argument validation.
    
    // We always do backward search, updating the forward interval (i.e.
    // reference taken forward)
    
    // Count up the number of characters < this base.
    int64_t start = bwt.getPC(c);
    
    // Get the rank among occurrences of the first instance of this base in
    // this slice.
    int64_t forwardStartRank = bwt.getOcc(c, range.getForwardStart() - 1);
    
    // Get the same rank for the last instance.
    int64_t forwardEndRank = bwt.getOcc(c, range.getForwardStart() + 
        range.getEndOffset()) - 1;

    // Fill in the forward-strand start position and range end offset for
    // the answer.
    range.setForwardStart(start + forwardStartRank);
    range.setEndOffset(forwardEndRank - forwardStartRank);
    
    // Leave the reverse interval alone

}

void FMDIndex::retractRightOnly(FMDPosition& range, 
    size_t newPatternLength) const {
    
    // Retract on the right in place. Can only extend left afterwards, since we
    // use the forward range and ignore the reverse range.
    
    // Get the bounds of our range in the LCP: Where our range starts, and 1
    // after it ends. (So 0 and 1 for a 1-element thing at 0).
    // The end offset is 0 for a 1-element range, so we fix it up.
    size_t rangeStart = range.getForwardStart();
    size_t rangeEnd = range.getForwardStart() + range.getEndOffset() + 1;
    
    Log::trace() << "Retracting from [" << rangeStart << ", " << rangeEnd << 
        ")" << std::endl;
    
    // rangeEnd may be actually past the end of the LCP array now. That is OK,
    // we just get 0 and an LCP NSV of the same position in that case.
        
    // Get the LCP value at each end
    size_t startLCP = getLCP(rangeStart);
    // Don't try looking off the end. Fill in an imaginary 0 to bound the root.
    size_t endLCP = rangeEnd < getBWTLength() ? getLCP(rangeEnd) : 0;
    
    // Figure out which end has the greater value, and what that value is, and
    // where in the LCP that value is. Default to using the start in ties,
    // because the start will always be at a real location, while the end can be
    // past the end of the LCP array.
    bool useStart = startLCP >= endLCP;
    size_t lcp = useStart ? startLCP : endLCP;
    size_t lcpIndex = useStart ? rangeStart : rangeEnd;
    
    // Now lcpIndex is guaranteed to be a real index in the LCP array, not off
    // the end.
    
    Log::trace() << "Parent node string depth: " << lcp << " at " << lcpIndex <<
        std::endl;
    
    // The larger LCP value cuts down to the string depth of the parent.
    
    if(lcp < newPatternLength) {
        // No reason to go anywhere if we don't get any more results. Only go up
        // to the parent node when the new pattern length will be as shoet as
        // its string depth (or shorter).
        return;
    } else {
        // The parent node string depth meets or excedes our new pattern length.
        // We need to be at that parent node or higher.
        
            
        // We'll update rangeStart and rangeEnd to the LCP array indices that
        // cut out the parent.
        rangeStart = getLCPPSV(lcpIndex);
        // Note that rangeEnd can be off the end of the LCP array now, if we
        // have moved up to the root node.
        rangeEnd = getLCPNSV(lcpIndex);
        
        // Update the Range object to describe this range, converting end
        // indices around. The range will never be empty so we don't have to
        // worry about negative end offsets.
        range.setForwardStart(rangeStart);
        range.setEndOffset(rangeEnd - rangeStart - 1);
        
        if(lcp > newPatternLength) {
            // We need to retract more to meet that depth target.
            retractRightOnly(range, newPatternLength);
        }
    }
}

size_t FMDIndex::retractRightOnly(FMDPosition& range) const {
    
    // Retract on the right in place. Can only extend left afterwards, since we
    // use the forward range and ignore the reverse range.
    
    // Goes all the way to the parent suffix tree node.
    
    // Get the bounds of our range in the LCP: Where our range starts, and 1
    // after it ends. (So 0 and 1 for a 1-element thing at 0).
    // The end offset is 0 for a 1-element range, so we fix it up.
    size_t rangeStart = range.getForwardStart();
    size_t rangeEnd = range.getForwardStart() + range.getEndOffset() + 1;
    
    Log::trace() << "Retracting from [" << rangeStart << ", " << rangeEnd << 
        ")" << std::endl;
    
    // rangeEnd may be actually past the end of the LCP array now. That is OK,
    // we just get 0 and an LCP NSV of the same position in that case.
        
    // Get the LCP value at each end
    size_t startLCP = getLCP(rangeStart);
    // Don't try looking off the end. Fill in an imaginary 0 to bound the root.
    size_t endLCP = rangeEnd < getBWTLength() ? getLCP(rangeEnd) : 0;
    
    // Figure out which end has the greater value, and what that value is, and
    // where in the LCP that value is. Default to using the start in ties,
    // because the start will always be at a real location, while the end can be
    // past the end of the LCP array.
    bool useStart = startLCP >= endLCP;
    size_t lcp = useStart ? startLCP : endLCP;
    size_t lcpIndex = useStart ? rangeStart : rangeEnd;
    
    // Now lcpIndex is guaranteed to be a real index in the LCP array, not off
    // the end.
    
    Log::trace() << "Parent node string depth: " << lcp << " at " << lcpIndex <<
        std::endl;
    
    // The larger LCP value cuts down to the string depth of the parent.
    
    // Go to the parent.
   
    // We'll update rangeStart and rangeEnd to the LCP array indices that
    // cut out the parent.
    rangeStart = getLCPPSV(lcpIndex);
    // Note that rangeEnd can be off the end of the LCP array now, if we
    // have moved up to the root node.
    rangeEnd = getLCPNSV(lcpIndex);
    
    // Update the Range object to describe this range, converting end
    // indices around. The range will never be empty so we don't have to
    // worry about negative end offsets.
    range.setForwardStart(rangeStart);
    range.setEndOffset(rangeEnd - rangeStart - 1);
    
    // Returnj the new pattern length (i.e. string depth)
    return lcp;
    
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
        extendFast(position, pattern[i], true);
    }
    
    // We either ran out of matching locations or finished the pattern.
    return position;

}

size_t FMDIndex::getLCP(size_t index) const {
    if(index >= getBWTLength()) {
        throw std::runtime_error("Looking at out-of-bounds LCP value!");
    }

    // Go get the longest common prefix length from the array.
    return lcpArray[index];
}
    
size_t FMDIndex::getLCPPSV(size_t index) const {
    
    if(index >= getBWTLength()) {
        throw std::runtime_error("Looking at out-of-bounds LCP PSV!");
    }
    
    // Go get the previous smaller value's index in the LCP array. Will
    // automatically handle if there isn't anything smaller.
    return lcpArray.getPSV(index);
}
    
size_t FMDIndex::getLCPNSV(size_t index) const {

    if(index >= getBWTLength()) {
        throw std::runtime_error("Looking at out-of-bounds LCP NSV!");
    }

    // Go get the next smaller value's index in the LCP array. Will
    // automatically handle if there isn't anything smaller.
    return lcpArray.getNSV(index);
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

int64_t FMDIndex::getContigEndIndex(size_t contig) const {
    // Looks a bit like the metadata functions from earlier. Actually pulls info
    // from the same file.
    return endIndices[contig];
}

char FMDIndex::display(int64_t index) const {
    // Just pull straight from the BWT string.
    return bwt.getChar(index);
}

char FMDIndex::display(size_t contig, size_t offset) const {
    // We need to loop through the contig from the back end until we get to the
    // right position.
    
    // How far from the back do we need to be?
    size_t backOffset = getContigLength(contig) - offset - 1;
    
    // Where are we in the BWT?
    int64_t bwtIndex = getContigEndIndex(contig);
    
    while(backOffset != (size_t) -1) {
        // Go left until we find the right letter.
        bwtIndex = getLF(bwtIndex);
        backOffset--;
    }
    
    return display(bwtIndex);
}

char FMDIndex::displayFirst(int64_t index) const {
    // Our BWT supports this natively.
    return bwt.getF(index);
}

std::string FMDIndex::displayContig(size_t index) const {
    // We can't efficiently un-locate, so we just use a vector of the last BWT
    // index in every contig. This works since there are no 0-length contigs.
    int64_t bwtIndex = getContigEndIndex(index);
    
    // Make a string to hold all the bases.
    std::string bases;

    for(size_t i = 0; i < getContigLength(index); i++) {
        // Until we have the right number of bases...
        
        // Grab this base
        bases.push_back(display(bwtIndex));
        
        // LF-map to the previous position in the contig (or off the front end).
        bwtIndex = getLF(bwtIndex);
    }
    
    // Flip the string around so it's front to front.
    std::reverse(bases.begin(), bases.end());
    
    // Give it back.
    return bases;
    
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

std::vector<Mapping> FMDIndex::map(const std::string& query,
    const GenericBitVector* mask, int minContext, int start, int length) const {
        
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }

    if(mask == NULL) {
        Log::debug() << "Mapping " << length << " bases to all genomes." <<
            std::endl;
    } else {
        Log::debug() << "Mapping " << length << " bases to one genome only." <<
            std::endl;
    }
    
    Log::debug() << "Mapping with minimum " << minContext << " context." <<
        std::endl;

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
        if(location.position.isEmpty(mask))
        {
            Log::debug() << "Starting over by mapping position " << i <<
                std::endl;
            // We do not currently have a non-empty FMDPosition to extend. Start
            // over by mapping this character by itself.
            location = this->mapPosition(query, i, mask);
        } else {
            Log::debug() << "Extending with position " << i << std::endl;
            // The last base either mapped successfully or failed due to multi-
            // mapping. Try to extend the FMDPosition we have to the right (not
            // backwards) with the next base.
            location.position = this->extend(location.position, query[i],
                false);
            location.characters++;
        }

        if(location.is_mapped && location.characters >= minContext &&
            location.position.getLength(mask) == 1) {
            
            // It mapped. We didn't do a re-start and fail, we have enough
            // context to be confident, and there's exactly one thing in our
            // interval.

            // Take the first (only) thing in the bi-interval's forward strand
            // side, not accounting for the mask.
            int64_t start = location.position.getForwardStart();
            
            if(mask != NULL) {
                // Account for the mask. The start position of the interval may
                // be masked out. Get the first 1 after (or at) the start,
                // instead of the start itself. Since the interval is nonempty
                // under the mask, we know this will exist.
                start = mask->valueAfter(start).first;
            }

            // Locate it, and then report position as a (text, offset) pair.
            // This will give us the position of the first base in the pattern,
            // which lets us infer the position of the last base in the pattern.
            TextPosition textPosition = locate(start);

            Log::debug() << "Mapped " << location.characters << "/" << 
                minContext << " context to text " << textPosition.getText() << 
                " position " << textPosition.getOffset() << std::endl;

            // Correct to the position of the last base in the pattern, by
            // offsetting by the length of the pattern that was used. A
            // 2-character pattern means we need to go 1 further right in the
            // string it maps to to find where its rightmost character maps.
            textPosition.setOffset(textPosition.getOffset() + 
                (location.characters - 1));

            // Add a Mapping for this mapped base, mapped on the right.
            Mapping mapped(textPosition);
            mapped.setMaxContext(0, location.characters);
            mappings.push_back(mapped);

            // We definitely have a non-empty FMDPosition to continue from

        } else {

            Log::debug() << "Failed (" << 
                location.position.getLength(mask) << " options for " <<
                location.characters << " context)." << std::endl;

            if(location.is_mapped && location.position.isEmpty(mask)) {
                // We extended right until we got no results. We need to try
                // this base again, in case we tried with a too-long left
                // context.

                Log::debug() << "Restarting from here..." << std::endl;

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

std::vector<Mapping> FMDIndex::mapRight(const std::string& query,
    const GenericBitVector* mask, int minContext) const {

    if(mask == NULL) {
        Log::debug() << "Mapping " << query.size() << 
            " bases to all genomes." << std::endl;
    } else {
        Log::debug() << "Mapping " << query.size() << 
            " bases to one genome only." << std::endl;
    }
    
    Log::debug() << "Mapping with minimum " << minContext << " context." <<
        std::endl;

    // We need a vector to return.
    std::vector<Mapping> mappings;

    // Start with the whole index selected.
    FMDPosition search = getCoveringPosition();
    
    // And with no characters searched
    size_t patternLength = 0;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--)
    {
        // For each position in the subrange we're mapping, stopping on
        // underflow
        
        // Try extending with that character
        FMDPosition extended(search);
        extendLeftOnly(extended, query[i]);
        
        
        while(extended.isEmpty(mask)) {
            // We would have no results if we extended with this character right
            // now.
            
            if(patternLength == 0) {
                // If you have no results after an extension, your original
                // pattern length has to be nonzero. We assume the index
                // contains at least one copy of every character.
            
                // TODO: Remove that assumption
                throw std::runtime_error("No results at zero pattern length! "
                    "Is a character not present in the index/genome?");
            }            
                        
            // Retract on the right until we have more results.
            patternLength = retractRightOnly(search);
            
            Log::debug() << "Retracted to length " << patternLength <<
                std::endl;
            
            // Try extending again until you do get results.
            extended = search;
            extendLeftOnly(extended, query[i]);
        
        }
        
        // Now you have some results. Adopt the new interval as your current
        // interval.
        search = extended;
        patternLength++;
        
        
        
        if(search.getLength(mask) == 1 && patternLength >= minContext) {
            // If you happen to have exactly one result with sufficient context,
            // record a mapping to it.
            
            // Take the first (only) search result.
            int64_t start = search.getForwardStart();
            
            if(mask != NULL) {
                // Account for the mask. The start position of the interval may
                // be masked out. Get the first 1 after (or at) the start,
                // instead of the start itself.
                start = mask->valueAfter(start).first;
            }

            // Locate it, and then report position as a (text, offset) pair.
            TextPosition textPosition = locate(start);

            Log::debug() << "Mapped " << patternLength << "/" << 
                minContext << " context to " << search << "; text " << 
                textPosition.getText() << " position " <<
                textPosition.getOffset() << std::endl;

            // Add a Mapping for this mapped base, with patternLength right
            // context.
            mappings.push_back(Mapping(textPosition, 0, patternLength));
        } else {
            // Otherwise record that this position is unmapped on the right.
            
            Log::debug() << "Failed: " << search.getLength(mask) << 
                " results for " << patternLength << "/" << minContext <<
                " context." << std::endl;
            
            mappings.push_back(Mapping());
        }
        
    }
        
    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());

    // Give back our answers.
    return mappings;

}

std::vector<Mapping> FMDIndex::mapRight(const std::string& query,
    int64_t genome, int minContext) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return mapRight(query, genome == -1 ? NULL : genomeMasks[genome],
        minContext);    
}

std::vector<Mapping> FMDIndex::mapLeft(const std::string& query,
    int64_t genome, int minContext) const {
    
    // Map the RC on the right
    std::vector<Mapping> mappings = mapRight(reverseComplement(query),
        genome, minContext);
        
    // Put them in proper base order.
    std::reverse(mappings.begin(), mappings.end());
    
    for(size_t i = 0; i < mappings.size(); i++) {
        // Go through all the mappings
        
        if(mappings[i].isMapped()) {
            // Flip the mapping for left semantics.
            mappings[i] = mappings[i].flip(getContigLength(getContigNumber(
                mappings[i].getLocation())));
        }
    }
    
    // Give back the fixed-up left-semantics mappings (so things that left-map
    // to the forward strand will be even texts).
    return mappings;  
}

std::vector<Mapping> FMDIndex::mapBoth(const std::string& query, int64_t genome, 
    int minContext) const {
    
    // Map it on the right. (text 0 = right-mapped to forward strand)
    std::vector<Mapping> right = mapRight(query, genome, minContext);
    
    // Map it on the left. (text 0 = left-mapped to forward strand)
    std::vector<Mapping> left = mapLeft(query, genome, 
        minContext);
    
    if(left.size() != right.size()) {
        throw std::runtime_error("Left and right size mismatch!");
    }
    
    for(size_t i = 0; i < left.size(); i++) {
        // Go through and disambiguate in place to resolve multi-mappings and
        // such. Make sure to read reverse backwards.
        
        left[i] = disambiguate(left[i], right[i]);
    }
    
    // Give back the disambiguated vector.
    return left;
    
}

std::vector<std::pair<int64_t,size_t>> FMDIndex::map(
    const GenericBitVector& ranges, const std::string& query, 
    const GenericBitVector* mask, int minContext, int addContext, 
    double multContext, double minCodingCost, int start, int length) const {
    
    // RIGHT-map to a range.
    
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }
    
    // Track the largest coding cost seen since restart.
    double codingCost;
    // And the state of the Markov model for calculating more coding costs.
    MarkovModel::iterator state = NULL;
    
    Log::debug() << "Mapping exactly with minimum " << minContext << 
        " and additional +" << addContext << ", *" << multContext << 
        " context, min coding cost " << minCodingCost << " bits." << std::endl;

    // We need a vector to return.
    std::vector<std::pair<int64_t,size_t>> mappings;

    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to trackour FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.
    MapAttemptResult location;
    // Make sure the scratch position is empty so we re-start on the first base
    location.position = EMPTY_FMD_POSITION;
    
    // Remember how many characters of context we have found after uniqueness.
    // Start it out at -1 so the first character we find making us unique brings
    // us to 0.
    int64_t extraContext = -1;

    for(int i = start + length - 1; i >= start; i--) {
        // Go from the end of our selected region to the beginning.

        Log::trace() << "On position " << i << " from " <<
            start + length - 1 << " to " << start << std::endl;

        // What range index does our current left-side position (the one we just
        // moved) correspond to, if any? Gets set in each branch of the if/else
        // below.
        int64_t range;
        
        // Keep track of whether that last thing we did was a restart or an
        // extend. This lets us know whether to try restarting from here again,
        // or just restart from the next base.
        bool lastRestarted;

        if(location.position.isEmpty(mask)) {
            Log::debug() << "Starting over by mapping position " << i <<
                std::endl;
            lastRestarted = true;
            
            // We do not currently have a non-empty FMDPosition to extend. Start
            // over by mapping this character by itself. This will return as
            // soon as it is unique.
            location = this->mapPosition(ranges, query, i, mask);
            
            // Reset the extra-context-after-uniqueness counter.
            extraContext = -1;
            
            // Look to see if we happen to be unique
            range = location.position.range(ranges, mask);
            
            if(location.is_mapped && !location.position.isEmpty(mask) &&
                range != -1) {
                
                // mapPosition finished happily and we have a unique result.
                extraContext = 0;
                
                while(i + location.characters < query.size()) {
                    
                    // We always need the maximal context later, and we can
                    // still extend on the right. Extend on the right more,
                    // making sure to count all these extensions as extra.
                    
                    Log::trace() << "Right extending with index " << 
                        i + location.characters << 
                        " to get maximal context (" << extraContext << "/" << 
                        addContext << " extra context)" << std::endl;
                    
                    // Where would we go if we extended?
                    FMDPosition nextPosition = this->extend(location.position, 
                        query[i + location.characters], false);
                        
                    if(!nextPosition.isEmpty(mask)) {
                        // Extension was successful
                        
                        location.position = nextPosition;
                        location.characters++;
                        extraContext++;
                        
                        // We can't do incremental coding cost updates here
                        // because we are going forwards.
                        
                    } else {
                        // Extension was not successful, so we've hit maximal
                        // context.
                        break;
                    }
                }
            }
            
            // Now we have a maximal context. Get its coding cost.
            
            if(markovModel != NULL) {
                // Get the coding cost of having what we've searched so far. We
                // searched from i out location.characters to the right. Also
                // updates the Markov state.
                codingCost = markovModel->backfill(state, reverseComplement(
                    query.substr(i, location.characters)));
            }
            
        } else {
            Log::debug() << "Extending with position " << i << std::endl;
            lastRestarted = false;
            // The last base either mapped successfully or failed due to multi-
            // mapping. Try to extend the FMDPosition we have to the left
            // (backwards) with the next base.
            location.position = this->extend(location.position, query[i], true);
            location.characters++;
            
            if(markovModel != NULL) {
                if(state == NULL) {
                    // We need to find a Markov state to start in.
                    codingCost = markovModel->backfill(state, reverseComplement(
                        query.substr(i, location.characters)));
                } else {
                    // We are already in a Markov state. Continue from there.
                    // Record the cost of this extension.
                    codingCost += markovModel->encodingCost(state, 
                        complement(query[i]));
                }
            }
            
            range = location.position.range(ranges, mask);
            if(range != -1 && extraContext == -1) {
                // We've just become unique. That will only ever happen here if
                // we aren't setting a minimum extra context, but we still need
                // it to be >= the default addContext of 0.
                extraContext = 0;
            }
            
        }
        
        if(location.is_mapped && !location.position.isEmpty(mask) &&
            range != -1 && location.characters >= minContext &&
            extraContext >= addContext && location.characters >= 
            (location.characters - extraContext) * multContext &&
            (markovModel == NULL || codingCost >= minCodingCost)) {
            // We have sufficient context to be confident (greater than the
            // minimum, with extraContext greater than the required additional
            // context, and the total context greater than the context
            // multiplier times the context taken to be unique), and our
            // interval is nonempty and subsumed by a range. And if we have a
            // Markov model the coding cost of everything in the context
            // (reverse complemented so it ends with what we just added) is
            // sufficiently high.

            Log::debug() << "Mapped " << location.characters << 
                " context (" << extraContext << "/" << addContext << 
                " extra) with multiplier " << multContext << 
                " to " << location.position << " in range #" << 
                range << std::endl;

            // Remember that this base mapped to this range
            mappings.push_back(std::make_pair(range, 
                location.characters));
            // We definitely have a non-empty FMDPosition to continue from
        } else {

            Log::debug() << "Failed at " << location.position << " (" << 
                location.position.ranges(ranges, mask) <<
                " options for " << location.characters << " context, " << 
                extraContext << " extra, " << 
                (location.characters - extraContext) * multContext << 
                " scaled, " << codingCost << " bits)." << std::endl;
                
            if(location.is_mapped && location.position.isEmpty(mask)) {
                // We extended right until we got no results. We need to try
                // this base again, in case we tried with a too-long left
                // context.

                Log::debug() << "Restarting from here..." << std::endl;

                // Move the loop index towards the end we started from (right)
                i++;

                // Since the FMDPosition is empty, on the next iteration we will
                // retry this base. The base is guaranteed to exist at least
                // once, so it won't be empty on the next pass.
            
                // TODO: check lastRestarted

            } else if(addContext > 0 && extraContext < addContext) {
                // Not enough context for additional context.
                
                // We need to have some amount of additional context when we do
                // map, and it needs to be out to the right of the position
                // we're mapping, after we get uniqueness going right. We didn't
                // have it this time, and we can't guess when we would have had
                // it by just extending left. So we need to restart to make sure
                // we can count it correctly.
                
                // We just restarted from here, so don't do it again.
                Log::debug() << 
                    "Restarting at next base to satisfy addContext..." << 
                    std::endl;
                    
                location.position = EMPTY_FMD_POSITION;
                
                // Incedentally, this position didn't map.
                mappings.push_back(std::make_pair(-1,0));
            
            } else if(multContext > 1 && location.characters < 
                (location.characters - extraContext) * multContext) {
                // Not enough context for multiplicative context scalar. We need
                // to restart to try to get more.
                
                if(lastRestarted) {
                    // We just restarted from here, so don't do it again.
                    Log::debug() << 
                        "Restarting at next base to satisfy multContext..." << 
                        std::endl;
                        
                    // This position didn't map.
                    mappings.push_back(std::make_pair(-1,0));
                        
                } else {
                    // We extended to get here, so we might find a closer
                    // uniqueness point and get more extra context if we restart
                    // from here instead.
                    Log::debug() << 
                        "Restarting from here to satisfy multContext..." << 
                        std::endl;
                        
                    // Try this position again
                    i++;
                }
                
                location.position = EMPTY_FMD_POSITION;
                
                
            } else {
                // It didn't map for some other reason:
                // - It was an initial mapping with too little right context to 
                //   be unique to a range.
                // - It was an initial mapping with a nonexistent right context
                // - It was an extension that was multimapped and still is

                // In none of these cases will re-starting from this base help
                // at all. If we just restarted here, we don't want to do it
                // again. If it was multimapped before, it had as much left
                // context as it could take without running out of string or
                // getting no results.

                // It didn't map. Say it corresponds to no range.
                mappings.push_back(std::make_pair(-1,0));

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

std::vector<int64_t> FMDIndex::mapRight(const GenericBitVector& ranges,
    const std::string& query, const GenericBitVector* mask, int minContext) const {
    
    // RIGHT-map to a range.
    
    Log::debug() << "Mapping with minimum " << minContext << " context." <<
        std::endl;

    // We need a vector to return.
    std::vector<int64_t> mappings;
    
    // Start with the whole index selected.
    FMDPosition search = getCoveringPosition();
    
    // And with no characters searched
    size_t patternLength = 0;
    
    for(size_t i = query.size() - 1; i != (size_t) -1; i--)
    {
        // For each position in the subrange we're mapping, stopping on
        // underflow
        
        // Try extending with that character
        FMDPosition extended = search;
        extendLeftOnly(extended, query[i]);
        
        
        while(extended.isEmpty(mask)) {
            // We would have no results if we extended with this character right
            // now.
            
            if(patternLength == 0) {
                // If you have no results after an extension, your original
                // pattern length has to be nonzero. We assume the index
                // contains at least one copy of every character.
            
                // TODO: Remove that assumption
                throw std::runtime_error("No results at zero pattern length! "
                    "Is a character not present in the index/genome?");
            }            
                        
            // Retract on the right. TODO: Integrate more deeply with
            // retractRightOnly because we will only ever want to retract to
            // points where we get more results.
            retractRightOnly(search, patternLength - 1);
            patternLength--;
            
            // Try extending again until you do get results.
            extended = search;
            extendLeftOnly(extended, query[i]);
        
        }
        
        // Now you have some results. Adopt the new interval as your current
        // interval.
        search = extended;
        patternLength++;
        
        // What range index does our current position correspond to, if any?
        int64_t range = search.range(ranges, mask);
        
        if(!search.isEmpty(mask) && range != -1 && patternLength >=
            minContext) {
            // If you happen to have results in exactly one range with
            // sufficient context, record a mapping to it.
            
            Log::debug() << "Mapped " << patternLength << " context to " << 
                search << " in range #" << range << std::endl;

            // Remember that this base mapped to this range
            mappings.push_back(range);
            
        } else {
            // Otherwise record that this position is unmapped on the right.
            
            Log::debug() << "Failed at " << search << " (" << 
                search.ranges(ranges, mask) <<
                " options for " << patternLength << " context)." << 
                std::endl;
            
            mappings.push_back(-1);
        }
        
    }
        
    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());

    // Give back our answers.
    return mappings;
}

std::vector<std::pair<int64_t,size_t>> FMDIndex::map(
    const GenericBitVector& ranges, const std::string& query, int64_t genome, 
    int minContext, int addContext, double multContext, double minCodingCost,
    int start, int length) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return map(ranges, query, genome == -1 ? NULL : genomeMasks[genome], 
        minContext, addContext, multContext, minCodingCost, start, length);    
}

std::vector<int64_t> FMDIndex::mapRight(const GenericBitVector& ranges, 
    const std::string& query, int64_t genome, int minContext) const {

    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return mapRight(ranges, query, genome == -1 ? NULL : genomeMasks[genome], 
        minContext);    
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
    size_t index, const GenericBitVector* mask) const {

    Log::trace() << "Mapping " << index << " in " << pattern << std::endl;
  
    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MapAttemptResult result;

    // Do a backward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.position = this->getCharPosition(pattern[index]);
    result.characters = 1;
    if(result.position.isEmpty(mask)) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if(result.position.getLength(mask) == 1) {
        // We've already mapped.
        result.is_mapped = true;
        return result;
    }

    if(index == 0) {
        // The rest of the function deals with characters to the left of the one
        // we start at. If we start at position 0 there can be none.
        return result;
    }

    Log::trace() << "Starting with " << result.position << std::endl;

    do {
        // Now consider the next character to the left.
        index--;

        // Grab the character to extend with.
        char character = pattern[index];

        Log::trace() << "Index " << index << " in " << pattern << " is " << 
            character << "(" << character << ")" << std::endl;

        // Backwards extend with subsequent characters.
        FMDPosition next_position = this->extend(result.position, character,
            true);

        Log::trace() << "Now at " << next_position << " after " << 
            pattern[index] << std::endl;
        if(next_position.isEmpty(mask)) {
            // The next place we would go is empty, so return the result holding
            // the last position.
            return result;
        } else if(next_position.getLength(mask) == 1) {
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

MapAttemptResult FMDIndex::mapPosition(const GenericBitVector& ranges, 
    const std::string& pattern, size_t index, const GenericBitVector* mask) const {
    
    
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
    if(result.position.isEmpty(mask)) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if (result.position.range(ranges, mask) != -1) {
        // We've already mapped.
        result.is_mapped = true;
        return result;
    }

    Log::trace() << "Starting with " << result.position << std::endl;
    
    for(index++; index < pattern.size(); index++) {
        // Forwards extend with subsequent characters.
        FMDPosition next_position = this->extend(result.position,
            pattern[index], false);

        Log::trace() << "Now at " << next_position << " after " << 
            pattern[index] << std::endl;
        if(next_position.isEmpty(mask)) {
            // The next place we would go is empty, so return the result holding
            // the last position.
            return result;
        }

        if(next_position.range(ranges, mask) != -1) {
            // We have successfully mapped to exactly one range. Update our
            // result to reflect the additional extension and our success, 
            // and return our result
          
            result.position = next_position;
            result.characters++;
            result.is_mapped = true;
            return result;
            
        } else {
            // Otherwise, we still map to a plurality of ranges. Record the
            // extension and loop again.
        
            result.position = next_position;
            result.characters++;
        
        }
    }
    
    return result;

}

Mapping FMDIndex::disambiguate(const Mapping& left, 
    const Mapping& right) const {
    
    // Make sure to disambiguate contexts along with mappings. Left/right
    // semantics agnostic, since both inputs must have the same which-text-is-
    // forward semantics.
    
    // Also make sure to combine contexts even if the result is not mapped.

    // Build a mapping that we will return.
    Mapping toReturn;

    if(!left.isMapped()) {
        // Copy the right mapping with its mapping location
        toReturn = right;
    } else if(!right.isMapped()) {
        // If right has nothing to say, use left
        toReturn = left;
    } else if(left.getLocation() == right.getLocation()) {
        // If they match, make sure to merge contexts, taking the left of left
        // and the right of right.
        toReturn = Mapping(right.getLocation());
    }
    
    // Else they disagree, so keep that default-constructed unmapped mapping.
    
    // Now fill in the contexts
    toReturn.setMinContext(left.getLeftMinContext(),
        right.getRightMinContext());
    toReturn.setMaxContext(left.getLeftMaxContext(),
        right.getRightMaxContext());
        
    return toReturn;
}

MisMatchAttemptResults FMDIndex::misMatchExtend(MisMatchAttemptResults& prevMisMatches,
        char c, bool backward, size_t z_max, const GenericBitVector* mask, bool startExtension, bool finishExtension) const {
    
    // Copy over to a new result
    MisMatchAttemptResults nextMisMatches;
    nextMisMatches.is_mapped = prevMisMatches.is_mapped;
    nextMisMatches.characters = prevMisMatches.characters;
    
    // Note that we do not flip parameters when !backward since
    // FMDIndex::misMatchExtend uses FMDIndex::extend which performs
    // this step itself
    
    if(prevMisMatches.positions.size() == 0) {
        throw std::runtime_error("Tried to extend when there are no positions to extend");
    }
    
    if(prevMisMatches.isEmpty(mask)) {
        throw std::runtime_error("Can't extend an empty position");
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
    
    std::pair<FMDPosition,size_t> m_position;
    std::pair<FMDPosition,size_t> m_position2;
        
    for(std::vector<std::pair<FMDPosition,size_t>>::iterator it =
          prevMisMatches.positions.begin(); it != prevMisMatches.positions.end(); ++it) {
        
        // Store each successive element as an m_position

        m_position.first = it->first;
        m_position.second = it->second;
                
        // extend m_position by correct base. Do not do this if the
        // finishExtension flag is true--in this case it's already been
        // done
    
        if(startExtension) {
                
            m_position2.first = extend(m_position.first, c, backward);
            m_position2.second = m_position.second;
        
            if(m_position2.first.getLength(mask) > 0) {
                nextMisMatches.positions.push_back(m_position2);        
            }
            
        } else if(finishExtension) {            
            // Extend by all mismatched positions
            if(m_position.second < z_max) {
                for(size_t base = 0; base < NUM_BASES; base++) {
                    if(BASES[base] != c) {
                        m_position2.first = extend(m_position.first, BASES[base], backward);
                        m_position2.second = m_position.second;
                        m_position2.second++;
                    
                        // If the position exists at all in the FMDIndex, place
                        // it in the results vector
                    
                        if(m_position2.first.getLength(mask) > 0) {
                            nextMisMatches.positions.push_back(m_position2);
                        }
                    }
                }
            }
        } else {
          
            m_position2.first = extend(m_position.first, c, backward);
            m_position2.second = m_position.second;
        
            if(m_position2.first.getLength(mask) > 0) {
                nextMisMatches.positions.push_back(m_position2);
                
            }
            
            if(m_position.second < z_max) {
                for(size_t base = 0; base < NUM_BASES; base++) {
                    if(BASES[base] != c) {
                        m_position2.first = extend(m_position.first, BASES[base], backward);
                        m_position2.second = m_position.second;
                        m_position2.second++;
                    
                        // If the position exists at all in the FMDIndex, place
                        // it in the results vector
                    
                        if(m_position2.first.getLength(mask) > 0) {
                                nextMisMatches.positions.push_back(m_position2);
                    
                        } 
                    }
                }
            }
        }                    
    }
    
    // If no results are found, place an empty FMDPosition in the
    // output vector
        
    if(nextMisMatches.positions.size() == 0) {
        nextMisMatches.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
    }
    
    // Return all matches
    
    // Or if there are matches, but not unique matches of at least
    // minimum context length, return the entire vector of positions
    // generated in this run, to use as starting material for the next
            
    return nextMisMatches;
    
}

std::vector<Mapping> FMDIndex::misMatchMap(
    const GenericBitVector& ranges, const std::string& query, 
    const GenericBitVector* mask, int minContext, int addContext,
    double multContext, double minCodingCost, size_t z_max, 
    bool keepIntermediates, int start, int length) const {
    
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }
        
    Log::debug() << "Mapping inexact (" << z_max << 
        " mismatches) with minimum " << minContext << " and additional +" << 
        addContext << ", *" << multContext << " context, min coding cost " << 
        minCodingCost << " bits." << std::endl;
        
    // We need a vector to return.
    std::vector<Mapping> mappings;
    
    // Start a set of mismatch search results with everything selected, no
    // characters searched and no mismatches.
    MisMatchAttemptResults search;
    search.positions.push_back(std::make_pair(getCoveringPosition(), 0));
    
    
    for(int i = start + length - 1; i >= start; i--) {
         // OK, I need to start at the right.
         
         // Loop invariant: search holds the result for searching up through the
         // previous position with all possible mismatches.
        
        Log::debug() << "On position " << i << " from " <<
            start + length - 1 << " to " << start << std::endl;
    
        // For each base, search out right until it becomes unique, and no
        // further. Don't consider any mismatches on the base itself.
        MisMatchAttemptResults mapUntilUnique = misMatchMapPosition(ranges, 
            query, i, z_max, false, false, mask);
        
        if(mapUntilUnique.is_mapped) {    
            Log::debug() << "Minimum right context: " << 
                mapUntilUnique.characters << std::endl;
        } else {
            Log::debug() << "No unique match." << std::endl;
        }
            
        // Set this if we have to throw away our current search and restart.
        bool restart = false;
            
        
        // Search interval can never be empty at the start of an iteration, by
        // the loop invariant.
        
        // Try extending left with the base we have. 
        MisMatchAttemptResults matchExtended = misMatchExtend(search, query[i],
            true, z_max, mask, true, false);
        matchExtended.characters++;
        
        if(matchExtended.range(ranges, mask) != -1) {
            // Flag it if it maps uniquely now.
            matchExtended.is_mapped = true;
        }
          
        // Make a Mapping to represent the result of this extension  
        Mapping mapping;
         
        if(!matchExtended.is_mapped) {
            Log::debug() << "Correct-base extension with " << query[i] << 
                " is not uniquely mapped" << std::endl;
            
            for(auto result : matchExtended.positions) {
                // Dump all the options.
                Log::trace() << "\t" << result.first << "\t" << result.second <<
                    std::endl;
            }
                
        } else if(matchExtended.characters < minContext) {
            Log::debug() << "Correct-base extension with " << query[i] << 
                " failed min context (" << matchExtended.characters << "/" << 
                minContext << ")" << std::endl;
                
        } else if(matchExtended.characters < 
            addContext + mapUntilUnique.characters) {
            
            Log::debug() << "Correct-base extension with " << query[i] << 
                " failed add context (" << matchExtended.characters << "/" << 
                addContext + mapUntilUnique.characters << ")" << std::endl;
        
        } else if(matchExtended.characters < multContext * 
            mapUntilUnique.characters) {
            
            Log::debug() << "Correct-base extension failed mult context (" <<
                matchExtended.characters << "/" << 
                multContext * mapUntilUnique.characters << ")" << std::endl;
        } else {
        
            Log::debug() << "Correct-base extension with " << query[i] <<
                " mapped and passed criteria" << std::endl;
        
            // Make a mapping from the right range. TODO: can this ever be
            // -1?
            mapping = Mapping(matchExtended.range(ranges, mask));
            
            // Keep the context that was the max we could go out to.
            mapping.setMaxContext(0, matchExtended.characters);
        }
        
        // Say there is some context that makes us unique, if there is one.
        mapping.setMinContext(0, mapUntilUnique.characters);
               
        // Then add the Mapping
        mappings.push_back(mapping);
        
        // Now compose the extension we will want to actually build on, by
        // extending with all characters.
        search = misMatchExtend(search, query[i], true, z_max, mask, false,
            false);
        search.characters++;
        
        Log::debug() << "Extending with all characters." << std::endl;
        
        for(auto result : search.positions) {
            // Dump all the options.
            Log::trace() << "\t" << result.first << "\t" << result.second <<
                std::endl;
        }
        
        
        if(search.isEmpty(mask)) {
            // Adding in this base kills off all the results we had. We need to
            // restart so we can try to find a non-empty interval. Go out as far
            // as you can without running out of results, and consider
            // mismatches on the original base.
            
            Log::debug() << "Search result set is empty. Restarting." <<
                std::endl;
            
            search = misMatchMapPosition(ranges, query, i, z_max, true, true,
                mask);
                
            if(search.isEmpty(mask)) {
                // We should always have at least some multimapping after a
                // restart.
                throw std::runtime_error("Restart got no results");
            }
            
        }
        
        // Repeat until the query string runs out.
        
    }

    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());

    // Give back our answers.
    return mappings;
}

std::vector<Mapping> FMDIndex::misMatchMap(
    const GenericBitVector& ranges, const std::string& query, int64_t genome, 
    int minContext, int addContext, double multContext, double minCodingCost,
    size_t z_max, bool keepIntermediates, int start, int length) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return misMatchMap(ranges, query, genome == -1 ? NULL : genomeMasks[genome], 
        minContext, addContext, multContext, minCodingCost, z_max,
        keepIntermediates, start, length);    
}

MisMatchAttemptResults FMDIndex::misMatchMapPosition(
    const GenericBitVector& ranges, const std::string& pattern, size_t index,
    size_t z_max, bool maxContext, bool allowFirstMismatch,
    const GenericBitVector* mask) const {
    
    // Reading right, note when(if ever) we become unique within z_max
    // mismatches. If maxContext is set, go all the way until we can't go
    // anymore instead. TODO: can we get both values out in that case? Or
    // continue an intermediate result right?
    
    // Not responsible for enforcing min context limits.
    
    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag that we
    // set if we are unique, and a variable counting the number of extensions
    // made to the FMDPosition.
    MisMatchAttemptResults result;
    
    // We aren't unique yet.
    result.is_mapped = false;
    // Start at the given index, and get the starting range for that character.
    // Add it in with 0 mismatches.
    result.positions.push_back(std::make_pair(getCharPosition(pattern[index]),
        0));
        
    if(allowFirstMismatch) {
        // We have to consider all the possible mismatches we could have at this
        // forst position as well.
        for(char other : BASES) {
            if(other == pattern[index]) {
                // Don't consider matches
                continue;
            }
            
            // Add in each other character's BWT range with 1 mismatch on it.
            result.positions.push_back(std::make_pair(getCharPosition(other),
                1));
        }
    }
        
    result.characters = 1;
    
    if(result.isEmpty(mask)) {
        // This character isn't even in it. Just return it as unmapped and let
        // the caller deal with it.
        return result;
    } else if (result.range(ranges, mask) != -1) {
        // We've already become unique.
        result.is_mapped = true;
        
        if(!maxContext) {
            // We only care about the minimum context for uniqueness
            // here.
            return result;
        }
        
        // In general we don't want to return yet; we need to go out to the
        // maximal context so that credit works correctly.
    }
    
    std::vector<std::pair<FMDPosition,size_t>> found_positions;
                
    for(index++; index < pattern.size(); index++) {
              
        // Forwards extend with subsequent characters.
      
        // Necessary here to create new result set since we're editing every
        // position of the last one. Is there a way around this? Don't think
        // so...
        
        // Extend right with all possibilities for the next character. This
        // tracks maximum mismatches for each FMDPosition and evicts the ones
        // that accumulate too many or otherwise become empty. Does not update
        // the character count.
        MisMatchAttemptResults new_result = this->misMatchExtend(result, 
            pattern[index], false, z_max, mask, false, false);
                        
        if(new_result.isEmpty(mask)) {
            // This is as far as we can go without running out of results,
            // whether we're mapped or not.
            return result;
        } else {
            // The new result is longer than ours but not empty, so take it.
            result = new_result;
            // Make sure to update the character count since the extend function
            // doesn't.
            result.characters++;
            
            // Have we uniquely mapped yet?
            // TODO: account for levels at which one range is not one node.
            if(result.range(ranges, mask) != -1) {
                // Flag it as having mapped
                result.is_mapped = true;
                
                if(!maxContext) {
                    // We only care about the minimum context for uniqueness
                    // here.
                    return result;
                }
            }
        }
    }
    
    // If we get here, we've run out of pattern, whether we became unique or
    // not.
    return result;

}


MisMatchAttemptResults FMDIndex::misMatchMapPosition(const GenericBitVector& ranges, 
    const std::string& pattern, size_t index, size_t minContext, 
    size_t addContext, double multContext, 
    int64_t* extraContext, size_t z_max, const GenericBitVector* mask, 
    bool maxContext) const {
    
    // We're going to right-map so ranges match up with the things we can map to
    // (downstream contexts)

    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MisMatchAttemptResults result;
    
    // To start with we haven't even become unique.
    *extraContext = -1;
            
    // Do a forward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.positions.push_back(std::pair<FMDPosition,size_t>(this->getCharPosition(pattern[index]),0));
    result.characters = 1;
    if(result.positions.front().first.isEmpty(mask)) {

        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        result.is_mapped = true;    
        return result;
    } else if (result.positions.front().first.range(ranges, mask) != -1) {
        // We've already become unique.

        *extraContext = 0;
        
        if(minContext <= 1 && addContext == 0 && multContext <= 1) {
            // And we just need to be unique to map
        
            result.is_mapped = true;
            
            if(!maxContext) {
                // We only care about the minimum context for uniqueness
                // here.
                return result;
            }
            
            // In general we don't want to return yet; we need to go out to the
            // maximal context so that credit works correctly.
        }
    }
    
    std::vector<std::pair<FMDPosition,size_t>> found_positions;
                
    for(index++; index < pattern.size(); index++) {
              
        // Forwards extend with subsequent characters.
      
        // Necessary here to create new result set since we're editing every
        // position of the last one. Is there a way around this? Don't think so...
            
        MisMatchAttemptResults new_result = this->misMatchExtend(result, pattern[index], false, z_max, mask, false, false);
                        
        if(new_result.positions.front().first.isEmpty(mask)) {
            // We can't extend any more and have results, so we have reached
            // maximal context length and need to return.
        
            if(result.positions.size() == 1 && 
                result.characters >= minContext && 
                *extraContext >= addContext && result.characters >= 
                (result.characters - *extraContext) * multContext) {
                
                Log::debug() << "Maximal context length reached" << std::endl;
                
                result.is_mapped = true;
                return result;
            } else {
                // Last position multimapped but this extension now
                // maps nowhere. Return an empty set of positions
              
                // Or we don't map anywhere with at least the minimum
                // context.
                
                Log::debug() << "Can't get sufficient context to map" << 
                    std::endl;
                
                result.positions.clear();
                result.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
                result.is_mapped = false;
                result.characters = 1;
                return result;
            }
        }
        
        Log::debug() << "Extra context: " << *extraContext << std::endl;
        
        if(!result.is_mapped && new_result.positions.front().first.range(ranges, mask) != -1 &&
          new_result.positions.size() == 1) {
            // We will successfully map to exactly one range. Update our result
            // to reflect the additional extension.
            
            // However, we may not have enough context to report this mapping,
            // so sit on it by default.
            
            // If this it the first time through, we have 0 extra context (from
            // -1). Otherwise, this counts as extra context even if we haven't
            // met the min context.
            (*extraContext)++;
            
            
            
            result.positions = new_result.positions;
            result.characters++;
            
            if(result.characters >= minContext) {
                // We do have enough context to report this mapping.
                
                result.is_mapped = true;
                found_positions = result.positions;
                
                Log::debug() << "Will become unique at " << index  << 
                    " with " << result.characters << " context" << std::endl;
                    
                if(!maxContext) {
                    // We only care about the minimum context for uniqueness
                    // here.
                    return result;
                }
                    
            } else {
                Log::debug() << 
                    "Will become unique with insufficient context (" << 
                    result.characters << ") at " << index << std::endl;
            }   
        } else if(result.is_mapped && new_result.positions.front().first.range(ranges, mask) != -1) {
            
            // We have mapped to a single range after previously having done
            // that.
            
            // Add 1 to the extra context after being unique.
            (*extraContext)++;
            
            result.positions = new_result.positions;
            result.characters++;
            
            Log::debug() << "Will add context (for " << result.characters << 
                ") at index " << index << std::endl;
            
        } else {
          
            // Otherwise, we still map to a plurality of ranges. Record the
            // extension and loop again.
            
            result.positions = new_result.positions;
            result.characters++;
            
            Log::debug() << "Will multimap at index " << index << " with " << 
                result.characters << " context" << std::endl;
        }          
   }

    if (result.is_mapped && *extraContext >= addContext && 
        result.characters >= 
        (result.characters - *extraContext) * multContext &&
        result.characters >= minContext) {
        
        Log::debug() << "Hit end of pattern with sufficient context" <<
            std::endl;
        
        result.positions = found_positions;
    } else {
        
        // If we get here, we ran out of downstream context and still map to
        // multiple ranges and/or don't meet the min context requirements. Send
        // back a result indicating no mapping.
        
        Log::debug() << "Hit end of pattern with insufficient context" <<
            std::endl;
        
        result.positions.clear();
        result.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
        result.is_mapped = false;
    }
        
    return result;

}

// In case anyone wants it later... the following two functions implement a
// mismatch extend which returns results sorted by number of mismatches

MisMatchAttemptResults FMDIndex::sortedMisMatchExtend(MisMatchAttemptResults& prevMisMatches,
        char c, bool backward, size_t z_max, const GenericBitVector* mask) const {
    MisMatchAttemptResults nextMisMatches;
    nextMisMatches.is_mapped = false;
    nextMisMatches.characters = prevMisMatches.characters;
    
    // Note that we do not flip parameters when !backward since
    // FMDIndex::misMatchExtend uses FMDIndex::extend which performs
    // this step itself
    
    if(prevMisMatches.positions.size() == 0) {
        throw std::runtime_error("Tried to extend zero length mismatch vector");
    }
    
    if(prevMisMatches.positions.front().first.isEmpty(mask)) {
        throw std::runtime_error("Can't extend an empty position");
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
    
    // z tracks how many mismatches each range has in order to
    // sort our range-holding data structure. We initialize z with
    // the minimum number of mismatches for which a context existed
    // the previous extension
        
    size_t z = prevMisMatches.positions.front().second;
    
    // To store our FMDPositions under consideration in order of
    // number z of mismatches so that we place them on the queue
    // in the correct order
    
    std::vector<std::pair<FMDPosition,size_t>> waitingMatches;
    std::vector<std::pair<FMDPosition,size_t>> waitingMisMatches;
    
    std::pair<FMDPosition,size_t> m_position;
    std::pair<FMDPosition,size_t> m_position2;
        
    // Output has flag to mark whether we have found a unique
    // hit at the most favourable z-level, which has our
    // minimum context length, in which case we want to terminate
    // search. If not flagged, we want to pass the entire queue
    // as working material for the next extension since we don't
    // know which "mismatch path" will find us our end result
    
        
    for(std::vector<std::pair<FMDPosition,size_t>>::iterator it =
        prevMisMatches.positions.begin(); it != prevMisMatches.positions.end(); ++it) {
        // Store each successive element as an m_position

        m_position.first = it->first;
        m_position.second = it->second;
        
        // Check if we exhausted the search over all sequences in the
        // queue with z mismatches. If so, search all exact-base
        // extensions of level-z sequences to see if there is a single unique
        // hit. Else add all such sequences to the queue, and next search
        // all mismatched extensions for a unique hit.
        
        if(m_position.second != z) {
          
            // Built-in check to make sure our data structure holding all extended
            // context ranges has actually been arranged in mismatch number order
            // by the previous iteration of processMisMatchPositions
          
            if(m_position.second < z) {
                throw std::runtime_error("Generated misordered mismatch list");
            }
            
            // Evaluate and search all extensions 
                                
            processMisMatchPositions(nextMisMatches, waitingMatches, waitingMisMatches, mask);
            
            // Else we need to see what we get at subsequent z-levels
            z++;
            
        }

        // extend m_position by correct base
                
        m_position2.first = extend(m_position.first, c, backward);
        m_position2.second = z;
        waitingMatches.push_back(m_position2);
        
        if(z < z_max) {
            for(size_t base = 0; base < NUM_BASES; base++) {
                if(BASES[base] != c) {
                    m_position2.first = extend(m_position.first, BASES[base], backward);
                    m_position2.second = z;
                    m_position2.second++;
                    waitingMisMatches.push_back(m_position2);
                
                }
            }
        }
    }
    
    
    // We have extended and searched the entire queue of matches from the
    // last level
        
    processMisMatchPositions(nextMisMatches, waitingMatches, waitingMisMatches, mask);

    
    // If there are matches, but not unique matches of at least
    // minimum context length, return the entire "heap"
    // generated in this run, to use as starting material for the next
    
    // If there are unique matches, these will also get passed. Don't
    // need the flag 
        
    if(nextMisMatches.positions.size() == 0) {
        nextMisMatches.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
        nextMisMatches.is_mapped = 0;
    }
            
    return nextMisMatches;
}

void FMDIndex::processMisMatchPositions(
                MisMatchAttemptResults& nextMisMatches,
                std::vector<std::pair<FMDPosition,size_t>>& waitingMatches,
                std::vector<std::pair<FMDPosition,size_t>>& waitingMisMatches,
                const GenericBitVector* mask) const {
                  
    while(!waitingMatches.empty()) {
        if(waitingMatches.back().first.getLength(mask) > 0) {
            nextMisMatches.positions.push_back(waitingMatches.back());
        }
      
        waitingMatches.pop_back();
    }
  
    while(!waitingMisMatches.empty()) {
        if(waitingMisMatches.back().first.getLength(mask) > 0) {
            nextMisMatches.positions.push_back(waitingMisMatches.back());
        }
      
        waitingMisMatches.pop_back();
  }
  
  return;
  }
