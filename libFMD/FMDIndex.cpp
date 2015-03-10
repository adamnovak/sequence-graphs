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

thread_local std::map<const FMDIndex*, std::map<size_t,
    std::map<size_t, std::string>::iterator>> FMDIndex::threadContigCache;

FMDIndex::FMDIndex(std::string basename, SuffixArray* fullSuffixArray): names(),
    starts(), lengths(), cumulativeLengths(), genomeAssignments(), endIndices(),
    genomeRanges(), genomeMasks(), bwt(basename + ".bwt"),
    suffixArray(basename + ".ssa"), fullSuffixArray(fullSuffixArray),
    lcpArray(basename + ".lcp"), contigCache(), contigCacheMutex() {
    
    // TODO: Too many initializers

    Log::info() << "Loading " << basename << std::endl;

    // We already loaded the index itself in the initializer. Go load the
    // length/order metadata.
    
    // Open the contig name/start/length/genome file for reading. Temporary
    // string survives until the end of the full expression.
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
    // Drop any entries for this object from the thread-local static cache
    // (which is supposed to be simulating thread-local members).
    clearThreadLocalCache();

    if(fullSuffixArray != NULL) {
        // If we were holding a full SuffixArray, throw it out.
        delete fullSuffixArray;
    }
    
    for(std::vector<GenericBitVector*>::iterator i = genomeMasks.begin(); 
        i != genomeMasks.end(); ++i) {
        
        // Also delete all the genome masks we loaded.
        delete (*i);
    }
}

size_t FMDIndex::getContigNumber(TextPosition base) const {
    // What contig corresponds to that text? Contigs all have both strands.
    return base.getContigNumber();
}

bool FMDIndex::getStrand(TextPosition base) const {
    // What strand corresponds to that text?
    return base.getStrand();
}

size_t FMDIndex::getContigOffset(TextPosition base) const {
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
    size_t offset = getContigOffset(base);
    
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
    return total + getContigOffset(base) - 1;
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

char FMDIndex::displayCached(const TextPosition& position) const {
    
    // Grab the reference contig in the cache
    const std::string& referenceContig = displayContigCached(
        position.getContigNumber());
        
    // Grab the text-local offset
    size_t offset = position.getOffset();
    
    if(position.getStrand()) {
        // We need to get counting from the other end and complement the base.
        return complement(referenceContig[referenceContig.size() - offset - 1]);
    } else {
        // Just count from the front.
        return referenceContig[offset];
    }
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

const std::string& FMDIndex::displayContigCached(size_t index) const {
    // Try to find the thing we are looking for
    auto found = threadContigCache[this].find(index);
    
    if(found == threadContigCache[this].end()) {
        // If we don't have a copy thread-local, we'll need to go to the all-
        // threads cache for this FMDIndex, so we need to synchronize.
        
        // Lock the cache
        std::unique_lock<std::mutex> lock(contigCacheMutex);
        
        // Look in it
        auto globalFound = contigCache.find(index);
        
        if(globalFound == contigCache.end()) {
            // Go get it and insert it, updating the iterator for finding it.
            globalFound = contigCache.insert({index,
                displayContig(index)}).first;
        }
        
        // Insert the iterator in the global cache into the local cache.
        found = threadContigCache[this].insert({index, globalFound}).first;
    }
    
    // Return a reference to the globally cached copy, which means we have to
    // dereference twice.
    return found->second->second;
}

void FMDIndex::clearThreadLocalCache() const {
    // Remove the entry for this object in the global-scope thread-local cache,
    // in case we later have to work with a different object at the same
    // address.
    threadContigCache.erase(this);
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

FMDIndex::iterator FMDIndex::begin(size_t depth, bool reportDeadEnds) const {
    // Make a new suffix tree iterator that automatically searches out the first
    // suffix of the right length.
    return FMDIndex::iterator(*this, depth, false, reportDeadEnds);
}
     
FMDIndex::iterator FMDIndex::end(size_t depth, bool reportDeadEnds) const {
    // Make a new suffix tree iterator that is just a 1-past-the-end sentinel.
    return FMDIndex::iterator(*this, depth, true, reportDeadEnds);
}
