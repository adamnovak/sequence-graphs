"""
fast.py
Parsers and writers for FASTA, FASTQ, and QUAL files.

All functions in this file work with Sequence objects (which carry a string 
"identifier", some "extraInfo" in a string, and a list "letters" of 
single-character strings). Some functions require or produce more specialized 
QualitySequence objects, which additionally have a list "qualities" of Phred 
quality scores for each base. The user of the module is responsible for making 
sure that Sequence and QualitySequence objects passed in are semantically valid.
A Sequence or QualitySequence obtained from this module will always be 
semantically valid unless modified by the user.

All readers are Python generators which take one or more streams, and 
potentially some optional control parameters, and yield Sequence or 
QualitySequence objects, as appropriate, representing the sequences read. All 
syntactic errors and some semantic errors are detected and reported with line 
numbers.

All writers take a Python iterable (referred to as a generator, but anything 
that can be looped over will work) that produces either Sequence or 
QualitySequence objects (as appropriate for the type(s) of file(s) being 
written), one or more streams to write to, and potentially some optional control
parameters. They write sequences obtained from the iterable to the specified 
stream(s) in the appropriate format. The writers will detect and report some 
semantic errors by sequence ID, but it is in general possible to get them to 
produce syntactically invalid output by providing them with flawed Sequence or 
QualitySequence objects.

Written for BME 205 Fall 2012
By Adam Novak <anovak1@ucsc.edu>
Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import re, itertools

class Sequence(object):
    """
    Represents a sequence, as might be found in a FASTA file.
    Carries an identifier, an info line and a list holding the sequence letters.
    """
    def __init__(self, identifier, extraInfo):
        """
        Make a new Sequence with the specified identifier and extra info.
        Letters will need to be added to it manually for it to be useful.
        """
        # The identifier (between > and first " " or , in a FASTA)
        self.identifier = identifier
        # The extra info (everything on the header line after the identifier and
        # separator)
        self.extraInfo = extraInfo
        # Sequence letters live here. It's a list so we can append efficiently.
        self.letters = []
        
    def __str__(self):
        """
        Produce a string describing this sequence.
        """
        return "<Sequence {}: {} letters>".format(self.identifier, 
            len(self.letters))

class QualitySequence(Sequence):
    """
    Represents a sequence with quality information.
    """
    def __init__(self, identifier, extraInfo):
        """
        Make a new Sequence with the specified identifier and extra info.
        Letters will need to be added to it manually for it to be useful.
        """
        # Inheritance!
        # See http://parand.com/say/index.php/2009/04/20
        # /python-simple-inheritance-example/
        super(QualitySequence, self).__init__(identifier, extraInfo)
        
        # Sequence quality scores live here (list of ints)
        self.qualities = []
        
    def __str__(self):
        """
        Produce a string describing this sequence.
        """
        return "<QualitySequence {}: {} letters>".format(self.identifier, 
            len(self.letters))

def read_fasta(stream, alphabet=r"[A-IK-Za-ik-z\-*]"):
    """
    Reads a FASTA file from the stream passed in in stream.
    Yields each sequence as a Sequence object. FASTA ;-flagged comments are 
    supported, and headers starting with > are required. All characters not in 
    the FASTA alphabet [A-IK-Za-ik-z\-*] are stripped from sequences, and 
    lower-case letters are converted to upper-case.
    
    The argument "alphabet" can be a regular expression string to override the 
    default alphabet of acceptable characters (e.g. include J).
    """
    
    # This holds the current sequence object if we are in the process of reading
    # a sequence, or None if we are not.
    currentSequence = None
    
    # This regex matches a header, and parses out the ID and the extra info
    # All headers start with >
    # Then any number of non-space, non-comma characters (greedy)
    # Then a space or a comma (optional)
    # Then anything
    # Nothing will be in the last part (extra info) unless there was a space or 
    # comma.
    headerRegex = re.compile(r"^>([^,\s]*)(?:,|\s)?(.*)") 
    
    # This regex matches each character that is acceptable in a sequence.
    alphabetRegex = re.compile(alphabet)
    
    # This regex matches comments
    commentRegex = re.compile(r"^;")
    
    # This tracks line number for error reporting
    lineNumber = 0
    
    for line in stream:
        lineNumber = lineNumber + 1
        
        # Skip comments
        if commentRegex.match(line):
            continue
        
        # Keep the regex match result, since it has the capture groups
        headerMatch = headerRegex.match(line)
        if headerMatch:
            # We have a new header!
            # If we were already parsing a sequence, that's done, so yield it
            if currentSequence:
                yield currentSequence
            
            # Pull out the capture groups and make a new sequence
            currentSequence = Sequence(headerMatch.group(1), 
                headerMatch.group(2))
            
        else:
            if currentSequence is None:
                raise Exception("FASTA missing header at "
                    "line {}".format(lineNumber))
                
            # This line goes in the current sequence
            for acceptableLetter in alphabetRegex.findall(line):
                # Be sure to capitalize all the lower-case letters
                currentSequence.letters.append(acceptableLetter.capitalize())
                
    # File over
    # If we got any sequences, we still need to yield the last one
    if currentSequence is not None:
        yield currentSequence  

def write_fasta(generator, stream):
    """
    Writes the Sequence objects yielded by the generator generator to the stream
    given in stream, in FASTA format.
    """
    
    for sequence in generator:
        stream.write(">{} {}\n".format(sequence.identifier, 
            sequence.extraInfo))
            
        # Holds letter count, so we can insert a line break every 80 letters
        letterCount = 0
        
        for letter in sequence.letters:
            stream.write(letter)
            letterCount = letterCount + 1
            if letterCount % 80 == 0:
                stream.write("\n")
                
        # End each sequence with a newline, even if we just got to 80 characters
        stream.write("\n")
           
def read_qualities(stream):
    """
    Reads a quality file (FASTA headers and whitespace-separated quality 
    scores). Ignores the headers and yields lists of ints.
    """         
    # This holds the current quality list if we are in the process of reading
    # a list, or None if we are not.
    currentList = None
    
    # This regex matches a header. We ignore the extra header info.
    headerRegex = re.compile(r"^>") 
    
    # This regex matches a positive integer, which is what the quality scores
    # are
    intRegex = re.compile(r"[0-9]+")
    
    # This tracks line number for error reporting
    lineNumber = 0
    
    for line in stream:
        lineNumber = lineNumber + 1
        if headerRegex.match(line):
            # We have a new header!
            # If we were already parsing a sequence, that's done, so yield it
            if currentList:
                yield currentList
            
            # Start a new list
            currentList = []
            
        else:
            if currentList is None:
                raise Exception("Quality file missing header "
                    "at line {}".format(lineNumber))
                
            # This line goes in the current sequence
            for scoreString in intRegex.findall(line):
                # Make sure score is non-negative
                # This holds the score as an integer
                phred = int(scoreString)
                if phred < 0:
                    raise Exception("Negative quality {} "
                        "at line {}".format(phred, lineNumber))
                currentList.append(phred)
                
    # File over
    # If we got any quality lists, we still need to yield the last one
    # List may be empty: compare to None explicitly
    if currentList is not None:
        yield currentList
        
def write_qualities(generator, stream):
    """
    Write the quality information for QualitySequences yielded by the generator
    to the given stream, in quality file format.
    """
    
    for sequence in generator:
        stream.write(">{} {}\n".format(sequence.identifier, 
            sequence.extraInfo))
            
        # Numbers can be at most 3 digits, or 4 characters per number with a 
        # separator. We limit line length by inserting a line break every 20 
        # numbers.
        
        # This variable keeps track of the number of numbers we have written so 
        # far.
        numberCount = 0
        
        for quality in sequence.qualities:
            stream.write(str(quality))
            numberCount = numberCount + 1
            if numberCount % 20 == 0:
                stream.write("\n")
            else:
                stream.write(" ")
        
        # End each sequence with a newline, even if we just got to 20 numbers
        stream.write("\n")
                
def read_fasta_with_quality(fastaStream, qualityStream):
    """
    Reads FASTA sequences from the fastaStream stream, and quality information 
    in quality file format from the qualityStream stream. Yields QualitySequence
    objects for each sequence from the two streams. Both streams are assumed to
    contain information for the same sequences in the same order, with matching
    header lines.
    """

    # This iterator yields Sequence objects from the FASTA
    fasta_iterator = read_fasta(fastaStream)
    
    # This iterator yields lists of qualities, corresponding to the above
    quality_iterator = read_qualities(qualityStream)

    # sequence holds the Sequence object, and qualities the corresponding 
    # quality list.
    # Extra entries in either file are dropped.
    for (sequence, qualities) in itertools.izip(fasta_iterator, quality_iterator):
        # This is the QualitySequence we will yield
        qualitySequence = QualitySequence(sequence.identifier, 
            sequence.extraInfo)
            
        qualitySequence.letters = sequence.letters
        qualitySequence.qualities = qualities
        
        if len(qualitySequence.letters) != len(qualitySequence.qualities):
            raise Exception("Sequence length and quality length do not match "
                "for sequence '{}'".format(qualitySequence.identifier))
        
        yield qualitySequence
        
def write_fasta_with_quality(generator, fastaStream, qualityStream):
    """
    Write a FASTA file and a corresponding quality file. generator must be a 
    generator that yields QualitySequence objects. fastaStream and qualityStream
    are output streams for FASTA sequence data and quality data, respectively.
    """
    
    # Problem: if we just tee the generator, and run write_fasta and 
    # write_qualities in sequence, we'll try to keep everything from the 
    # generator in memory.
    
    # Solution 1: munge iterators by hand and repeatedly call write functions
    # Solution 2: multithreading
    
    # Here we implement solution 1.
    # The ugliness here will hopefully be worth it when everything has a 
    # consistent API later.
    
    for qualitySequence in generator:
        if len(qualitySequence.letters) != len(qualitySequence.qualities):
            raise Exception("Sequence length and quality length do not match "
                "for sequence '{}'".format(qualitySequence.identifier))
                
        # Lists count as generators for our purposes
        write_fasta([qualitySequence], fastaStream)
        write_qualities([qualitySequence], qualityStream)
        
def read_fastq(stream, qualityOffset=33):
    """
    Takes in an input stream with FASTQ sequence-and-quality data. Yields
    QualitySequence objects. qualityOffset stores the value to subtract from 
    quality character ASCII values to get the quality score (typically 33, 
    sometimes 64).
    """
    
    # It would be quite hard to parse this without some sort of state machine,
    # so let's make that explicit.
    # parsingState can be "header", "sequence", or "quality", and
    # tells what the parser is currently trying to read.
    parsingState = "header"
    
    # This holds the current QualitySequence object we are populating.
    currentSequence = None
 
    # This regex matches a header, and parses out the ID and the extra info
    # All headers start with @
    # Then any number of non-space, non-comma characters (greedy)
    # Then a space or a comma (optional)
    # Then anything
    # Nothing will be in the last part (extra info) unless there was a space or 
    # comma.
    headerRegex = re.compile(r"^@([^,\s]*)(?:,|\s)?(.*)") 
    
    # This regex matches the sequence/quality separator line
    # We ignore the extra optional information on this line.
    separatorRegex = re.compile(r"^\+")
    
    # This regex matches each character that is acceptable in sequence data.
    # All letters (both cases) but J, and - and *. "." is not allowed and we  
    # should ignore it.
    alphabetRegex = re.compile(r"[A-IK-Za-ik-z\-*]")
    
    # This regex matches valid quality characters.
    # Any character but whitespace is accepted
    qualityRegex = re.compile(r"\S")
    
    # This variable stores the current line number for error reporting
    lineNumber = 0
    
    for line in stream:
        lineNumber = lineNumber + 1
        if parsingState == "header":
            headerMatch = headerRegex.match(line)
            if headerMatch:
                # We have a new header!
                # Pull out the capture groups and make a new sequence
                currentSequence = QualitySequence(headerMatch.group(1), 
                    headerMatch.group(2))
                    
                # Go to "sequence" state
                parsingState = "sequence"
            else:
                raise Exception("Expected header at "
                    "fastq line {}".format(lineNumber))   
        elif parsingState == "sequence":
            if separatorRegex.match(line):
                # Go to "quality" state
                parsingState = "quality"
            else:
                for acceptableLetter in alphabetRegex.findall(line):
                    # Be sure to capitalize all the lower-case letters
                    currentSequence.letters.append(
                        acceptableLetter.capitalize())
        elif parsingState == "quality":
            for qualityCharacter in qualityRegex.findall(line):
                # Convert from single byte to integer
                # phred stores the actual integer quality Phred score
                phred = ord(qualityCharacter) - qualityOffset
                if phred < 0:
                    raise Exception("Negative quality score "
                        "{} decoded at fastq line {}".format(phred, lineNumber))
                currentSequence.qualities.append(phred)
            if len(currentSequence.qualities) == len(currentSequence.letters):
                # Done parsing all the qualities for this sequence
                yield currentSequence
                currentSequence = None
                parsingState = "header"
            elif len(currentSequence.qualities) > len(currentSequence.letters):
                raise Exception("More quality scores than sequence letters "
                    "at fastq line {}".format(lineNumber))
        else:
            raise Exception("Invalid parser state {} "
                "at fastq line {}".format(parsingState, lineNumber))
    
    # If we reach the end of the file and we aren't in "header" state, that's a 
    # bad file
    if parsingState != "header":
        raise Exception("File ended when expecting "
            "{} at fastq line {}".format(parsingState, lineNumber))

def write_fastq(generator, stream, qualityOffset=33):
    """
    Takes QualitySequence objects yielded by the generator, and writes them in 
    FASTQ format to the output stream in stream. FASTQ sequence and quality data
    for each sequence are written as single lines. qualityOffset controls the 
    offset used to convert Phred scores into printable ASCII character values. 
    It should usually be 33, but can sometimes be 64.
    """
    
    for qualitySequence in generator:
        # Header
        stream.write("@{} {}\n".format(qualitySequence.identifier, 
            qualitySequence.extraInfo))
            
        # Body
        map(stream.write, qualitySequence.letters)
        stream.write("\n+\n")
        for phred in qualitySequence.qualities:
            # Pre-calculate the char value so we can check it for validity
            charValue = phred + qualityOffset
            if charValue > 255:
                raise Exception("Out of range quality value "
                "{} in sequence '{}'".format(phred, qualitySequence.identifier))
            stream.write(chr(phred + qualityOffset))
        stream.write("\n")
        
