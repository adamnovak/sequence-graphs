#!/usr/bin/env python2.7
"""
unitigStats.py: Read a (possibly gzipped) unitig .mag file, and report some
statistics about it.

MAG format is a variant of FASTQ format, with specially formatted header lines
and quality scores.

The header of a record contains four fields, separated by tab characters. The
first field gives the numerical identifiers for the two ends of the sequence,
separated by a colon. The second field gives the number of reads that went into
the record. The third field, if it is not ".", gives semicolon-terminated end-
id,offset pairs detailing the end with which the left end of this record
overlaps, and the number of bases of overlap. The fourth field is the same, but
for the right end of this sequence.

Example record header line:

    @3056104052:2484561266	1040	1855531145,98;	743607148,86;183408,88;

The qualities in the FASTQ record give the number of reads covering each base,
in chr(33 + value) format.

Note that the + line for qualities will not contain the same header information.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, zlib

import Bio.SeqIO


def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("in_file", type=argparse.FileType("r"),
        default=sys.stdin,
        help="MAG-format file of the unitig graph to read")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

class TransparentUnzip(object):
    """
    A class that represents a transparently-un-gzipping stream.
    
    Right now only supports readline
    """
    
    def __init__(self, stream):
        """
        Encapsulate the given stream in an unzipper.
        """
        
        # Keep the stream
        self.stream = stream
        
        # Make a decompressor with the accept a gzip header flag.
        # See <http://stackoverflow.com/a/22311297/402891>
        self.decompressor = zlib.decompressobj(zlib.MAX_WBITS + 16)
        
        # We need to do lines ourselves, so we need to keep a buffer
        self.line_buffer = ""
        
    def readline(self):
        """
        Return the next line with trailing "/n", or "" if there is no next line.
        """
        
        print("Buffered: {} bytes".format(len(self.line_buffer)))
        
        # See if we have a line to spit out
        newline_index = self.line_buffer.find("\n")
        
        if newline_index == -1:
            # No line is in the buffer
            print("No line")
            
            # Go get more data from the stream (in 16 k blocks)
            compressed = self.stream.read(16 * 2 ** 10)
            
            #print("Compressed: {}".format(compressed))
            
            if compressed == "":
                print("No compressed data")
                if self.decompressor is not None:
                    # No mor einput data; flush out the output data.
                    self.line_buffer += self.decompressor.flush()
                    self.decompressor = None
                    # Spit out a line from that
                    print("Flushed decompressor")
                    return self.readline()
                
                # We didn't find a newline, and there's no more data, and
                # nothing to flush. Take what we have (with no trailing \n)
                line = self.line_buffer
                # Clear the buffer so next time we return "" as desired.
                self.line_buffer = ""
                
                return line
                
            # Otherwise we found more data
            
            try:
                # Decompress each block if needed
                decompressed = self.decompressor.decompress(compressed)
            except zlib.error:
                # Just skip decompressing; it's probably not actually
                # compressed.
                decompressed = compressed
            
            print("Got {} bytes compressed, {} bytes expanded".format(
                len(compressed), len(decompressed)))
            
            # Stick it in the buffer
            self.line_buffer += decompressed
            
            # Try again
            return self.readline()
        
        else:
            print("Returning line")
            # We have a line. Grab it.
            line = self.line_buffer[0:newline_index + 1]
            
            # Pop it off
            self.line_buffer = self.line_buffer[newline_index + 1:]
            
            # Return it
            return line
        
def parse_mag(stream):
    """
    Yield MAG records from the given (possibly compressed) FASTQ stream.
    
    """
    
    for fastq_record in Bio.SeqIO.parse(stream, "fastq"):
        # For each FASTQ record
        
        # Parse the header line
        parts = fastq_record.description.split("\t")
        
        yield parts
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    for record in parse_mag(TransparentUnzip(options.in_file)):
        print(record)
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

