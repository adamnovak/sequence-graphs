#!/usr/bin/env python2.7
"""
shred.py: chop the sequences in a FASTA into fake reads.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    parser.add_argument("--fastaIn", type=argparse.FileType("r"), 
        default=sys.stdin, 
        help="FASTA file of sequence to process")
    parser.add_argument("--fastaOut", type=argparse.FileType("w"), 
        default=sys.stdout,
        help="FASTA file of reads to write")
    parser.add_argument("--size", type=int, default=200,
        help="read length")
    parser.add_argument("--spacing", type=int, default=100,
        help="spacing between read start points")
    
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def shred(records_in, size, spacing):
    """
    Given an iterator of SeqRecords, a read size, and a start position spacing,
    yield SeqRecords for the reads.
    """
    
    for input_record in records_in:
        # Get each record
        for start_pos in xrange(0, len(input_record), spacing):
            # Start at each start position
            
            if len(input_record) - start_pos < size:
                # Skip any too-small trailing reads
                continue
                
            # Subset out the part we want to use as a read
            read = input_record[start_pos:start_pos + size]
            # Add the coordinates to the read's FASTA ID
            read.id += ":{}-{}".format(start_pos, start_pos + size)
            read.description = ""
            
            # Give out the read we've made
            yield read
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    

    # Stream records in, shred them, and stream the results out.
    SeqIO.write(shred(SeqIO.parse(options.fastaIn, "fasta"), options.size, 
        options.spacing), options.fastaOut, "fasta")
            
    
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
