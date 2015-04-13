#!/usr/bin/env python2.7
"""
evaluateMappings.py: evaluate a HAL alignment with a collection of evaluations
and metrics.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections, csv

import tsv

from Bio import AlignIO, SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# We need yet another comprehensive bio library in Python since BioPython hasn't
# got interval trees.
from bx.intervals import IntervalTree



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
    parser.add_argument("hal",
        help="HAL file to evaluate")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def metric_halstats(hal_filename):
    """
    A metric that runs halStats on the given HAL, and returns the results.
    
    Returns a list of dicts from halStats column header (GenomeName,
    NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments) to
    string or int value as appropriate.
    
    The halStats program must be on the PATH.
    
    """
    
    # Build the halStats command line
    args = ["halStats", hal_filename]
    
    # Open the process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
        
    # Get a reader and skip the first 3 lines
    reader = csv.reader(process.stdout)
    for i in xrange(2):
        next(reader)
        
    # Read the table header
    header = next(reader)
    
    # This is the list of dicts we will populate
    records = []
        
    for data in csv.reader(process.stdout):
        if len(data) == 0:
            # We finished the table
            break
        
        elif len(data) != len(header):
            raise RuntimeError("Table width mismatch!")
            
        # We'll fill in a dict with all the type-auto-detected values for each
        # column. We convert to int if a field is all digits (which is fine
        # since nothing will be negative).
        record = dict(zip(header, [int(x) if x.isdigit() else x for x in data]))
            
        # Put it in the list
        records.append(record)
        
    if process.wait() != 0:
        raise RuntimeError("halStats failed")
            
    # We are done with this process.
    process.stdout.close()

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
                

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
