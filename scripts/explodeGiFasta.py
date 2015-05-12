#!/usr/bin/env python2.7
"""
explodeGiFasta.py: take a FASTA containing several records with proper Entrez gi
headers and turn it into a directory of single-record FASTAs.

The first record in the FASTA will be renamed to "ref" and saved as ref.fa in
the target directory, while the remaining records will be renamed to
GI<whatever> according to their GI numbers and saved as GI<whatever>.fa in the
target directory.

This is intended to be used to process the BRCA1/2 gene FASTAs generated from
<https://github.com/nouyang-curoverse/GA4GH_regions> so that they can be used
with createIndex.

"""

import argparse, sys, os, os.path, doctest, logging, random, subprocess, shutil
import itertools, collections

from Bio import AlignIO, SeqIO, Align
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
    parser.add_argument("fasta",
        help="FASTA file to explode")
    parser.add_argument("directory",
        help="directory to place individual FASTAs in")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
    
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
    
    if not os.path.isdir(options.directory):
        # If the output directory isn't good, set it up.
        os.makedirs(options.directory)
    
    # We treat the first record specially
    is_first = True
    
    for record in SeqIO.parse(options.fasta, "fasta"):
        # For every record in the file
        if is_first:
            # It's the first, so call it "ref"
            new_name = "ref"
            is_first = False
        else:
            # We need to name it based on its GI number. The IDs are like
            # "gi|568815581:43044294-43125482" to start with.
            gi_number = int(record.id.split("|")[1].split(":")[0])
            new_name = "GI{}".format(gi_number)
        
        # Set the record ID
        record.id = new_name
        
        # Save the record
        SeqIO.write(record, "{}/{}.fa".format(options.directory, record.id),
            "fasta")
        
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
