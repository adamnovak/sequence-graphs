#!/usr/bin/env python2.7
"""
evaluateMappings.py: evaluate a HAL alignment with a collection of evaluations
and metrics.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections, csv, tempfile

import doctest, pprint

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
    parser.add_argument("--beds", nargs="*",
        help=".bed file(s) of genes on the genomes in the HAL")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def parse_halstats_output(stream):
    """
    Parse the output of halStats on a HAL file, from the given iterable.
    
    Returns a list of dicts from halStats column header (GenomeName,
    NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments) to
    string or int value as appropriate.
    
    >>> lines = ["",
    ... "hal v2.1",
    ... "(GI568335879,ref)rootSeq;",
    ... "",
    ... "GenomeName, NumChildren, Length, NumSequences, NumTopSegments, "
    ... "NumBottomSegments",
    ... "rootSeq, 10, 5524595, 1, 0, 67141",
    ... "GI568335879, 0, 4672374, 1, 25067, 0",
    ... "ref, 0, 4970457, 1, 49220, 0",
    ... ""]
    >>> parsed = list(parse_halstats_output(lines))
    >>> len(parsed)
    3
    >>> parsed[1]["GenomeName"]
    'GI568335879'
    >>> parsed[0]["Length"]
    5524595

    """
    
    # Get a reader and skip the first 4 lines
    reader = csv.reader(stream)
    for i in xrange(4):
       next(reader)
        
    # Read the table header with the field names, which is the next nonempty
    # line.
    header = next(reader)
    
    # Throw out whitespace
    header = [x.strip() for x in header]
    
    for data in reader:
        if len(data) == 0:
            # We finished the table
            break
        
        elif len(data) != len(header):
            raise RuntimeError("Got {} fields for {} columns".format(len(data),
                len(header)))
                
        # Drop all the whitespace
        data = [x.strip() for x in data]
            
        # We'll fill in a dict with all the type-auto-detected values for each
        # column. We convert to int if a field is all digits (which is fine
        # since nothing will be negative).
        record = dict(zip(header, [int(x) if x.isdigit() else x for x in data]))
            
        # Yield each record as we parse it
        yield record
        
def metric_halstats(hal_filename):
    """
    A metric that runs halStats on the given HAL, and returns the results.
    
    The halStats program must be on the PATH.
    
    """
    
    # Build the halStats command line
    args = ["halStats", hal_filename]
    
    # Open the process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
        
    # Parse the results
    parsed = list(parse_halstats_output(process.stdout))
            
    if process.wait() != 0:
        raise RuntimeError("halStats failed")
            
    # We are done with this process.
    process.stdout.close()
    
    # Return the results
    return parsed
    
def hal2maf(hal_filename):
    """
    Convert the HAL at the given filename to a MAF in a temporary file. Returns
    the name of the temp file, which must be deleted by the caller.
    
    hal2maf must be installed and on the PATH.
    
    """
    
    # Make a temporary file
    handle, maf_filename = tempfile.mkstemp()
    os.close(handle)
    
    # Run the conversion
    subprocess.check_call(["hal2maf", hal_filename, maf_filename])
    
    # Return the filename of the MAF
    return maf_filename
    
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
    
    # Print the halStats output
    pprint.pprint(metric_halstats(options.hal))
    
    if options.beds is not None:

        # We want to do the checkGenes stuff, so we need a MAF
        maf_filename = hal2maf(options.hal)
        
        import checkGenes
        class_counts, gene_sets = checkGenes.check_genes(maf_filename,
            options.beds)
            
        # Clean up the MAF
        os.unlink(maf_filename)
        
        # Print the output
        pprint.pprint(class_counts)
        pprint.pprint(gene_sets)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
