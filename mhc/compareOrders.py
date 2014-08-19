#!/usr/bin/env python2.7
"""
compareOrders.py: build reference structures from the MHC haplotypes in
different orders and compare statistics gathered about each.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
from xml.etree import ElementTree
import tsv

import jobTree.scriptTree.target
import jobTree.scriptTree.stack
import sonLib.bioio

from targets import *

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
    parser.add_argument("coverageFile", 
        help="filename to save the coverage stats output to")
    parser.add_argument("agreementFile",
        help="filename to save the alignment agreement stats output to")
    parser.add_argument("fastas", nargs="+",
        help="FASTA files to index")
    parser.add_argument("--seed", default=random.getrandbits(256),
        help="seed for a particular deterministic run")
    parser.add_argument("--samples", type=int, default=1,
        help="number of samples to take")
    parser.add_argument("--trueMaf", default=None,
        help="filename of a MAF to compare all generated MAFs against")
    parser.add_argument("--truthFile", default=None,
        help="filename to output comparison against the truth to")
        
    
    
    # Add the jobTree options
    jobTree.scriptTree.stack.Stack.addJobTreeOptions(parser)
        
    
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
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    if options.trueMaf is not None and options.truthFile is None:
        # Make sure they gave us a place to put it if they want a truth
        # comparison.
        raise Exception("--truthFile is required if --trueMaf is used.")
    
    # Make sure we've given everything an absolute module name.
    # Don't try to import * because that's illegal.
    if __name__ == "__main__":
        from compareOrders import ReferenceStructureTarget, \
            StructureAssessmentTarget, ConcatenateTarget, \
            AlignmentSetComparisonTarget, AlignmentTruthComparisonTarget, \
            AlignmentComparisonTarget, RunTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(StructureAssessmentTarget(
        options.fastas, options.trueMaf, options.seed, options.coverageFile, 
        options.agreementFile, options.truthFile, options.samples))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
