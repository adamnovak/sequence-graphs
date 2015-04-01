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

class StructureAssessmentTarget(SchemeUsingTarget):
    """
    A target that builds several reference structures, collates their index
    vs. coverage results, and compares their MAF alignments.
    
    """
    
    def __init__(self, fasta_list, true_maf, seed, coverage_filename, 
        agreement_filename, truth_filename, num_children):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, and comparing against the given truth MAF, using the
        specified RNG seed, and writing coverage statistics to the specified
        file, and alignment agreement statistics to the other specicied file,
        and a comparison against the true MAF to the other other specified file.
        
        The coverage statistics are, specifically, alignment coverage of a
        genome vs. the order number at which the genome is added for all the
        child targets, and they are saved in a <genome number>\t<coverage
        fraction> TSV.
        
        The alignment statistics are just a <precision>\t<recall> TSV.
        
        TODO: This target is getting a bit unweildy and probably should be
        broken up to use an output directory and a configuration object or
        something.
        
        true_maf and truth_filename may be None, but if one is None then both
        need to be None.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(StructureAssessmentTarget, self).__init__(memory=2147483648)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the random seed
        self.seed = seed
        
        # Save the concatenated coverage stats file name to use
        self.coverage_filename = coverage_filename
        
        # And the alignment agreement stats filename
        self.agreement_filename = agreement_filename
        
        # Save the number of child targets to run
        self.num_children = num_children
        
        # Save the filename of a MAF to compare all our MAFs against (or None)
        self.true_maf = true_maf
        
        # And the filename to send the results of that comparison to (which also
        # may be None)
        self.truth_filename = truth_filename
        
        self.logToMaster(
            "Creating StructureAssessmentTarget with seed {}".format(seed))
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting StructureAssessmentTarget")
        
        # Make a temp file for each of the children to write coverage stats to
        stats_filenames = [sonLib.bioio.getTempFile(
            rootDir=self.getGlobalTempDir()) for i in xrange(self.num_children)]
            
        # And another one to hold each child's MAF alignment
        maf_filenames = [sonLib.bioio.getTempFile(
            rootDir=self.getGlobalTempDir()) for i in xrange(self.num_children)]
            
        # Seed the RNG after sonLib does whatever it wants with temp file names
        random.seed(self.seed)
        
        for stats_filename, maf_filename in zip(stats_filenames, maf_filenames):
            # Make a child to produce this reference structure, giving it a
            # 256-bit seed.
            self.addChildTarget(ReferenceStructureTarget(self.fasta_list,
                random.getrandbits(256), stats_filename, maf_filename, 
                extra_args=["--coverage", "100"]))
        
        # We need a few different follow-on jobs.
        followOns = []
                
        # Make a follow-on job to merge all the child coverage outputs and
        # produce our coverage output file.
        followOns.append(ConcatenateTarget(stats_filenames, 
            self.coverage_filename))
        
        # But we also need a follow-on job to analyze all those alignments
        # against each other.
        followOns.append(AlignmentSetComparisonTarget(maf_filenames, 
            random.getrandbits(256), self.agreement_filename))
            
        if self.true_maf is not None:
            # We also need another target for comparing all these MAFs against
            # the truth, which we have been given.
            followOns.append(AlignmentTruthComparisonTarget(self.true_maf, 
            maf_filenames, random.getrandbits(256), self.truth_filename))
            
            
        # So we need to run both of those in parallel after this job
        self.setFollowOnTarget(RunTarget(followOns))
                
        self.logToMaster("StructureAssessmentTarget Finished")

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
