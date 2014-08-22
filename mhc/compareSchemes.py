#!/usr/bin/env python2.7
"""
compareSchemes.py: align pairs of MHC haplotypes with different LR merging
schemes (exact, inexact, credit, no credit), and compare the resulting
alignments to the official GRC alt loci alignments.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
from xml.etree import ElementTree
import tsv

import jobTree.scriptTree.target
import jobTree.scriptTree.stack
import sonLib.bioio

from targets import *

class SchemeAssessmentTarget(jobTree.scriptTree.target.Target):
    """
    A target that builds several reference structures by different schemes,
    collates their coverage results, and compares their MAF
    alignments to the ground truth.
    
    """
    
    def __init__(self, fasta_list, true_maf, seed, coverage_basename,
        agreement_basename, spectrum_basename, truth_basename):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, and comparing against the given truth MAF, using the
        specified RNG seed, and writing coverage statistics to files with
        specified coverage base name, a set of adjacency component size spectra
        to files with the specified spectrum basename, stats for agreement
        between the schemes to files anmed after the agreement basename, and a
        comparison against the true MAF to files named after the truth_basename.
        
        The coverage statistics are just alignment coverage for each pair of
        genomes, in <genome number>\t<coverage fraction> TSVs per scheme.
        
        The alignment statistics are just <precision>\t<recall> TSVs per scheme.
        
        true_maf and truth_filename may be None, but if one is None then both
        need to be None.
        
        Runs each subsequent FASTA against the first in each scheme.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(SchemeAssessmentTarget, self).__init__(memory=2147483648)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the random seed
        self.seed = seed
        
        # Save the concatenated coverage stats file name to use
        self.coverage_basename = coverage_basename
        
        # And the agreement file name
        self.agreement_basename = agreement_basename
        
        # And the concatenated frequency spectrum basename.
        # TODO: reduce or sum instead of concatenating?
        # TODO: just use an output directory and dump in it.
        self.spectrum_basename = spectrum_basename
        
        # Save the filename of a MAF to compare all our MAFs against (or None)
        self.true_maf = true_maf
        
        # And the filename to send the results of that comparison to (which also
        # may be None)
        self.truth_basename = truth_basename
        
        self.logToMaster(
            "Creating SchemeAssessmentTarget with seed {}".format(seed))
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting SchemeAssessmentTarget")
        
        # We need a few different follow-on jobs.
        followOns = []    
        
        # We'll keep track of our random state so we can seed without messing up
        # temp filenames.
        random_state = None
        
        # This will hold lists of alignments, in FASTA pair order, by scheme
        # name.
        alignments_by_scheme = {}
        
        for mismatch, credit in itertools.product([False, True], repeat=2):
            # Decide if we want mismatches and credit for this run.
            
            self.logToMaster("Preparing for mismatch: {} credit: {}".format(
                mismatch, credit))

            # Prepare the extra args that we want to send to createIndex to tell
            # it to use this scheme. Start out just setting context length
            extra_args = ["--context", "100"]
            # And give it a name to stick on our output files
            scheme = "Exact"
            if mismatch:
                extra_args.append("--mismatch")
                extra_args.append("--mismatches")
                extra_args.append("1")
                scheme = "Inexact"
            if credit:
                extra_args.append("--credit")
                scheme += "Credit"
            
            # Grab the reference FASTA (first in the list)
            reference_fasta = self.fasta_list[0]
            
            # How many children will we have for this scheme? One per FASTA to
            # compare against the reference FASTA
            num_children = len(self.fasta_list) - 1
            
            if random_state is not None:
                # Pull out our unseeded random state so filenames on successive
                # iterations don't conflict.
                random.setstate(random_state)
            
            # Make a temp file for each of the children with this scheme to
            # write coverage stats to.
            stats_filenames = [sonLib.bioio.getTempFile(
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's MAF alignment
            maf_filenames = [sonLib.bioio.getTempFile(
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # Save these so we can compare them to other schemes later.
            alignments_by_scheme[scheme] = maf_filenames
                
            # And another one to hold each child's adjacency component size
            # spectrum
            spectrum_filenames = [sonLib.bioio.getTempFile(
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
            
            # Save the RNG state before clobbering it with the seed.
            random_state = random.getstate()
            
            # Seed the RNG after sonLib does whatever it wants with temp file
            # names. Each scheme will get runs with the same seeds, but that is
            # probably OK.
            random.seed(self.seed)
            
            for i, other_fasta in enumerate(self.fasta_list[1:]):
                # For each other FASTA to compare against
                
                # Pull out the files that this child should output.
                stats_filename = stats_filenames[i]
                maf_filename = maf_filenames[i]
                spectrum_filename = spectrum_filenames[i]
                
                
                # Make a child to produce those, giving it a seed. Make sure to
                # give it only two FASTAs, reference first, so that when it
                # shuffles the non-reference ones it doesn't do anything. Also
                # make sure to tell it to use the spectrum output file.
                self.addChildTarget(ReferenceStructureTarget(
                    [reference_fasta, other_fasta], random.getrandbits(256), 
                    stats_filename, maf_filename,
                    spectrum_filename=spectrum_filename, extra_args=extra_args))
        
        
                
            # Make a follow-on job to merge all the child coverage outputs and
            # produce our coverage output file for this scheme.
            followOns.append(ConcatenateTarget(stats_filenames, 
                self.coverage_basename + "." + scheme))
                
            # Make a follow-on job to merge all the child spectrum outputs and
            # produce our spectrum output file for this scheme.
            followOns.append(ConcatenateTarget(spectrum_filenames, 
                self.spectrum_basename + "." + scheme))
            
            if self.true_maf is not None:
                # We also need another target for comparing all these MAFs
                # against the truth, which we have been given.
                followOns.append(AlignmentTruthComparisonTarget(self.true_maf, 
                maf_filenames, random.getrandbits(256), 
                self.truth_basename + "." + scheme))
                
        # Now we have all the followons for concatenating our stats and
        # comparing against the truth, we need a followon for comparing
        # corresponding MAFs from the same FASTA pair with different schemes,
        # for all combinations of schemes.
        followOns.append(AlignmentSchemeAgreementTarget(alignments_by_scheme, 
            random.getrandbits(256), self.agreement_basename))
            
            
        # So we need to run all of the follow-ons in parallel after this job
        self.setFollowOnTarget(RunTarget(followOns))
                
        self.logToMaster("SchemeAssessmentTarget Finished")

class AlignmentSchemeAgreementTarget(jobTree.scriptTree.target.Target):
    """
    A target that, given lists of corresponding alignments in a dict by scheme,
    compares corresponding alignments across all pairs of schemes and collates
    the results.
    
    """
    
    def __init__(self, alignments_by_scheme, seed, agreement_basename):
        """
        Takes a dict of lists of alignments by scheme (arranged so alignments at
        the same indices are of the same FASTAs), and a file basename to name
        all results after. Runs all pairwise agreement comparisons.
        
        Uses the seed to generate mafComparator seeds and temporary file names,
        so don't run two copies with the same seed.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AlignmentSchemeAgreementTarget, self).__init__(memory=2147483648)
        
        self.seed = seed
        
        # Save the alignment names
        self.alignments_by_scheme = alignments_by_scheme
        
        # And the agreement file name
        self.agreement_basename = agreement_basename
        
        self.logToMaster(
            "Creating AlignmentSchemeAgreementTarget with seed {}".format(seed))
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting AlignmentSchemeAgreementTarget")
        
        # Seed the RNG
        random.seed(self.seed)
        
        # We need a few different follow-on jobs.
        followOns = []    
        
        for scheme1, scheme2 in itertools.combinations(
            self.alignments_by_scheme.iterkeys(), 2):
            
            # For each pair of schemes to compare
            
            # What files will we concatenate? Get one for each pair of
            # corresponding alignments. Just use our seeded random state for
            # this; we only seeded once. So don't run this target twice with the
            # same seed.
            comparison_filenames = [sonLib.bioio.getTempFile(
                rootDir=self.getGlobalTempDir()) 
                for i in xrange(len(self.alignments_by_scheme[scheme1]))]
            
            for comparison_filename, alignment1, alignment2 in itertools.izip(
                comparison_filenames, self.alignments_by_scheme[scheme1], 
                self.alignments_by_scheme[scheme2]):
                
                # For corresponding alignments
                
                # Compare them
                self.addChildTarget(AlignmentComparisonTarget(alignment1, 
                    alignment2, random.getrandbits(256), comparison_filename))
                    
            # Concatenate all the comparisons into a file named after both
            # schemes.
            followOns.append(ConcatenateTarget(comparison_filenames,
                self.agreement_basename + "." + scheme1 + "." + scheme2))
        
        
        # We need to run all of the follow-ons in parallel after this job
        self.setFollowOnTarget(RunTarget(followOns))
                
        self.logToMaster("AlignmentSchemeAgreementTarget Finished")

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
    parser.add_argument("coverageBasename", 
        help="filename prefix to save the coverage stats output to")
    parser.add_argument("agreementBasename",
        help="filename prefix to save the alignment agreement stats output to")
    parser.add_argument("spectrumBasename", 
        help="filename prefix to save the adjacency component sizes to")
    parser.add_argument("fastas", nargs="+",
        help="FASTA files to index")
    parser.add_argument("--seed", default=random.getrandbits(256),
        help="seed for a particular deterministic run")
    parser.add_argument("--trueMaf", default=None,
        help="filename of a MAF to compare all generated MAFs against")
    parser.add_argument("--truthBasename", default=None,
        help="filename prefix to output comparison against the truth to")
        
    
    
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
    
    if (options.trueMaf is None) != (options.truthBasename is None):
        # Make sure they gave us a place to put it if they want a truth
        # comparison.
        raise Exception("--truthBasename and --trueMaf are corequisites.")
    
    # Make sure we've given everything an absolute module name.
    # Don't try to import * because that's illegal.
    if __name__ == "__main__":
        from compareSchemes import SchemeAssessmentTarget, \
            AlignmentSchemeAgreementTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(SchemeAssessmentTarget(
        options.fastas, options.trueMaf, options.seed, options.coverageBasename, 
        options.agreementBasename, options.spectrumBasename,
        options.truthBasename))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
