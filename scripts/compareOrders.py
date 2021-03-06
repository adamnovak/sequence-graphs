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
    
    The comparison is made for several schemes.
    
    """
    
    def __init__(self, fasta_list, true_maf, seed, stats_dir, num_children):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, and comparing against the given truth MAF, using the
        specified RNG seed, and writing coverage statistics, alignment agreement
        statistics, and a comparison against the true MAF to files in the
        specified stats directory.
        
        The coverage statistics are, specifically, alignment coverage of a
        genome vs. the order number at which the genome is added for all the
        child targets, and they are saved in a <genome number>\t<coverage
        fraction> TSV.
        
        The alignment statistics are just a <precision>\t<recall> TSV.
        
        true_maf and truth_filename may be None, but if one is None then both
        need to be None.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(StructureAssessmentTarget, self).__init__(memory=2147483648)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the random seed
        self.seed = seed
        
        # Make sure we have an RNG
        self.rng = random.Random()
        
        # Save the stats directory
        self.stats_dir = stats_dir
        
        # Save the number of child targets to run
        self.num_children = num_children
        
        # Save the filename of a MAF to compare all our MAFs against (or None)
        self.true_maf = true_maf
        
        self.logToMaster(
            "Creating StructureAssessmentTarget with seed {}".format(seed))
            
    def getSchemePlan(self):
        """
        Return a collection of tuples describing schemes.
        
        """
        
        # Plan out all the schemes. They all start with map_type. For "natural"
        # schemes, we have map_type, min context, credit, mismatches for credit,
        # min edit distance bound, max edit distance. For "zip" schemes we have
        # map_type, min context, max range count.
        return set([
            # Do varying edit counts
            ("zip", 100, 0, 100),
            ("zip", 100, 5, 100),
            ("zip", 100, 10, 100),
            ("zip", 100, 15, 100),
            ("zip", 50, 0, 100),
            ("zip", 50, 5, 100),
            ("zip", 50, 10, 100),
            ("zip", 50, 15, 100)
        ])
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting StructureAssessmentTarget")
        
        if not os.path.exists(self.stats_dir):
            # Make our out directory exist
            os.makedirs(self.stats_dir)
          
        if not os.path.exists(self.stats_dir + "/hals"):
            # And our directory for HALs
            os.makedirs(self.stats_dir + "/hals")
            
        if not os.path.exists(self.stats_dir + "/mapStats"):
            # And our directory for mapping scheme stats
            os.makedirs(self.stats_dir + "/mapStats")
            
        # We need a few different follow-on jobs, and we need them for all the
        # schemes put together.
        followOns = []
            
        # What genomes are we using? The first is the reference.
        genomes = [os.path.splitext(os.path.basename(fasta))[0] for fasta in
            self.fasta_list]
        
        for scheme, extra_args in self.generateSchemes():
         
            # Each scheme gets its own output files.
         
            # Where should we save the coverage stats?
            coverage_filename = "{}/orderCoverage.{}".format(self.stats_dir,
                scheme)
            
            # Where should we put the agreement precision recall stats?
            agreement_filename = "{}/agreement.{}".format(self.stats_dir,
                scheme)
                
            # Where should we save the comparison to the truth's precision and
            # recall?
            truth_filename = "{}/truth.{}".format(self.stats_dir, scheme)
            
            # Where should we put the spectrum data for rearrangement sizes?
            spectrum_filename = "{}/spectrum.{}".format(self.stats_dir, scheme)
            
            # Where should we put the tandem duplication data?
            tandem_filename = "{}/tandem.{}".format(self.stats_dir, scheme)
            
            # Where should the reference coverage data go (coverage only vs. the
            # reference)
            reference_coverage_filename = "{}/refCoverage.{}".format(
                self.stats_dir, scheme)
                
            # And we want the runtimes of all the schemes in seconds.
            time_filename = "{}/runtime.{}".format(self.stats_dir, scheme)
            
            # Make a temp file for each of the children to write coverage stats
            # to
            stats_filenames = [sonLib.bioio.getTempFile(
                rootDir=self.getGlobalTempDir()) for i in xrange(
                self.num_children)]
                
            # And another one to hold each child's MAF alignment
            maf_filenames = [sonLib.bioio.getTempFile(
                rootDir=self.getGlobalTempDir()) for i in xrange(
                self.num_children)]
                
            # And we need the spectrum data from each run. 
            spectrum_filenames = [sonLib.bioio.getTempFile(suffix=".spectrum",
                rootDir=self.getGlobalTempDir()) for i in xrange(
                self.num_children)]
                
            # And we need the tandem duplication data from each run. 
            tandem_filenames = [sonLib.bioio.getTempFile(suffix=".tandem",
                rootDir=self.getGlobalTempDir()) for i in xrange(
                self.num_children)]
                
            # Each reference structure will record the number of seconds it took
            # to build in files in this list. TODO: convert the other filename
            # list/concatenate target patterns to work like this now that the
            # temp file names don't break the seeded RNG. TODO: Make
            # ConcatenateTarget spit out input file names.
            time_filenames = []
                
            # Seed the RNG so we get the same orders with different schemes.
            self.rng.seed(self.seed)
            
            # Keep track of which child is which within a scheme.
            child_number = 0
            
            for (stats_filename, maf_filename, spectrum_data_filename, 
                tandem_data_filename) in zip(stats_filenames, maf_filenames, 
                spectrum_filenames, tandem_filenames):
                
                # Where should we put the HAL file? We're going to make one for
                # every ordering we run, in every scheme.
                hal_filename = "{}/hals/{}.{}.hal".format(self.stats_dir,
                    scheme, child_number)
                
                # Where should we put the mapping scheme stats? We collect this
                # for every ordering and scheme too. TODO: Collate with a
                # summing target somehow across orderings.
                map_stats_filename = "{}/mapStats/{}.{}.tsv".format(
                    self.stats_dir, scheme, child_number)
                
                # Where should we record our runtime?
                time_data_filename = sonLib.bioio.getTempFile(suffix=".time",
                    rootDir=self.getGlobalTempDir())
                time_filenames.append(time_data_filename)    
                
                
                # Make a child to produce this reference structure, giving it a
                # 256-bit seed. TODO: Make it so we can set fields not all at
                # once in the constructor to wire up the data flow.
                self.addChildTarget(ReferenceStructureTarget(self.fasta_list,
                    self.rng.getrandbits(256), stats_filename, maf_filename,
                    spectrum_filename=spectrum_data_filename,
                    tandem_filename=tandem_data_filename,
                    hal_filename=hal_filename, time_filename=time_data_filename,
                    map_stats_filename=map_stats_filename,
                    extra_args=extra_args))
                    
                # Next child for this scheme gets a different number.
                child_number += 1
            
            # Make a follow-on job to merge all the child coverage outputs and
            # produce our coverage output file.
            followOns.append(ConcatenateTarget(stats_filenames, 
                coverage_filename))
                
            # Make a follow-on job to merge all the child coverage outputs and
            # produce our coverage output file.
            followOns.append(ConcatenateTarget(spectrum_filenames, 
                spectrum_filename))
            
            # But we also need a follow-on job to analyze all those alignments
            # against each other.
            followOns.append(AlignmentSetComparisonTarget(maf_filenames, 
                self.rng.getrandbits(256), agreement_filename))
                
            # And we have one to get the coverage against the reference for all
            # the genomes in all the MAFs
            followOns.append(AlignmentSetCoverageTarget(maf_filenames,
                genomes[0], genomes[1:], reference_coverage_filename))
                
            # Concatenate all the runtimes for this scheme together.
            followOns.append(ConcatenateTarget(time_filenames, time_filename))
                
            if self.true_maf is not None:
                # We also need another target for comparing all these MAFs
                # against the truth, which we have been given.
                followOns.append(AlignmentTruthComparisonTarget(self.true_maf, 
                    maf_filenames, self.rng.getrandbits(256), truth_filename))
            
            
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
    parser.add_argument("outDir", 
        help="directory to fill with statistics files and hubs")
    parser.add_argument("fastas", nargs="+",
        help="FASTA files to index")
    parser.add_argument("--seed", default=random.getrandbits(256),
        help="seed for a particular deterministic run")
    parser.add_argument("--samples", type=int, default=1,
        help="number of samples to take")
    parser.add_argument("--trueMaf", default=None,
        help="filename of a MAF to compare all generated MAFs against")
        
    
    
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
    
    if __name__ == "__main__":
        # Sure would be nice if Benedict could fix this.
        from compareOrders import StructureAssessmentTarget
    
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(StructureAssessmentTarget(
        options.fastas, options.trueMaf, options.seed, options.outDir,
        options.samples))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
