#!/usr/bin/env python2.7
"""
compareGraphs.py: align sets of haplotypes from different regions in different
orders, and compare the resulting graphs against each other and against the GRC
alignments for the regions (if present).

Takes a list of region directories, with the following structure:

<region>
<region>/ref.fa (optional, used first)
<region>/<genome>.fa
<region>/genes (optional, used for gene ortholog/paralog stats)
<region>/genes/<genome>/<whatever>.bed
<region>/<something>.maf (optional, used as true alignment)


"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest

import jobTree.scriptTree.target
import jobTree.scriptTree.stack
import sonLib.bioio

from targets import *

class GraphGenerationTarget(SchemeUsingTarget):
    """
    A target that decides what graphs to build on a given region, and builds the
    graphs.
    
    """
    
    def __init__(self, region_dir, output_dir, seed=0):
        """
        Make a new target for running a number of mapping/merging schemes
        against the FASTAs in the given directory.
        
        Output will be placed in the given output directory.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(GraphGenerationTarget, self).__init__(memory=2147483648)
        
        # Save arguments
        self.region_dir = region_dir
        self.output_dir = output_dir
        self.seed = seed
        
        # Make sure we have an RNG
        self.rng = random.Random()
        
        self.logToMaster("Creating {}".format(self.__class__.__name__))
        
   
    def getSchemePlan(self):
        """
        Return a collection of tuples describing schemes.
        
        """
        
        # Plan out all the schemes. They all start with map_type. For "natural"
        # schemes, we have map_type, min context, credit, mismatches for credit,
        # min edit distance bound, max edit distance. For "zip" schemes we have
        # map_type, min context, min edit distance, max range count, credit
        # mismatches.
        return set([
            # Do zip with min edit distance
            ("zip", 20, 1, None, 100, None),
            ("zip", 20, 2, None, 100, None),
            ("zip", 20, 3, None, 100, None),
            ("zip", 20, 4, None, 100, None),
            ("zip", 20, 5, None, 100, None),
            # And with 2-mismatch credit
            ("zip", 20, 1, None, 100, 2),
            ("zip", 20, 2, None, 100, 2),
            ("zip", 20, 3, None, 100, 2),
            ("zip", 20, 4, None, 100, 2),
            ("zip", 20, 5, None, 100, 2),
            # And with 1-mismatch tolerance and credit
            ("zip", 20, 1, 1, 100, 2),
            ("zip", 20, 2, 1, 100, 2),
            ("zip", 20, 3, 1, 100, 2),
            ("zip", 20, 4, 1, 100, 2),
            ("zip", 20, 5, 1, 100, 2)
        ])
        
    def run(self):
        """
        Set up all the child targets needed to actually run the schemes.
        
        """
        
        self.logToMaster("Starting {}".format(self.__class__.__name__))
        
        # Seed the RNG.
        self.rng.seed(self.seed)
        
        if not os.path.exists(self.output_dir):
            # Make our out directory exist
            os.makedirs(self.output_dir)
            
        # Get the list of FASTAs
        fastas = []
        for fasta in glob.glob(self.region_dir + "/*.fa"):
            # We'll stick in each FASTA one at a time.
            
            if os.path.basename(fasta) == "ref.fa":
                # The reference has to be first
                fastas = [fasta] + fastas
            else:
                # Otherwise just stick it on the end.
                fastas.append(fasta)
                
        # Look for true MAFs
        true_mafs = glob.glob(self.region_dir + "/*.maf")
        # Save the one we get, or None
        true_maf = None if len(true_mafs) == 0 else true_mafs[0]
        
        # Get all the gene beds. If there are no genes this will be empty.
        gene_beds = glob.glob(self.region_dir + "/genes/*/*.bed")
            
        for scheme, extra_args in self.generateSchemes():
            # We need to try each parameter set.
            
            # We'll run a sequence of targets for each scheme in the region:
            # make the graph, and compute per-graph statistics.
            sequence = []
            
            # Pre-calculate some intermediate file names
            maf_filename = "{}/{}.maf".format(self.output_dir, scheme)
            hal_filename="{}/{}.hal".format(self.output_dir, scheme)
            
            # We're going to make a ReferenceStructure for each scheme. This
            # will generate HAL and MAF files, and an itnernal coverage
            # assessment.
            sequence.append(ReferenceStructureTarget(fastas,
                self.rng.getrandbits(256),
                "{}/coverage.{}".format(self.output_dir, scheme), 
                maf_filename,
                hal_filename=hal_filename,
                extra_args=extra_args))
                
            # We'll make a bunch of evaluations to run in parallel
            evaluations = []
                
            if true_maf is not None:
                # We'll do a comparison against the true MAF
                evaluations.append(AlignmentTruthComparisonTarget(true_maf, 
                    [maf_filename], self.rng.getrandbits(256), 
                    "{}/truth.{}".format(self.output_dir, scheme)))
                    
            if len(gene_beds) > 0:
                # We'll do a gene-based evaluation.
                evaluations.append(MafGeneCheckerTarget(maf_filename, 
                    gene_beds, "{}/checkgenes.{}".format(self.output_dir,
                    scheme), "{}/checkgenes-genesets.{}".format(self.output_dir,
                    scheme)))
                    
            # TODO: Add a halStats evaluation, copying from
            # evaluateAlignment.py.
            
            # TODO: Set this up to be able to work with Cactus.
                    
            if len(evaluations) > 0:
                # If we have any evaluations to do, do them in parallel after we
                # make the graph.
                sequence.append(RunTarget(evaluations))
            
            # Run the graph generation and any evaluations as a child.
            self.addChildTarget(SequenceTarget(sequence))
        
        self.logToMaster("{} finished".format(self.__class__.__name__))

class MultiRegionTarget(jobTree.scriptTree.target.Target):
    """
    A target that runs GraphGenerationTargets on several regions.
    
    """
    
    def __init__(self, region_dirs, output_dir, seed=0):
        """
        Make a new target for makign graphs for each of the given directories.
        
        Output will be placed in the given output directory.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(MultiRegionTarget, self).__init__(memory=2147483648)
        
        # Save arguments
        self.region_dirs = region_dirs
        self.output_dir = output_dir
        self.seed = seed
        
        # Make sure we have an RNG
        self.rng = random.Random()
        
        self.logToMaster("Creating {}".format(self.__class__.__name__))
        
   
    def run(self):
        """
        Set up all the child targets for the different regions.
        
        """
        
        self.logToMaster("Starting {}".format(self.__class__.__name__))
        
        # Seed the RNG.
        self.rng.seed(self.seed)
        
        if not os.path.exists(self.output_dir):
            # Make our out directory exist
            os.makedirs(self.output_dir)
            
        for region_dir in self.region_dirs:
            # For every region, get just its name
            region_name = os.path.basename(region_dir)
            
            # Make all the graph types for that region, with per-region output
            # dirs.
            self.addChildTarget(GraphGenerationTarget(region_dir, 
                self.output_dir + "/" + region_name,
                seed=self.rng.getrandbits(256)))
        
        self.logToMaster("{} finished".format(self.__class__.__name__))

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
    parser.add_argument("out_dir", 
        help="directory to fill with alignments, statistics files, and hubs")
    parser.add_argument("in_dirs", nargs="+",
        help="region directories to use")
    parser.add_argument("--seed", default=random.getrandbits(256),
        help="seed for a particular deterministic run")
    
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
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Make sure we've given everything an absolute module name.
    # Don't try to import * because that's illegal.
    if __name__ == "__main__":
        from compareGraphs import GraphGenerationTarget, MultiRegionTarget
        
    # Make a stack of jobs to run, starting with all our arguments.
    stack = jobTree.scriptTree.stack.Stack(MultiRegionTarget(options.in_dirs,
        options.out_dir, options.seed))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

