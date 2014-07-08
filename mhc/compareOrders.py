#!/usr/bin/env python2.7
"""
compareOrders.py: build reference structures from the MHC haplotypes in
different orders and compare statistics gathered about each.

"""

import argparse, sys, random, subprocess
import tsv

import jobTree.scriptTree.target
import jobTree.scriptTree.stack
import sonLib.bioio

class ReferenceStructureTarget(jobTree.scriptTree.target.Target):
    """
    A target that builds a reference structure.
    
    """
    
    def __init__(self, fasta_list, seed, output_filename):
        """
        Make a new Target for building a reference structure from the given
        FASTAs, using the specified RNG seed, and writing statistics to the
        specified file.
        
        Those statistics are, specifically, alignment coverage of a genome vs.
        the order number at which the genome is added, and they are saved in a
        <genome number>\t<coverage fraction> TSV.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(CompareManyTarget, self).__init__(memory=2147483648)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the random seed
        self.seed = seed
        
        # Save the output file name to use
        self.output_filename = output_filename
        
        self.logToMaster(
            "Creating ReferenceStructureTarget with seed {}".format(seed))
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting ReferenceStructureTarget")
        
        # Seed the RNG after sonLib does whatever it wants with temp file names
        random.seed(self.seed)
            
        # Generate an order for the FASTAs
        random.shuffle(self.fasta_list)
        
        # Make a temp directory for the index
        indexDir = sonLib.bioio.getTempFile(rootDir=self.getLocalTempDir())
        
        # Make an index with the FASTAs in that order, so we can read all the
        # logging output.
        process = subprocess.Popen(["../createIndex/createIndex", "--context", 
            "100", "--scheme", "greedy", indexDir] + self.fasta_list,
            stdout=subprocess.PIPE)

        # Make a writer for the statistics
        writer = tsv.TsvWriter(open(self.output_filename))
    
        for line in process.stdout:            
            # Collect and parse the log output, and get the coverage vs. genome
            # number data.
            
            if "Coverage from alignment of genome" in line:
                # Grab the genome number
                genome = int(line.split(":")[0].split("genome")[1])
                # Grab the coverage fraction
                coverage = float(line.split("=")[1])
            
                # Save the coverage vs. genome number data as a reasonable TSV.
                writer.line(genome, coverage)
               
        # Close up the output file. 
        writer.close()
                
        self.logToMaster("ReferenceStructureTarget Finished")

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
    parser.add_argument("outputFile", 
        help="filename to save the output to")
    parser.add_argument("fastas", nargs="+",
        help="FASTA files to index")
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
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Make sure we've given everything an absolute module name.
    # Don't try to import * because that's illegal.
    if __name__ == "__main__":
        from compareOrders import ReferenceStructureTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(ReferenceStructureTarget(
        options.fastas, options.seed, options.outputFile))
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
