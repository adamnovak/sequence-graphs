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

class RunTarget(jobTree.scriptTree.target.Target):
    """
    A target that runs other targets as its children.
    
    """
    
    def __init__(self, targets):
        """
        Make a new Target which, when run, adds all the targets in the given
        list as children.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(RunTarget, self).__init__(memory=2147483648)
        
        # Save the list of targets
        self.targets = targets
        
        
    def run(self):
        """
        Send off all the child targets we were given on construction.
        """
        
        self.logToMaster("Starting RunTarget")
        
        for target in self.targets:
            # Make each child a child.
            self.addChildTarget(target)
                
        self.logToMaster("RunTarget Finished")

class ReferenceStructureTarget(jobTree.scriptTree.target.Target):
    """
    A target that builds a reference structure.
    
    """
    
    def __init__(self, fasta_list, seed, coverage_filename, alignment_filename):
        """
        Make a new Target for building a reference structure from the given
        FASTAs, using the specified RNG seed, and writing coverage statistics to
        the specified file, and the alignment MAF to the other specified file.
        
        Each FASTA file must contain exactly one contig, named after the FASTA
        without the extension.
        
        Those coverage statistics are, specifically, alignment coverage of a
        genome vs. the order number at which the genome is added, and they are
        saved in a <genome number>\t<coverage fraction> TSV.
        
        """
        
        # Make the base Target. Ask for 20gb of memory since this is kinda hard.
        # Also ask for 4 CPUs because we can probably use as many as we can get.
        super(ReferenceStructureTarget, self).__init__(memory=(1024 ** 3) * 20,
            cpu=4)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the random seed
        self.seed = seed
        
        # Save the coverage file name to use
        self.coverage_filename = coverage_filename
        
        # And the alignemnt filename to use
        self.alignment_filename = alignment_filename
        
        self.logToMaster(
            "Creating ReferenceStructureTarget with seed {}".format(seed))
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting ReferenceStructureTarget")
        self.logToMaster("Running in {}".format(os.getcwd()))
        
        # Seed the RNG after sonLib does whatever it wants with temp file names
        random.seed(self.seed)
            
        # Generate an order for the FASTAs
        random.shuffle(self.fasta_list)
        
        # Make a temp directory for the index
        index_dir = sonLib.bioio.getTempFile(rootDir=self.getLocalTempDir())
        
        # Make a file for the cactus2hal
        c2h_filename = sonLib.bioio.getTempFile(rootDir=self.getLocalTempDir())
        
        # Make a file for the FASTA that we need to use the cactus2hal
        fasta_filename = sonLib.bioio.getTempFile(
            rootDir=self.getLocalTempDir())
        
        # Put together the arguments to invoke
        args = ["../createIndex/createIndex", "--context", 
            "100", "--scheme", "greedy", "--alignment", c2h_filename, 
            "--alignmentFasta", fasta_filename, index_dir] + self.fasta_list
            
        # Announce our command we're going to run
        self.logToMaster("Invoking {}".format(" ".join(args)))
        
        # Make an index with the FASTAs in that order, so we can read all the
        # logging output.
        process = subprocess.Popen(args,
            stdout=subprocess.PIPE)

        # Make a writer for the statistics
        writer = tsv.TsvWriter(open(self.coverage_filename, "w"))
    
        # Track how many lines of indexer output we process
        lines = 0
    
        for line in process.stdout:            
            # Collect and parse the log output, and get the coverage vs. genome
            # number data.
            
            if "WARNING" in line or "ERROR" in line:
                # Log things that came out at high logging priorities (TODO: or
                # are false positives) to our own log.
                print line
            
            # Record that we processed a line
            lines += 1
            
            if "Coverage from alignment of genome" in line:
                # Grab the genome number
                genome = int(line.split(":")[-2].split("genome")[1])
                # Grab the coverage fraction
                coverage = float(line.split("=")[1])
            
                # Save the coverage vs. genome number data as a reasonable TSV.
                writer.line(genome, coverage)
                
        # Close up the output file. 
        writer.close()
        
        # Clean up temporary index
        shutil.rmtree(index_dir)
        
        # Say how much info we got out of the indexer.
        self.logToMaster("Processed {} lines of indexer output".format(lines))
        
        if(process.wait() != 0):
            # If the indexing process died, complain.
            raise Exception("createIndex failed with code {}".format(
                process.returncode))
                
        # Now build the HAL file
        
        # What tree should we use? Assume the genome names are the same as the
        # FASTA names without their extensions. TODO: This tightly couples us to
        # both the genome name choosing logic and the particular test case we
        # are working on. Fix that.
        tree = "(" + ",".join([os.path.splitext(os.path.split(fasta)[1])[0] 
            for fasta in self.fasta_list]) + ")rootSeq;"
        
        # Where should we save it?
        hal_filename = sonLib.bioio.getTempFile(rootDir=self.getLocalTempDir())
        
        self.logToMaster("Creating HAL with tree: {}".format(tree))
        
        # Make sure HAL doesn't exist when we try to make it
        os.unlink(hal_filename)
        
        # Turn the c2h into a HAL with the given tree.
        subprocess.check_call(["halAppendCactusSubtree", c2h_filename, 
            fasta_filename, tree, hal_filename])
            
        # Clean up the temp files we needed to make the HAL
        os.unlink(c2h_filename)
        os.unlink(fasta_filename)
            
        self.logToMaster("Creating MAF")
            
        # Now we need to make that HAL into a MAF, saving to the file where
        # we're supposed to send out output
        subprocess.check_call(["hal2maf", hal_filename,
            self.alignment_filename])
            
        # Get rid of the intermediate HAL
        os.unlink(hal_filename)
        
        self.logToMaster("ReferenceStructureTarget Finished")
        
class StructureAssessmentTarget(jobTree.scriptTree.target.Target):
    """
    A target that builds several reference structures, collates their index
    vs. coverage results, and compares their MAF alignments.
    
    """
    
    def __init__(self, fasta_list, seed, coverage_filename, agreement_filename,
        num_children):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, using the specified RNG seed, and writing coverage
        statistics to the specified file, and alignment agreement statistics to
        the other specicied file.
        
        The coverage statistics are, specifically, alignment coverage of a
        genome vs. the order number at which the genome is added for all the
        child targets, and they are saved in a <genome number>\t<coverage
        fraction> TSV.
        
        The alignment styatistics are just a bunch of concatenated mafComparator
        XML files.
        
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
                random.getrandbits(256), stats_filename, maf_filename))
                
        # Make a follow-on job to merge all the child coverage outputs and
        # produce our coverage output file.
        coverageFollowOn = ConcatenateTarget(stats_filenames, 
            self.coverage_filename)
        
        # But we also need a follow-on job to analyze all those alignments.
        alignmentFollowOn = AlignmentSetComparisonTarget(maf_filenames, 
            random.getrandbits(256), self.agreement_filename)
            
        # So we need to run both of those in parallel after this job
        self.setFollowOnTarget(RunTarget([coverageFollowOn, alignmentFollowOn]))
                
        self.logToMaster("StructureAssessmentTarget Finished")
        
class AlignmentComparisonTarget(jobTree.scriptTree.target.Target):
    """
    A target that compares two MAF alignments.
    
    """
    
    def __init__(self, maf_a, maf_b, seed, output_filename):
        """
        Compare the two MAFs referred to by the given input filenames, and write
        a digest of mafComparator results to the given output filename. Uses a
        random seed to generate the mafComparator seed.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AlignmentComparisonTarget, self).__init__(memory=2147483648)
        
        # Save the parameters
        self.maf_a = maf_a
        self.maf_b = maf_b
        self.seed = seed
        self.output_filename = output_filename
        
        self.logToMaster("Creating AlignmentComparisonTarget")
        
        
    def run(self):
        """
        Run mafComparator and grab its results.
        """
        
        self.logToMaster("Starting AlignmentComparisonTarget")
        
        # Make a temp file for the XML output
        xml_filename = sonLib.bioio.getTempFile(rootDir=self.getLocalTempDir())
        
        # Seed the RNG after sonLib does whatever it wants with temp file names
        random.seed(self.seed)
        
        # Generate a seed for mafComparator that's a C-ish integer (not 256
        # bits)
        seed = random.getrandbits(32)
        
        # Run mafComparator and generate the output XML
        subprocess.check_call(["mafComparator", "--maf1", self.maf_a, "--maf2",
            self.maf_b, "--out", xml_filename, "--seed", str(seed)])
        
        # Now parse the XML output down to the statistics we actually want.
        tree = ElementTree.parse(xml_filename)
        
        # Grab the nodes with the aggregate results for each direction
        stats = tree.findall(
            "./homologyTests/aggregateResults/all")
            
        # Grab the averages
        averages = (stat.attrib["average"] for stat in stats)
        
        # Open the output file to save them to.
        writer = tsv.TsvWriter(open(self.output_filename, "w"))
        
        for average in averages:
            # Dump averages as a 1-column TSV.
            # TODO: Get more stuff to justify TSV-nes
            writer.line(average)
            
        writer.close()
        
        # Clean up our raw XML
        os.unlink(xml_filename)
        
        self.logToMaster("AlignmentComparisonTarget Finished")

class AlignmentSetComparisonTarget(jobTree.scriptTree.target.Target):
    """
    A target that compares a set of MAF alignments.
    
    """
    
    def __init__(self, maf_filenames, seed, output_filename):
        """
        Compare the MAFs referred to by the given input filenames, and write
        collated results to the given output filename. Uses a random seed
        to generate mafComparator seeds.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AlignmentSetComparisonTarget, self).__init__(memory=2147483648)
        
        # Save the parameters
        self.maf_filenames = maf_filenames
        self.seed = seed
        self.output_filename = output_filename
        
        self.logToMaster("Creating AlignmentSetComparisonTarget")
        
        
    def run(self):
        """
        Send off all the child targets to do alignment comparisons, and
        integrate their results.
        
        """
        
        self.logToMaster("Starting AlignmentSetComparisonTarget")
        
        # What pairs do we compare? TODO: This could get kinda big
        pairs = list(itertools.combinations(self.maf_filenames, 2))
        
        # Make a temp file for each of the children to write stats to
        comparison_filenames = [sonLib.bioio.getTempFile(
            rootDir=self.getGlobalTempDir()) for p in pairs]
        
        # Seed the RNG after sonLib does whatever it wants with temp file names
        random.seed(self.seed)
        
        for pair, comparison_filename in zip(pairs, comparison_filenames):
            # Make a child to compare these two MAF files and produce this
            # output file, giving it a 256-bit seed.
            self.addChildTarget(AlignmentComparisonTarget(pair[0], pair[1], 
                random.getrandbits(256), comparison_filename))
            
        # When we're done, concatenate all those comparison files together into
        # our output file.
        self.setFollowOnTarget(ConcatenateTarget(comparison_filenames,
            self.output_filename))
        
        self.logToMaster("AlignmentSetComparisonTarget Finished")
        
class ConcatenateTarget(jobTree.scriptTree.target.Target):
    """
    A target that concatenates several files together.
    
    """
    
    def __init__(self, input_filenames, output_filename):
        """
        Concatenate all the input files into the given output file. Deletes the
        input files.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(ConcatenateTarget, self).__init__(memory=2147483648)
        
        # Save the parameters
        self.input_filenames = input_filenames
        self.output_filename = output_filename
        
        self.logToMaster("Creating ConcatenateTarget")
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting ConcatenateTarget")
        
        output = open(self.output_filename, "w")
        
        for input_filename in self.input_filenames:
            # For each input file
            for line in open(input_filename):
                # Copy each line of the file to the output.
                output.write(line)
            
        output.close()
        
        self.logToMaster("ConcatenateTarget Finished")

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
        from compareOrders import ReferenceStructureTarget, \
            StructureAssessmentTarget, ConcatenateTarget, \
            AlignmentSetComparisonTarget, AlignmentComparisonTarget, RunTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(StructureAssessmentTarget(
        options.fastas, options.seed, options.coverageFile, 
        options.agreementFile, options.samples))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
