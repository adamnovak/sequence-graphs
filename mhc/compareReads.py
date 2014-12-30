#!/usr/bin/env python2.7
"""
compareReads.py: align fake reads from all the MHC alts to the reference with
BWA and with my new read aligning tool, and compare the quality of gene
mappings.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
from xml.etree import ElementTree
import tsv

import jobTree.scriptTree.target
import jobTree.scriptTree.stack
import sonLib.bioio

from targets import *

class TargetQueuer(object):
    """
    Class we can use to build complicated layouts of targets to execute.
    
    """
    
    def __init__(self, target):
        """
        Make a new TargetQueuer to queue up targets. No subtasks exist by
        default.
        
        Takes the parent target so it can use the temp directory.
        
        TODO: Just make this into a target?
        
        """
        
        # We have a stack of (mode, list of targets) tuples.
        self.stack = []
        
        # What target will we run to run everything?
        self.root = None
        
        # What target can we steal temp files from
        self.target = target
        
    def subtask_serial(self, name=None):
        """
        Start a new subtask of things that happen in serial. It may have a name,
        which isn't used.
        
        """
        
        self.stack.append(("serial", []))
    
    def subtask_parallel(self, name=None):
        """
        Start a new subtask of things that happen in parallel. It may have a
        name, which isn't used.
        
        """
        
        self.stack.append(("parallel", []))
        
    def end_subtask(self):
        """
        Finish the most recently added subtask, and add it as a target to its
        parent subtask.
        
        """
        
        mode, targets = self.stack.pop()
        
        if mode == "serial":
            # Make a target that runs all those tasks in serial.
            self.append(SequenceTarget(targets))
        elif mode == "parallel":
            # Make a target that runs all those tasks in parallel.
            self.append(RunTarget(targets))
        else:
            # We broke ourselves.
            raise Exception("Invalid mode: {}".format(mode))
            
    def append(self, target):
        """
        Add a new target to the current subtask.
        
        """
        
        if len(self.stack) > 0:
            # We have a subtask.
            # Get the top thing on the stack, find its list, and append to that.
            self.stack[-1][1].append(target)
        else:
            # No subtask active. We have finished the last one. Make it the
            # root. TODO: This should never be called except through
            # end_subtask.
            
            if self.root is not None:
                raise Exception("Trying to set root twice!")
            
            self.root = target
            
    def get_root(self):
        """
        Get the root task that runs all the subtasks. Might be None.
        """
        
        if self.root is None:
            # We don't have a root.
            raise Exception(
                "Did not pop final subtask, or no subtasks were made")
        
        return self.root
        
    def tempfile(self):
        """
        Quickly and easily grab a temp file name.
        """
        
        return sonLib.bioio.getTempFile(rootDir=self.target.getGlobalTempDir())
            
            
        

class AlignerAssessmentTarget(jobTree.scriptTree.target.Target):
    """
    A target that assesses each of seevral schemes along with BWA.
    
    """
    
    def __init__(self, fasta_list, gene_beds, out_dir):
        """
        Make a target that assesses each of several schemes along with BWA.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AlignerAssessmentTarget, self).__init__(memory=2147483648)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the genes we want to compare against
        self.gene_beds = gene_beds
        
        # Make output in the output directory
        self.out_dir = out_dir
        
        self.logToMaster("Creating AlignerAssessmentTarget")
        
   
    def getSchemePlan(self):
        """
        Return a list of tuples describing schemes.
        
        """
        
        # Plan out all the schemes as mismatch, credit, unstable, min_context,
        # add_context, mult_context, ignore_below, hamming_bound, hamming_max,
        # map_type
        return set([
            # Exact credit (tolerating 1 mismatch) with Hamming bound, but no
            # mismatches.
            (True, True, False, None, None, None, None, 1, None, "old"),
            # Natural,  Hamming bound 6, 5 mismatches in gaps.
            (True, True, False, None, None, None, None, 6, 5, "old"),
            # Natural, Hamming bound 3, 2 mismatches.
            (True, True, False, None, None, None, None, 3, 2, "old"), 
            # Natural Hamming bound 1
            (True, True, False, None, None, None, None, 1, None, "old"),
            # Unstable versions
            (True, True, True, None, None, None, None, 1, None, "old"),
            (True, True, True, None, None, None, None, 6, 5, "old"),
            (True, True, True, None, None, None, None, 3, 2, "old"),
            (True, True, True, None, None, None, None, 1, None, "old"),
            # Unstable new versions
            (True, True, True, None, None, None, None, 1, None, "natural"),
            (True, True, True, None, None, None, None, 6, 5, "natural"),
            (True, True, True, None, None, None, None, 3, 2, "natural"),
            (True, True, True, None, None, None, None, 1, None, "natural")
        ])

        
   
    def generateSchemes(self):
        """
        Yield tuples of scheme names and the extra arguments mapReads needs to
        use those mapping schemes.
        
        """
        
        # Get a big set of tuples succinctly describing all the schemes we want
        # to run.
        scheme_plan = self.getSchemePlan()
        
        for mismatch, credit, unstable, min_context, add_context, \
            mult_context, ignore_below, hamming_bound, hamming_max, map_type \
            in scheme_plan:
            # Unpack each planned scheme
            
            # Start out with no configuration arguments
            extra_args = []
            
            # Give the scheme a name to stick on our output files
            scheme_name = "E"
            if mismatch:
                # Add the args and scheme name component for mismatch
                extra_args.append("--mismatches")
                extra_args.append("1")
                scheme_name = "I"
            if credit:
                # Add the args and scheme name component for credit
                extra_args.append("--credit")
                scheme_name += "C"
            else:
                # Add an N for no credit
                scheme_name += "N"
            
            # Handle mapping types
            scheme_name += map_type
            extra_args.append("--mapType")
            extra_args.append(map_type)
            
            if min_context is not None:
                # Require a min context
                extra_args.append("--context")
                extra_args.append(str(min_context))
                
                # Include min conext in the name if in use
                scheme_name += "Min{}".format(min_context)
                
            if add_context is not None:
                # Require an additional context
                extra_args.append("--addContext")
                extra_args.append(str(add_context))
                
                # Include additional context if in use. No need for padding
                # since we can now natural sort.
                scheme_name += "Add{}".format(add_context)
                
            if mult_context is not None:
                # Require a context multiplier
                extra_args.append("--multContext")
                extra_args.append(str(mult_context))
                
                # Include multiplicative context if in use. Don't let any .s
                # into the scheme name since it needs to be a valid file
                # extension.
                scheme_name += "Mult{}".format(mult_context).replace(".", "p")
                
            if ignore_below is not None:
                # For the natrual mapping scheme, require a minimum length to
                # even count a maximum unique match (so short by-chance ones
                # can't cause conflicts that blacklist bases.)
                extra_args.append("--ignoreMatchesBelow")
                extra_args.append(str(ignore_below))
                
                # Mention it in the scheme name
                scheme_name += "Ign{}".format(ignore_below)
                
            if hamming_bound is not None:
                # For the natural mapping scheme, don't map on a maximum unique
                # match run unless we can get a lower bound on its Hamming
                # distance from all other reference locations that is at least
                # this high.
                extra_args.append("--minHammingBound")
                extra_args.append(str(hamming_bound))
                
                # Mention it in the scheme name
                scheme_name += "Ham{}".format(hamming_bound)
                
            if hamming_max is not None:
                # For the natural mapping scheme, allow this many mismatches in
                # maximal unique match runs.
                extra_args.append("--maxHammingDistance")
                extra_args.append(str(hamming_max))
                
                # Mention it in the scheme name
                scheme_name += "Mis{}".format(hamming_max)
                
            if unstable:
                # Allow some instability
                extra_args.append("--unstable")
                scheme_name += "U"
                
            # Yield the name with the args.
            yield (scheme_name, extra_args)
        
        
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting AlignerAssessmentTarget")
        
        if not os.path.exists(self.out_dir):
            # Make sure our output directory is real.
            os.makedirs(self.out_dir)
        
        # Make a TargetQueuer so we can phrase this top-down control logic in
        # terms of subtasks.
        queuer = TargetQueuer(self)
        
        # Start out saying things to do in serial.
        queuer.subtask_serial("shredThenMap")
        
        # Make a file for each query's reads
        read_files = [queuer.tempfile() for _ in self.fasta_list[1:]]
           
        # Shred all the reads in parallel
        queuer.subtask_parallel("mapInParallel")
        for contig, reads in itertools.izip(self.fasta_list[1:], read_files):
            queuer.append(ShredTarget(contig, reads))
        queuer.end_subtask()
        
        # Map reads with each scheme in parallel
        queuer.subtask_parallel("mapAllSchemes")
        
        # And we need a serial subtask for BWA: do all the mappings, then merge
        # them.
        queuer.subtask_serial("mapConcatBWA")
        
        # Do all the mappings and check all the genes.
        queuer.subtask_parallel("mapAllGenomes")
            
        # This holds the checkGenes class counts files from each genome mapped
        class_count_files = []
        
        # This holds the gene set files
        gene_set_files = []
            
        for reads in read_files:
            # Make a file to save the BWA mappings from this set of reads in
            bwa_mappings = queuer.tempfile()
                
            # We need to BWA map and check genes in serial.
            queuer.subtask_serial("mapThenCheck")
                
            # Make a target to do the mapping
            queuer.append(BWATarget(self.fasta_list[0], reads, bwa_mappings))
            
            # Make a file to put the class counts in
            class_count_file = queuer.tempfile()
            class_count_files.append(class_count_file)
            
            # Make a file to put the gene sets in
            gene_set_file = queuer.tempfile()
            gene_set_files.append(gene_set_file)
            
            # And one to do the gene checking
            queuer.append(TsvGeneCheckerTarget(bwa_mappings, self.gene_beds, 
                class_count_file, gene_set_file))
                
            queuer.end_subtask()
            
        # We finished all the mapping and gene checking
        queuer.end_subtask()
        
        # Now we need to merge the checkGenes output TSVs
        queuer.subtask_parallel("concat")
        # Concatenate the class counts
        queuer.append(ConcatenateTarget(class_count_files, 
            self.out_dir + "/BWA.counts"))
        # And the gene set files
        queuer.append(ConcatenateTarget(gene_set_files, 
            self.out_dir + "/BWA.genes"))
        # Output TSVs are merged
        queuer.end_subtask()
            
        # End the subtask of the serial steps for the BWA scheme
        queuer.end_subtask()
        
        for scheme_name, extra_args in self.generateSchemes():
            # Now we need to do all the other schemes. Each will have a mapping
            # step and a concatenate step in order.
            queuer.subtask_serial("mapConcatScheme")
            
            # Do all the mappings and check all the genes.
            queuer.subtask_parallel("mapAllScheme")
                
            # This holds the checkGenes class counts files from each genome
            # mapped
            class_count_files = []
            
            # This holds the gene set files
            gene_set_files = []
                
            for reads in read_files:
                # Make a file to save the mapReads mappings from this set of
                # reads
                mappings = queuer.tempfile()
                    
                # We need to map and check genes in serial for each genome.
                queuer.subtask_serial("mapThenCheck")
                    
                # Make a target to do the mapping
                queuer.append(MapReadsTarget(self.fasta_list[0], reads,
                    mappings, extra_args=extra_args))
                
                # Make a file to put the class counts in
                class_count_file = queuer.tempfile()
                class_count_files.append(class_count_file)
                
                # Make a file to put the gene sets in
                gene_set_file = queuer.tempfile()
                gene_set_files.append(gene_set_file)
                
                # And one to do the gene checking
                queuer.append(TsvGeneCheckerTarget(mappings, self.gene_beds, 
                    class_count_file, gene_set_file))
                    
                queuer.end_subtask()
                
            # We finished all the mapping and gene checking
            queuer.end_subtask()
            
            # Now we need to merge the checkGenes output TSVs
            queuer.subtask_parallel("concat")
            # Concatenate the class counts
            queuer.append(ConcatenateTarget(class_count_files, 
                self.out_dir + "/" + scheme_name + ".counts"))
            # And the gene set files
            queuer.append(ConcatenateTarget(gene_set_files, 
                self.out_dir + "/" + scheme_name + ".genes"))
            # Output TSVs are merged
            queuer.end_subtask()
            
            # And that was the last serial step for this scheme.
            queuer.end_subtask()
            
        # End the subtask of doing all the schemes in parallel
        queuer.end_subtask()
        
        # And end the shred-then-map process.
        queuer.end_subtask()
        
        # Run all the tasks we put in the queuer.
        self.addChildTarget(queuer.get_root())
            
        
        self.logToMaster("AlignerAssessmentTarget Finished")

class ShredTarget(jobTree.scriptTree.target.Target):
    """
    A target that shreds a FASTA into reads
    
    """
    
    def __init__(self, contig_fasta, read_fasta, errors=False):
        """
        Make a target that shreds the contig_fasta into partly-overlapping reads
        with no Ns, producing the read_fasta.
        
        If errors is true, introduce some errors.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is super easy.
        super(ShredTarget, self).__init__(memory=2147483648)
        
        # Save everything
        self.contig_fasta = contig_fasta
        self.read_fasta = read_fasta
        self.erors = errors
        
        self.logToMaster("Creating ShredTarget")
        
   
    def run(self):
        """
        Map the reads.
        """
        
        self.logToMaster("Starting ShredTarget")
        
        # Make a temp directory for the index
        index_dir = sonLib.bioio.getTempDirectory(
            rootDir=self.getLocalTempDir())
        
        # Put together the arguments
        args = ["./shred.py", "--fastaIn", self.contig_fasta, "--fastaOut",
            self.read_fasta]
            
        if self.errors:
            args.append("--errors")
        
        # Run the shred script         
        check_call(self, args)
            
        self.logToMaster("ShredTarget Finished")

class MapReadsTarget(jobTree.scriptTree.target.Target):
    """
    A target that maps the reads in a FASTA to another FASTA
    
    """
    
    def __init__(self, reference, query, mappings_out, extra_args=[]):
        """
        Make a target that maps the query to the reference, using mapReads.
        Saves a TSV of mappings in mappings_out. Passes the optional extra_args
        to mapReads.
        
        """
        
        # Make the base Target. Ask for 8gb of memory since this is kinda hard,
        # and a bunch of CPUs so we can run in parallel.
        super(MapReadsTarget, self).__init__(memory=8589934592, cpu=32)
        
        # Save everything
        self.reference = reference
        self.query = query
        self.mappings_out = mappings_out
        self.extra_args = extra_args
        
        self.logToMaster("Creating MapReadsTarget")
        
   
    def run(self):
        """
        Map the reads.
        """
        
        self.logToMaster("Starting MapReadsTarget")
        
        # Make a temp directory for the index
        index_dir = sonLib.bioio.getTempDirectory(
            rootDir=self.getLocalTempDir())
        
        # Put together the arguments to invoke. Make sure to specify we want
        # lots of threads.
        args = (["../createIndex/mapReads", index_dir, self.reference, 
            self.query, "--alignment", self.mappings_out, "--threads", "32"] +
            self.extra_args)
            
        # Invoke the mapper
        check_call(self, args)
            
        self.logToMaster("MapReadsTarget Finished")
        
class BWATarget(jobTree.scriptTree.target.Target):
    """
    A target that maps the reads in a FASTA to another FASTA. The reference must
    be already indexed, because we don't want to fight over who gets to index
    it.
    
    """
    
    def __init__(self, reference, query, mappings_out):
        """
        Make a target that maps the query to the reference, using BWA. Saves a
        TSV of mappings in mappings_out.
        
        """
        
        # Make the base Target. Ask for 8gb of memory since this is kinda hard.
        super(BWATarget, self).__init__(memory=8589934592)
        
        # Save everything
        self.reference = reference
        self.query = query
        self.mappings_out = mappings_out
        
        self.logToMaster("Creating BWATarget")
        
   
    def run(self):
        """
        Map the reads.
        """
        
        self.logToMaster("Starting BWATarget")
        
        # TODO: we're assuming the reference is indexed. Check and complain to
        # the user if it isn't.
        
        # Get a SAM file to save as, temporarily.
        sam_file = sonLib.bioio.getTempFile(suffix=".sam",
            rootDir=self.getLocalTempDir())
        
        # Prepare the subprocess arguments
        args = ["bwa", "mem", self.reference, self.query]
        
        # Say we're going to do it
        self.logToMaster("Invoking {}".format(" ".join(args)))
        print("Invoking {}".format(" ".join(args)))
        
        # Start BWA
        process = subprocess.Popen(args, stdout=open(sam_file, "w"))
            
        if(process.wait() != 0):
            # If the BWA process died, complain.
            raise Exception("bwa failed with code {}".format(
                process.returncode))
                
        # Now we need to make the TSV output our parent expects.
        check_call(self, ["./sam2tsv.py", "--samIn", sam_file, "--tsvOut", 
            self.mappings_out, "--reference", self.reference])
            
        self.logToMaster("BWATarget Finished")
        
class TsvGeneCheckerTarget(jobTree.scriptTree.target.Target):
    """
    A target that checks a TSV of mappings against a set of BED files, and
    reports the number of mappings and the genes in various categories of
    correctness.
    
    """
    
    def __init__(self, mapping_file, gene_bed_files, class_count_file, 
        gene_set_file, class_beds={}):
        """
        Reads the mapping_file and the genes from the gene_bed_files list.
        Classifies every mapping in the alignment. Saves class counts to
        class_count_file ans a <class>\t<count> TSV. Saves genes with any
        mappings of each class in them to gene_set_file as a <class>\t<gene
        name> TSV.
        
        For each classification to BED filename mapping in class_beds, saves a
        BED file of mapped locations that are placed in that class.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(TsvGeneCheckerTarget, self).__init__(memory=2147483648)
        
        # Save the arguments
        self.mapping_file = mapping_file
        self.gene_bed_files = gene_bed_files
        self.class_count_file = class_count_file
        self.gene_set_file = gene_set_file
        self.class_beds = class_beds
        
        self.logToMaster("Creating TsvGeneCheckerTarget")
        
        
    def run(self):
        """
        Run the checks.
        """
        
        self.logToMaster("Starting TsvGeneCheckerTarget")
        
        # Prepare arguments
        args = (["./checkGenes.py", "--tsv", self.mapping_file, "--beds"] + 
            self.gene_bed_files + ["--classCounts", self.class_count_file,
            "--geneSets", self.gene_set_file])
            
        for classification, bed_file in self.class_beds.iteritems():
            # We need to pass each mapping in this dict along. Since the options
            # are predictably named, we generate them.
            args.append("--" + classification + "Bed")
            args.append(bed_file)
        
        # Make the call
        check_call(self, args)
        
        self.logToMaster("TsvGeneCheckerTarget Finished")


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
        help="directory to fill with gene sets and category counts per scheme")
    parser.add_argument("fastas", nargs="+",
        help="FASTA files to use, reference first")
    parser.add_argument("--geneBeds", nargs="+",
        help="BED files of genes to check against")
        
    
    
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
        from compareReads import AlignerAssessmentTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(AlignerAssessmentTarget(
        options.fastas, options.geneBeds, options.outDir))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
