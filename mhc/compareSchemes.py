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
        agreement_basename, spectrum_basename, truth_basename, hub_root):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, and comparing against the given truth MAF, using the
        specified RNG seed, and writing coverage statistics to files with
        specified coverage base name, a set of adjacency component size spectra
        to files with the specified spectrum basename, stats for agreement
        between the schemes to files anmed after the agreement basename, a
        comparison against the true MAF to files named after the truth_basename,
        and a directory full of assembly hubs for different pairs of genomes and
        schemes.
        
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
        
        # And the place to put the assembly hubs
        self.hub_root = hub_root
        
        self.logToMaster(
            "Creating SchemeAssessmentTarget with seed {}".format(seed))
        
   
    def generateSchemes(self):
        """
        Yield tuples of (scheme name, extra args list) for all the schemes we
        want to do.
        
        """
        
        for mismatch, credit in itertools.product([False, True], repeat=2):
            # For all combinations of mismatch and credit
            for min_context in [50, 100]:
                # And min context length
                for add_context in [0, 30]:
                    # And additional context
                    
                    # Start out with the context args
                    extra_args = ["--context", str(min_context), "--addContext",
                        str(add_context)]
            
                    # And give it a name to stick on our output files
                    scheme_base = "Exact"
                    if mismatch:
                        # Add the args and scheme name component for mismatch
                        extra_args.append("--mismatch")
                        extra_args.append("--mismatches")
                        extra_args.append("1")
                        scheme_base = "Inexact"
                    if credit:
                        # Add the args and scheme name component for credit
                        extra_args.append("--credit")
                        scheme_base += "Credit"
                        
                    # Put together the full scheme name and yield it with the
                    # args.
                    yield ("{}Min{}Add{}".format(scheme_base, min_context, 
                        add_context), extra_args)
        
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting SchemeAssessmentTarget")
        
        # We need a few different follow-on jobs.
        followOns = []    
        
        # We need a directory to save out tree of BED files under for making the
        # assembly hubs.
        bed_root = sonLib.bioio.getTempDirectory(rootDir=self.getGlobalTempDir())
        
        # We'll keep track of our random state so we can seed without messing up
        # temp filenames.
        random_state = None
        
        # This will hold lists of MAF alignments, in FASTA pair order, by scheme
        # name.
        alignments_by_scheme = {}
        
        # And a similar structure for HALs
        hals_by_scheme = {}
        
        for scheme, extra_args in self.generateSchemes():
            # Work out all the schemes we want to run.
            
            self.logToMaster("Preparing for scheme {}...".format(scheme))

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
                
            # And another one to hold each child's HAL alignment
            hal_filenames = [sonLib.bioio.getTempFile(suffix="hal",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            for hal in hal_filenames:
                # Make sure the HALs don't exist yet.
                os.unlink(hal)
                
            # Save these so we can compare them to other schemes later.
            alignments_by_scheme[scheme] = maf_filenames
            hals_by_scheme[scheme] = hal_filenames
                
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
                hal_filename = hal_filenames[i]
                spectrum_filename = spectrum_filenames[i]
                
                # Make a child to produce those, giving it a seed. Make sure to
                # give it only two FASTAs, reference first, so that when it
                # shuffles the non-reference ones it doesn't do anything. Also
                # make sure to tell it to use the spectrum output file.
                self.addChildTarget(ReferenceStructureTarget(
                    [reference_fasta, other_fasta], random.getrandbits(256), 
                    stats_filename, maf_filename, hal_filename=hal_filename,
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
        
        
        
        # What genome pairs did we run, in order? Make sure to strip extensions.
        genome_pairs = [(os.path.splitext(self.fasta_list[0])[0], 
            os.path.splitext(other)[0]) 
            for other in self.fasta_list[1:]]
                
        # Now we have all the followons for concatenating our stats and
        # comparing against the truth, we need a followon for comparing
        # corresponding MAFs from the same FASTA pair with different schemes,
        # for all combinations of schemes.
        agreement_target = AlignmentSchemeAgreementTarget(alignments_by_scheme,
            genome_pairs, random.getrandbits(256), self.agreement_basename,
            bed_root)
            
        # After that though, we need to take the BED files and the HAL files and
        # make assembly hubs in hub_root.
        hubs_target = AssemblyHubsTarget(hals_by_scheme, genome_pairs, bed_root,
            self.hub_root)
            
        # Do those two things in order.
        followOns.append(SequenceTarget([agreement_target, hubs_target]))
            
            
        # So we need to run all of the follow-ons in parallel after this job
        self.setFollowOnTarget(RunTarget(followOns))
                
        self.logToMaster("SchemeAssessmentTarget Finished")

class AlignmentSchemeAgreementTarget(jobTree.scriptTree.target.Target):
    """
    A target that, given lists of corresponding alignments in a dict by scheme,
    compares corresponding alignments across all pairs of schemes and collates
    the results.
    
    """
    
    def __init__(self, alignments_by_scheme, genome_names, seed, 
        agreement_basename, bed_root):
        """
        Takes a dict of lists of alignments by scheme (arranged so alignments at
        the same indices are of the same FASTAs), a list of pairs of genomes
        which are compared by each alignment. Sends output to a file basename to
        name all results after, and a directory in which to build a fairly
        oragnized tree of BED files.
        
        Compares each scheme to each other scheme on corresponding genome pairs.
        
        Uses the seed to generate mafComparator seeds and temporary file names,
        so don't run two copies with the same seed.
        
        Each (scheme1, scheme2, genome pair) combination produces a pair of BED
        files of positions where the two schemes differ. These are saved in:
            
        <bed_root>/<scheme1>-<scheme2>-<genome1>-<genome2>/<genome>/<genome>.bed
        
        where scheme 1 is the lexicographically smaller scheme,  genome 1 and
        genome 2 are the genomes as specified in the pair, and <genome> is the
        genome on which the BED coordinates are given.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AlignmentSchemeAgreementTarget, self).__init__(memory=2147483648)
        
        self.seed = seed
        
        # Save the alignment names
        self.alignments_by_scheme = alignments_by_scheme
        
        # And the genome name pairs
        self.genome_names = genome_names
        
        # And the agreement file name
        self.agreement_basename = agreement_basename
        
        # And the BED root
        self.bed_root = bed_root
        
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
            
            for comparison_filename, alignment1, alignment2, compared_genomes \
                in itertools.izip(comparison_filenames, 
                self.alignments_by_scheme[scheme1], 
                self.alignments_by_scheme[scheme2], self.genome_names):
                
                # For corresponding alignments
                
                # Find a place for a temporary BED file
                temp_bed = sonLib.bioio.getTempFile(
                    rootDir=self.getGlobalTempDir()) 
                
                # Work out where the split up BED files should land
                bed_dir = self.bed_root + "/{}-{}-{}-{}".format(scheme1, 
                    scheme2, compared_genomes[0], compared_genomes[1])
                
                # Make sure that directory exists.
                os.makedirs(bed_dir)
                
                # Compare the alignments
                self.addChildTarget(AlignmentComparisonTarget(alignment1, 
                    alignment2, random.getrandbits(256), comparison_filename,
                    bed=temp_bed))
                    
                # Split the BED file from this comparison out by genome.
                # Override the feature names to be meaningful and list the
                # scheme we're comparing against.
                followOns.append(BedSplitTarget(temp_bed, compared_genomes, 
                    bed_dir, feature_name="{}-{}".format(scheme1, scheme2)))
                    
            # Concatenate all the comparisons into a file named after both
            # schemes.
            followOns.append(ConcatenateTarget(comparison_filenames,
                self.agreement_basename + "." + scheme1 + "." + scheme2))
        
        
        # We need to run all of the follow-ons in parallel after this job
        self.setFollowOnTarget(RunTarget(followOns))
                
        self.logToMaster("AlignmentSchemeAgreementTarget Finished")
        
class BedSplitTarget(jobTree.scriptTree.target.Target):
    """
    A target which splits a BED up by genome. Also fixes them up a bit.
    
    """

    def __init__(self, bed_file, genomes, out_dir, feature_name=None):
        """
        Split the given bed file per genome into <genome>/<genome>.bed in the
        given output directory, using the given list of genomes.
        
        out_dir must exist.
        
        If a feature_name is specified, it is used to rename all the BED
        feratures.
        
        Does not leave a trailing newline.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(BedSplitTarget, self).__init__(memory=2147483648)
        
        # Save the arguments
        self.bed_file = bed_file
        self.genomes = genomes
        self.out_dir = out_dir
        self.feature_name = feature_name
        
            
    def run(self):
        """
        Run this target and do the splitting.
        
        """
    
        self.logToMaster("Starting BedSplitTarget")
        
        # Work out output file names
        out_filenames = [self.out_dir + "/{}/{}.bed".format(genome, genome) 
            for genome in self.genomes]
            
        # Work out temp files to write to first
        temp_filenames = [sonLib.bioio.getTempFile(
            rootDir=self.getLocalTempDir()) for _ in out_filenames]
        
        # Open the temp BED file for each genome
        temp_files = [open(temp_filename, "w")
            for temp_filename in temp_filenames]
            
        
        # We collapse adjacent identical things as we merge. So this dict holds
        # the (start, end, name) tuple last written for each genome.
        last_interval = {}
            
        # We need to keep track of whether anything has been written to each
        # file yet, so we know when to put newlines, so we don't leave trailing
        # newlines. We put a genome in this set when we have written to its
        # file.
        content_written = set()
            
        for line in open(self.bed_file):
            # For each BED record, grab all the parts.
            parts = line.strip().split()
            
            if self.feature_name is not None:
                if len(parts) > 3:
                    # Rename features with names
                    parts[3] = self.feature_name
                elif len(parts) == 3:
                    # Add names to features without them.
                    parts.append(self.feature_name)
            
            if len(parts) == 1 and parts[0] == "":
                # Skip any blank lines
                continue
            
            for genome, temp_file in itertools.izip(self.genomes, temp_files):
                if parts[0] == genome:
                    # Send the line to the file corresponding to the matching
                    # genome. TODO: use a dict or something instead of scanning
                    # for matching names.
                    
                    if genome in content_written:
                        # Terminate the previous line
                        temp_file.write("\n")
                    # Write the line
                    temp_file.write("\t".join(parts))
                    # Remember to end it if we need to write another.
                    content_written.add(genome)
                    
        for temp_file in temp_files:
            # Close up all our output files.
            temp_file.close()
            
        for genome, temp_filename, out_filename in itertools.izip(self.genomes, 
            temp_filenames, out_filenames):
            
            # Fix up non-empty BEDs with bedtools, and don't pass on empty ones.
            
            if genome in content_written:
                # We have content
                
                if not os.path.exists(self.out_dir + "/" + genome):
                    # Make the directory for the final BED file.
                    os.mkdir(self.out_dir + "/" + genome)
            
                # We need an intermediate file for sorting.
                intermediate = sonLib.bioio.getTempFile(
                    rootDir=self.getLocalTempDir())
                
                # Sort the BED file    
                handle = subprocess.Popen(["sort", "-k1,1", "-k2,2n"], 
                    stdin=open(temp_filename), stdout=open(intermediate, "w"))
                if handle.wait() != 0:
                    raise RuntimeError("Could not sort " + temp_filename)
                
                # Merge the BED file and write to the pre-calculate output file
                # (in the above directory) for this genome.
                handle = subprocess.Popen(["bedtools", "merge", "-i", 
                    intermediate], stdout=open(out_filename, "w"))
                if handle.wait() != 0:
                    raise RuntimeError("Could not merge " + intermediate) 
                    
            # Otherwise, don't bother making the actual output file.
        
        
        self.logToMaster("BedSplitTarget Finished")
        
class AssemblyHubsTarget(jobTree.scriptTree.target.Target):
    """
    A target that makes assembly hubs showing a bunch of HALs and how they
    compare to each other.
    
    """
    
    def __init__(self, hals_by_scheme, genome_names, bed_root, hub_root):
        """
        Takes a dict of lists of HAL alignments by scheme (arranged so
        alignments at the same indices are of the same FASTAs), and a list of
        the pairs of genomes used to make each hal (in the same order).
        
        Also takes the root of a directory tree full of BED files:
            
        <bed_root>/<scheme1>-<scheme2>-<genome1>-<genome2>/<genome>/<genome>.bed
        
        where scheme 1 is the lexicographically smaller scheme,  genome 1 and
        genome 2 are the genomes as specified in the pair, and <genome> is the
        genome on which the BED coordinates are given.
        
        Produces a set of assembly hubs in:
        
        <hub_root>/<scheme>/<genome1>-<genome2>
        
        with one for each genome pair under each scheme. hub_root will be
        created if it does not exist.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AssemblyHubsTarget, self).__init__(memory=2147483648)
        
        # Save the alignment names
        self.hals_by_scheme = hals_by_scheme
        
        # And the genome name pairs
        self.genome_names = genome_names
        
        # And the BED root
        self.bed_root = bed_root
        
        # And the hub root
        self.hub_root = hub_root
        
        self.logToMaster("Creating AssemblyHubsTarget")
        
        
    def run(self):
        """
        Make the assembly hubs.
        """
        
        self.logToMaster("Starting AssemblyHubsTarget")
        
        for scheme in self.hals_by_scheme.iterkeys():
            
            # Figure out where hubs of this scheme go
            scheme_root = self.hub_root + "/" + scheme
            if not os.path.exists(scheme_root):
                os.makedirs(scheme_root)
            
            for i, (genome1, genome2) in enumerate(self.genome_names):
                # Grab pairs of genomes, and their indices.
                
                # Work out its directory name
                genome_string = "{}-{}".format(genome1, genome2)
                
                # Determine the directory for the assembly hub.
                hub_dir = scheme_root + "/" + genome_string
                
                # Make a child target to make this actual hub
                self.addChildTarget(AssemblyHubTarget(
                    self.hals_by_scheme[scheme][i], scheme, 
                    self.hals_by_scheme.keys(), (genome1, genome2), 
                    self.bed_root, hub_dir))
                
        self.logToMaster("AssemblyHubsTarget Finished")
        
class AssemblyHubTarget(jobTree.scriptTree.target.Target):
    """
    A target that makes a single assembly hub showing a HAL and some BEDs.
    
    """
    
    def __init__(self, hal, scheme, scheme_names, genome_pair, bed_root, hub):
        """
        Takes a hal alignment and a scheme name, along with the names of all the
        other schemes to compare against (possibly including itself, which will
        be skipped).
        
        Also takes the pair (as a tuple) of genomes we are making the hub for.
        
        Also takes the root of a directory tree full of BED files:
            
        <bed_root>/<scheme1>-<scheme2>-<genome1>-<genome2>/<genome>/<genome>.bed
        
        where scheme 1 is the lexicographically smaller scheme,  genome 1 and
        genome 2 are the genomes as specified in the pair, and <genome> is the
        genome on which the BED coordinates are given.
        
        Produces an assembly hub in the hub directory.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AssemblyHubTarget, self).__init__(memory=2147483648)
        
        # Save the alignment file
        self.hal = hal
        
        # And the scheme
        self.scheme = scheme
        
        # And all the other schemes to compare against
        self.scheme_names = scheme_names
        
        # And the genome name pair
        self.genome_pair = genome_pair
        
        # And the BED root
        self.bed_root = bed_root
        
        # And the hub directory
        self.hub = hub
        
        self.logToMaster("Creating AssemblyHubsTarget")
        
        
    def run(self):
        """
        Make the assembly hubs.
        """
        
        self.logToMaster("Starting AssemblyHubTarget")
        
                
        # Make a temporary directory for the jobTree tree
        tree_dir = sonLib.bioio.getTempDirectory(
            rootDir=self.getLocalTempDir())
            
        # Make sure it doesn't exist yet
        os.rmdir(tree_dir)
        
        # Get all the pairs of schemes involving this one
        possible_bed_pairs = ([(self.scheme, other) 
            for other in self.scheme_names
            if other != self.scheme] + 
            [(other, self.scheme) 
            for other in self.scheme_names
            if other != self.scheme])
            
        # Get the directory each pair of schemes would produce for this pair of
        # genomes
        bed_dirs = [self.bed_root + "/{}-{}-{}".format(scheme1, scheme2,
            self.genome_pair) for (scheme1, scheme2) in possible_bed_pairs]
            
        # Keep the ones that exist and make a string of them.
        bed_dirs_string = ",".join([directory for directory in bed_dirs
            if os.path.exists(directory)])
        
        # We want to make an assembly hub like so: 
                       
        # hal2assemblyHub.py data/alignment1.hal data/alignment1.hub
        # --jobTree data/tree --bedDirs data/beddir --hub=nocredit
        # --shortLabel="No Credit" --lod --cpHalFileToOut --noUcscNames
        check_call(self, ["hal2assemblyHub.py", 
            self.hal, self.hub, "--jobTree", 
            tree_dir, "--bedDirs", bed_dirs_string, "--shortLabel", 
            self.scheme, "--lod", "--cpHalFileToOut", "--noUcscNames"])
        
        self.logToMaster("AssemblyHubTarget Finished")

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
    parser.add_argument("hubRoot", 
        help="directory to populate with assembly hubs")
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
            AlignmentSchemeAgreementTarget, BedSplitTarget, \
            AssemblyHubsTarget, AssemblyHubTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(SchemeAssessmentTarget(
        options.fastas, options.trueMaf, options.seed, options.coverageBasename, 
        options.agreementBasename, options.spectrumBasename,
        options.truthBasename, options.hubRoot))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
