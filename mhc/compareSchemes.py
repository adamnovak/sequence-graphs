#!/usr/bin/env python2.7
"""
compareSchemes.py: align pairs of MHC haplotypes with different LR merging
schemes (exact, inexact, credit, no credit), and compare the resulting
alignments to the official GRC alt loci alignments.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
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
    
    def __init__(self, fasta_list, true_maf, markov_model, gene_bed_dir, seed,
        stats_dir, hub_root):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, and comparing against the given truth MAF, using the given
        Markov model to measure coding costs, using the specified RNG seed, and
        writing statistics to one directory, and assembly hubs for different
        pairs of genomes and schemes to another.
        
        The BED files in gene_bed_dir will be used to evaluate alignment quality
        in light of gene annotations.
        
        The coverage statistics are just alignment coverage for each pair of
        genomes, in <genome number>\t<coverage fraction> TSVs per scheme.
        
        The alignment statistics are just <precision>\t<recall> TSVs per scheme.
        
        true_maf and markov_model may be None.
        
        Runs each subsequent FASTA against the first in each scheme.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(SchemeAssessmentTarget, self).__init__(memory=2147483648)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the filename of a MAF to compare all our MAFs against (or None)
        self.true_maf = true_maf
        
        # Save the Markov model file
        self.markov_model = markov_model
        
        self.gene_bed_dir = gene_bed_dir
        
        # Save the random seed
        self.seed = seed
        
        # Save the stats directory to populate
        self.stats_dir = stats_dir
        
        # And the place to put the assembly hubs
        self.hub_root = hub_root
        
        self.logToMaster(
            "Creating SchemeAssessmentTarget with seed {}".format(seed))
        
   
    def getSchemePlan(self):
        """
        Return a list of tuples describing schemes.
        
        """
        
        # Plan out all the schemes as mismatch, credit, min_context,
        # add_context, mult_context, min_coding_cost
        return [
            # Exact no credit min 100
            (False, False, 100, 0, 0, 0),
            # MultContext with and without min
            #(True, True, 0, 0, 4.0, 0),
            #(True, True, 60, 0, 4.0, 0),
            (True, True, 0, 0, 8.0, 0),
            #(True, True, 120, 0, 8.0, 0),
            # MultContext sans credit
            #(True, False, 0, 0, 8.0, 0),
        ]

    def getMonotonicSchemePlan(self):
        """
        Return a list of different tuples describing schemes.
        
        We want to look at addContext schemes without credit to try and catch
        nonmonotonicity.
        
        """
        
        # Plan out all the schemes as mismatch, credit, min_context,
        # add_context, mult_context, min_coding_cost
        return [
            # No credit
            (True, False, 0, 25, 0, 0),
            (True, False, 0, 50, 0, 0),
            (True, False, 0, 75, 0, 0),
            (True, False, 0, 100, 0, 0),
            # Credit
            (True, True, 0, 25, 0, 0),
            (True, True, 0, 50, 0, 0),
            (True, True, 0, 75, 0, 0),
            (True, True, 0, 100, 0, 0)
        ]

        
   
    def generateSchemes(self):
        """
        Whatever schemes I want to test at the moment
        
        """
        
        # Get a big list of tuples succinctly describing all the schemes we want
        # to run.
        scheme_plan = self.getSchemePlan()
        
        for mismatch, credit, min_context, add_context, mult_context, \
            min_coding_cost in scheme_plan:
            # Unpack each planned scheme
            
            # Start out with the context args
            extra_args = ["--context", str(min_context), "--addContext",
                str(add_context), "--multContext", str(mult_context)]
            
            # And give it a name to stick on our output files
            scheme_name = "E"
            if mismatch:
                # Add the args and scheme name component for mismatch
                extra_args.append("--mismatch")
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
                
            if min_context > 0:
                # Include min conext in the name if in use
                scheme_name += "Min{}".format(min_context)
            
            if add_context > 0:
                # Include additional context if in use. No need for badding
                # since we can now natural sort.
                scheme_name += "Add{}".format(add_context)
                
            if mult_context > 0:
                # Include multiplicative context if in use. Don't let any .s
                # into the scheme name since it needs to be a valid file
                # extension.
                scheme_name += "Mult{}".format(mult_context).replace(".", "p")
                
            if min_coding_cost > 0:
                # Say we're going to use a min coding cost
                extra_args.append("--minCodingCost")
                extra_args.append(str(min_coding_cost))
                    
                # Make sure we have a Markov model
                if self.markov_model is None:
                    raise Exception("We need a Markov model for coding cost!")
                    
                # Send it along to the merger
                extra_args.append("--markovModel")
                extra_args.append(self.markov_model)
            
                # Also, replace the decimal in the scheme name
                scheme_name += "Bits{}".format(min_coding_cost).replace(".",
                    "p")
                
            # Yield the name with the args.
            yield (scheme_name, extra_args)
        
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting SchemeAssessmentTarget")
        
        if not os.path.exists(self.stats_dir):
            # Make our out directory exist
            os.makedirs(self.stats_dir)
            
        if not os.path.exists(self.hub_root):
            # And our hubs directory
            os.makedirs(self.hub_root)
        
        # Look in the gene bed dir and grab all the bed files for alignment vs.
        # gene annotations comparison.
        if self.gene_bed_dir is not None:
            # This will hold all the gene BED files found
            gene_beds = list(glob.glob(self.gene_bed_dir + "/*.bed"))
        else:
            # No gene BED files will be used.
            gene_beds = []
        
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
        
        # This holds pairs of c2h and FASTA files
        c2h_fasta_pairs = []
        
        # This holds left and right context wiggle pairs. TODO: rename this to
        # max_context_pairs, and all the plain context stuff to max context
        # stuff.
        context_pairs = []
        
        # This holds left and right min context wiggle pairs
        min_context_pairs = []
        
        # This holds the genomes to which those pairs correspond
        pair_genomes = []
        
        # And this holds the schemes they are for
        pair_schemes = []
        
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
            coverage_filenames = [sonLib.bioio.getTempFile(suffix=".coverage",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's MAF alignment
            maf_filenames = [sonLib.bioio.getTempFile(suffix=".maf",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
            
            # And another one to hold each child's c2h alignment
            c2h_filenames = [sonLib.bioio.getTempFile(suffix=".c2h",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's fasta output for the c2h
            fasta_filenames = [sonLib.bioio.getTempFile(suffix=".fa",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # Save these so we can compare them to other schemes later.
            alignments_by_scheme[scheme] = maf_filenames
                
            # And another one to hold each child's adjacency component size
            # spectrum
            spectrum_filenames = [sonLib.bioio.getTempFile(suffix=".spectrum",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's indel lengths
            indel_filenames = [sonLib.bioio.getTempFile(suffix=".indels",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold the tandem duplication counts.
            # TODO: these are probably just wrong at the moment.
            tandem_filenames = [sonLib.bioio.getTempFile(suffix=".tandem",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one for left context length
            left_context_filenames = [sonLib.bioio.getTempFile(suffix=".wig",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one for right context length
            right_context_filenames = [sonLib.bioio.getTempFile(suffix=".wig",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And then for min contexts
            left_min_context_filenames = [sonLib.bioio.getTempFile(
                suffix=".wig", rootDir=self.getGlobalTempDir())
                for i in xrange(num_children)]
            right_min_context_filenames = [sonLib.bioio.getTempFile(
                suffix=".wig", rootDir=self.getGlobalTempDir())
                for i in xrange(num_children)]
            
            # Save the RNG state before clobbering it with the seed.
            random_state = random.getstate()
            
            # Seed the RNG after sonLib does whatever it wants with temp file
            # names. Each scheme will get runs with the same seeds, but that is
            # probably OK.
            random.seed(self.seed)
            
            for i, other_fasta in enumerate(self.fasta_list[1:]):
                # For each other FASTA to compare against
                
                # Pull out the files that this child should output.
                coverage_filename = coverage_filenames[i]
                maf_filename = maf_filenames[i]
                spectrum_filename = spectrum_filenames[i]
                indel_filename = indel_filenames[i]
                tandem_filename = tandem_filenames[i]
                c2h_filename = c2h_filenames[i]
                fasta_filename = fasta_filenames[i]
                left_context_filename = left_context_filenames[i]
                right_context_filename = right_context_filenames[i]
                left_min_context_filename = left_min_context_filenames[i]
                right_min_context_filename = right_min_context_filenames[i]
                
                # Make a child to produce those, giving it a seed. Make sure to
                # give it only two FASTAs, reference first, so that when it
                # shuffles the non-reference ones it doesn't do anything. Also
                # make sure to tell it to use the spectrum output file.
                self.addChildTarget(ReferenceStructureTarget(
                    [reference_fasta, other_fasta], random.getrandbits(256),
                    coverage_filename, maf_filename,
                    spectrum_filename=spectrum_filename,
                    indel_filename=indel_filename,
                    tandem_filename=tandem_filename, c2h_filename=c2h_filename,
                    fasta_filename=fasta_filename,
                    left_context_filename=left_context_filename,
                    right_context_filename=right_context_filename,
                    left_min_context_filename=left_min_context_filename,
                    right_min_context_filename=right_min_context_filename,
                    extra_args=extra_args))
                    
            
            # We need to check every MAF against the gene annotations.
            
            # What file will we use for each?
            checkgenes_filenames = [sonLib.bioio.getTempFile(suffix=".tsv",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # What target will we use for each? We can just give the whole list
            # of BEDs to each.
            checkgenes_targets = [MafGeneCheckerTarget(maf, gene_beds, out) 
                for maf, out in itertools.izip(maf_filenames, 
                checkgenes_filenames)]
            
            # Check each MAF against all the genes, and then concatenate the
            # answers for this scheme.
            followOns.append(SequenceTarget([RunTarget(checkgenes_targets), 
                ConcatenateTarget(checkgenes_filenames, 
                self.stats_dir + "/checkgenes." + scheme)]))
                
            # Make a follow-on job to merge all the child coverage outputs and
            # produce our coverage output file for this scheme.
            followOns.append(ConcatenateTarget(coverage_filenames, 
                self.stats_dir + "/coverage." + scheme))
                
            # Make a follow-on job to merge all the child spectrum outputs and
            # produce our spectrum output file for this scheme.
            followOns.append(ConcatenateTarget(spectrum_filenames, 
                self.stats_dir + "/spectrum." + scheme))
                
            # Make a follow-on job to merge all the child indel length outputs.
            followOns.append(ConcatenateTarget(indel_filenames, 
                self.stats_dir + "/indels." + scheme))
                
            # Make a follow-on job to merge all the child tandem duplication
            # count outputs.
            followOns.append(ConcatenateTarget(tandem_filenames, 
                self.stats_dir + "/tandem." + scheme))
                
            # Make a follow-on job to merge all the files with counts of how
            # many times a mappings are good/bad with respect to gene
            # annotations.
            followOns.append(ConcatenateTarget(checkgenes_filenames, 
                self.stats_dir + "/checkgenes." + scheme))
            
            if self.true_maf is not None:
                # We also need another target for comparing all these MAFs
                # against the truth, which we have been given.
                followOns.append(AlignmentTruthComparisonTarget(self.true_maf, 
                maf_filenames, random.getrandbits(256), 
                self.stats_dir + "/truth." + scheme))
                
            # Save the c2h and FASTA files we made along with the scheme and
            # genomes we used.
            c2h_fasta_pairs += zip(c2h_filenames, fasta_filenames)
            # And the genomes they were for.
            pair_genomes += [os.path.splitext(fasta)[0] 
                for fasta in self.fasta_list[1:]]
            # And record that each came from this scheme. Needs to keep the same
            # lenght as pair_genomes.
            pair_schemes += [scheme for _ in self.fasta_list[1:]]
            
            # Also save the context wiggle pairs.
            context_pairs += zip(left_context_filenames,
                right_context_filenames)
                
            # And the min context wiggle pairs.
            min_context_pairs += zip(left_min_context_filenames,
                right_min_context_filenames)
            
        
        # What genome pairs did we run, in order? Make sure to strip extensions.
        genome_pairs = [(os.path.splitext(self.fasta_list[0])[0], 
            os.path.splitext(other)[0]) 
            for other in self.fasta_list[1:]]
                
        # We want to combine all our c2h files into one massive hub.
        
        # What HAL will we use?
        merged_hal = sonLib.bioio.getTempFile(suffix=".hal",
                rootDir=self.getGlobalTempDir())
        os.unlink(merged_hal)
        
        # What merged c2h and fasta files will we use?
        merged_c2h = sonLib.bioio.getTempFile(suffix=".c2h",
                rootDir=self.getGlobalTempDir())
        merged_fasta = sonLib.bioio.getTempFile(suffix=".fa",
                rootDir=self.getGlobalTempDir())
                
        # Where should we keep our left and right max context wiggles?
        left_context_dir = sonLib.bioio.getTempDirectory(
                rootDir=self.getGlobalTempDir()) + "/leftMaxContext"
        right_context_dir = sonLib.bioio.getTempDirectory(
                rootDir=self.getGlobalTempDir()) + "/rightMaxContext"
                
        # Where should we keep our left and right max context wiggles?
        left_min_context_dir = sonLib.bioio.getTempDirectory(
                rootDir=self.getGlobalTempDir()) + "/leftMinContext"
        right_min_context_dir = sonLib.bioio.getTempDirectory(
                rootDir=self.getGlobalTempDir()) + "/rightMinContext"
        
        # What suffixes should we put on genomes?
        suffixes = [scheme for scheme in pair_schemes]
        
        # And what are the genomes with suffixes?
        genomes_with_suffixes = [genome + suffix for genome, suffix in 
            itertools.izip(pair_genomes, suffixes)]
        
        # What genome is the reference?
        reference = os.path.splitext(self.fasta_list[0])[0]
        
        # Make a target to make the merged c2h/fasta files
        merged_c2h_target = C2hMergeTarget(c2h_fasta_pairs, suffixes, 
            merged_c2h, merged_fasta)
            
        # And a target to make the hal from them
        merged_hal_target = HalTarget(merged_c2h, merged_fasta, reference, 
            genomes_with_suffixes, merged_hal)
            
        # And a target to sort out all of the left wiggles
        left_wiggle_target = WiggleCollateTarget(
            [pair[0] for pair in context_pairs], pair_genomes, suffixes,
            left_context_dir)
        # And the right ones
        right_wiggle_target = WiggleCollateTarget(
            [pair[1] for pair in context_pairs], pair_genomes, suffixes,
            right_context_dir)
            
        # And for the minimum contexts for uniqueness
        left_min_wiggle_target = WiggleCollateTarget(
            [pair[0] for pair in min_context_pairs], pair_genomes, suffixes,
            left_min_context_dir)
        right_min_wiggle_target = WiggleCollateTarget(
            [pair[1] for pair in min_context_pairs], pair_genomes, suffixes,
            right_min_context_dir)
            
        # And a target to make a hub from that
        merged_hub_target = AssemblyHubOnlyTarget(merged_hal, self.hub_root + 
            "/all", wiggle_dirs=[left_context_dir, right_context_dir, 
            left_min_context_dir, right_min_context_dir])
            
        # Do those last ones in order. TODO: express actual order restrictions
        # with a few more intermediate targets.
        followOns.append(SequenceTarget([left_wiggle_target, 
            right_wiggle_target, left_min_wiggle_target, 
            right_min_wiggle_target, merged_c2h_target, merged_hal_target,
            merged_hub_target]))
            
        # So we need to run all of the follow-ons in parallel after this job
        self.setFollowOnTarget(RunTarget(followOns))
                
        self.logToMaster("SchemeAssessmentTarget Finished")

class MafGeneCheckerTarget(jobTree.scriptTree.target.Target):
    """
    A target that checks a MAF alignment against a set of BED files, and reports
    the number of mappings in various categories of correctness.
    
    """
    
    def __init__(self, maf_file, bed_files, output_file):
        """
        Given a MAF filename and a list of BED filenames, check the MAF against
        the genes in the BEDs and write the output statistics to the given
        output file.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(MafGeneCheckerTarget, self).__init__(memory=2147483648)
        
        # Save the arguments
        self.maf_file = maf_file
        self.bed_files = bed_files
        self.output_file = output_file
        
        self.logToMaster("Creating MafGeneCheckerTarget")
        
        
    def run(self):
        """
        Run the checks.
        """
        
        self.logToMaster("Starting MafGeneCheckerTarget")
        
        # Prepare arguments
        args = (["./checkGenes.py", "--maf", self.maf_file, "--beds"] + 
            self.bed_files + ["--out", self.output_file])
        
        # Make the call
        check_call(self, args)
        
        self.logToMaster("MafGeneCheckerTarget Finished")

class C2hMergeTarget(jobTree.scriptTree.target.Target):
    """
    A target that merges a list of c2h/FASTA pairs on a shared reference, where
    each is a unary tree of the given reference and a genome from the given list
    of genomes. The corresponding suffix from the list of suffixes is assigned
    to each genome. Saves a merged c2h/fasta pair.
    
    """
    
    def __init__(self, c2h_fasta_pairs, suffixes, merged_c2h, merged_fasta):
        """
        Make a merged c2h/fasta pair from the given pairs. The reference is the
        genome with the given name, and the other genome in each pair is in the
        corresponding entry in genome_names. Suffixes is a list of suffixes to
        add to the leaves from each file when merging it in.
        
        """
        
        # Make the base Target. Ask for 8gb of memory since the intermediate
        # graphs can get big.
        super(C2hMergeTarget, self).__init__(memory=8589934592)
        
        # Save the files to merge
        self.c2h_fasta_pairs = c2h_fasta_pairs
        
        # And the suffixes to apply to all the genomes
        self.suffixes = suffixes
        
        # And the merged files to write
        self.merged_c2h = merged_c2h
        self.merged_fasta = merged_fasta
        
        self.logToMaster("Creating C2hMergeTarget")
        
        
    def run(self):
        """
        Do all the merges.
        """
        
        self.logToMaster("Starting C2hMergeTarget")
        
        if len(self.c2h_fasta_pairs) == 0:
            raise Exception("No alignments at all!")
            
        # Start preparing arguments
        args = ["../createIndex/cactusMerge", self.merged_c2h, self.merged_fasta]
        
        for (c2h, fasta), suffix in itertools.izip(self.c2h_fasta_pairs,
            self.suffixes):
            
            # Make args for each thing you want to merge in
            args.append("--c2h")
            args.append(c2h)
            args.append("--fasta")
            args.append(fasta)
            args.append("--suffix")
            args.append(suffix)
        
        # Make the call
        check_call(self, args)
        
        self.logToMaster("C2hMergeTarget Finished")

class C2hChangeRootTarget(jobTree.scriptTree.target.Target):
    """
    A target that changes the root sequence of star tree a c2h/fasta pair.
    
    """
    
    def __init__(self, c2h_file, fasta_file, new_event, new_sequence, 
        new_root_fasta, old_root_offset, c2h_out, fasta_out):
        """
        Change the root sequence in the c2h/fasta pair to have the new event and
        sequence names, and add the new root fasta to the paired fasta. Also
        adapts block coordinates on the root to accomplish a lift-over.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(C2hChangeRootTarget, self).__init__(memory=2147483648)
        
        # Save everything
        # C2H file to adjust
        self.c2h_file = c2h_file
        # And its associated fasta
        self.fasta_file = fasta_file
        # The new event name to use for the root
        self.new_event = new_event
        # The new sequence name to use for the root (we assume one root
        # sequence, and that it is first)
        self.new_sequence = new_sequence
        # The FASTA of the sequence for the new root, with ID event.sequence
        self.new_root_fasta = new_root_fasta
        # The offset of the old root's first base in the new root, in 0-based
        # coordinates
        self.old_root_offset = old_root_offset
        # The c2h file to write
        self.c2h_out = c2h_out
        # The fasta file to go with it
        self.fasta_out = fasta_out
        
        self.logToMaster("Creating C2hChangeRootTarget")
        
        
    def run(self):
        """
        Do the liftover and rewriting.
        """
        
        self.logToMaster("Starting C2hChangeRootTarget")
        
        # Put the old FASTA as the output
        
        # Copy the new one, counting bases
        
        # Open up the old c2h file
        
        # The first line gets ammended with our new event and sequnce names in
        # columns 1 and 2
        
        # Then we do a bottom block up to the offset we are supposed to add
        
        # Then for each subsequent bottom block we up the start position
        
        # Then at the end position of the last one (next line is another s line)
        # we put a block out to the end of the FASTA sequence we read before.
        
        # Then we copy all remaining lines as normal.
                
        self.logToMaster("C2hChangeRootTarget Finished")
        
class HalTarget(jobTree.scriptTree.target.Target):
    """
    A target that makes a c2h/fasta pair with a root and a list of leaves into a
    hal file.
    
    """
    
    def __init__(self, c2h_file, fasta_file, reference_name, genome_names, hal):
        """
        Make a hal file from a star tree c2h/fasta pair, given the name of the
        root and the names of all the leaves.
        
        """
        
        # Make the base Target. Ask for 16gb of memory since this is hard.
        super(HalTarget, self).__init__(memory=17179869184)
        
        # Save the files to merge
        self.c2h_file = c2h_file
        self.fasta_file = fasta_file
        
        # And the name of the reference to root the hal with.
        self.reference = reference_name
        
        # And the genome names that are used
        self.genomes = genome_names
        
        # And the hal to write
        self.hal = hal
        
        self.logToMaster("Creating HalTarget")
        
        
    def run(self):
        """
        Make a hal from a c2h and fasta.
        """
        
        self.logToMaster("Starting HalTarget")
        
        # Compose the tree.      
        tree = "(" + ",".join(self.genomes) + ")" + self.reference + ";"
                
        # Make the HAL. Do it in memory to be faster.
        check_call(self, ["halAppendCactusSubtree", "--inMemory", self.c2h_file,
            self.fasta_file, tree, self.hal])
            
                
        self.logToMaster("HalTarget Finished")
        
        
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
                # (in the above directory) for this genome. Make sure to
                # preserve feature names with
                # <https://www.biostars.org/p/109041/#109046>
                handle = subprocess.Popen(["bedtools", "merge", "-i", 
                    intermediate, "-c", "4", "-o", "distinct"], 
                    stdout=open(out_filename, "w"))
                if handle.wait() != 0:
                    raise RuntimeError("Could not merge " + intermediate) 
                    
            # Otherwise, don't bother making the actual output file.
        
        
        self.logToMaster("BedSplitTarget Finished")
        
class WiggleCollateTarget(jobTree.scriptTree.target.Target):
    """
    A target which collates wiggle files by suffixed genome and applies
    suffixes.
    
    Produces an output directory of wiggles of the form 
    <out_dir>/<genome><suffix>/<genome><suffix>.wig
    
    """

    def __init__(self, wiggle_files, genomes, suffixes, out_dir):
        """
        Given a list of wiggle files, a list of the genome name for each wiggle
        file, and a list of suffixes to be applied to the genome names, apply
        the suffixes and organize the wiggle files in out_dir for use by
        hal2assemblyHub.py.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(WiggleCollateTarget, self).__init__(memory=2147483648)
        
        # Save the arguments
        self.wiggle_files = wiggle_files
        self.genomes = genomes
        self.suffixes = suffixes
        self.out_dir = out_dir
        
            
    def run(self):
        """
        Run this target and do the splitting.
        
        """
    
        self.logToMaster("Starting WiggleCollateTarget")
        
        # Work out the suffixed genome nemae
        suffixed_names = [genome + suffix 
            for (genome, suffix) in itertools.izip(self.genomes, self.suffixes)]
            
        for suffixed in suffixed_names:
            # Make a directory for each genome with its suffix
            os.makedirs(self.out_dir + "/" + suffixed)
        
        # Work out what each genome's wiggle should be called.
        out_filenames = [self.out_dir + "/{}/{}.wig".format(suffixed, suffixed) 
            for suffixed in suffixed_names]
            
        for in_filename, original_genome, new_genome, out_filename in \
            itertools.izip(self.wiggle_files, self.genomes, suffixed_names, 
            out_filenames):
            
            # For every wiggle we have to copy (and update)
            
            # Open the file to write
            out_file = open(out_filename, "w")
            
            for line in open(in_filename):
                # For each line we want to copy, change the genome name
                line = line.replace(original_genome, new_genome)
                # Then write it
                out_file.write(line)
            
        
        self.logToMaster("WiggleCollateTarget Finished")
        
        
class AssemblyHubOnlyTarget(jobTree.scriptTree.target.Target):
    """
    A target that makes a single assembly hub showing just a HAL.
    
    """
    
    def __init__(self, hal, hub, wiggle_dirs=None):
        """
        Takes a hal alignment and makes a hub.
        
        """
        
        # Make the base Target. Ask for 8gb of memory and several CPUs since
        # this is hard.
        super(AssemblyHubOnlyTarget, self).__init__(memory=8589934592, cpu=32)
        
        # Save the alignment file
        self.hal = hal
        
        # And the hub directory
        self.hub = hub
        
        # And the wiggle directory list
        self.wiggle_dirs = wiggle_dirs
        
        self.logToMaster("Creating AssemblyHubOnlyTarget")
        
        
    def run(self):
        """
        Make the assembly hubs.
        """
        
        self.logToMaster("Starting AssemblyHubOnlyTarget")
        
                
        # Make a temporary directory for the jobTree tree. Make sure it's
        # absolute.
        tree_dir = os.path.abspath(sonLib.bioio.getTempDirectory(
            rootDir=self.getLocalTempDir()))
            
        # Make sure it doesn't exist yet
        os.rmdir(tree_dir)
        
        # We need to do the hal2assemblyHub.py call in its own directory since
        # it makes temp files in the current directory.
        working_directory = sonLib.bioio.getTempDirectory(
            rootDir=self.getLocalTempDir())
            
        # Where should we come back to when done?
        original_directory = os.getcwd()
            
        # Turn our arguments into absolute paths.
        # May or may not be necessary.
        hal_abspath = os.path.abspath(self.hal)
        hub_abspath = os.path.abspath(self.hub)
        
        if self.wiggle_dirs is not None:
            # Make wiggle directories absolute paths if needed.
            self.wiggle_dirs = [os.path.abspath(path) 
                for path in self.wiggle_dirs]
        
        # Assemble the command
        command = ["hal2assemblyHub.py", 
            hal_abspath, hub_abspath, "--jobTree", 
            tree_dir, "--maxThreads", "32", # No LOD
            "--cpHalFileToOut", "--noUcscNames"]
            
        if self.wiggle_dirs is not None:
            # Add the option to use these wiggle dirs. TODO: these directories
            # can't have commas in their names
            command.append("--wigDirs")
            command.append(",".join(self.wiggle_dirs))
        
        # Go in the temp directory and run the script
        os.chdir(working_directory)
        check_call(self, command)
        
        # Go back to the original directory
        os.chdir(original_directory)
            
        
        self.logToMaster("AssemblyHubOnlyTarget Finished")

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
    parser.add_argument("--trueMaf", default=None,
        help="filename of a MAF to compare all generated MAFs against")
    parser.add_argument("--markovModel", default=None,
        help="filename of a Markov model to model query sequences with")
    parser.add_argument("--geneBedDir", default=None,
        help="directory in which to find gene annotations to judge alignments")
        
    
    
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
        from compareSchemes import SchemeAssessmentTarget, BedSplitTarget, \
            MafGeneCheckerTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(SchemeAssessmentTarget(
        options.fastas, options.trueMaf, options.markovModel, 
        options.geneBedDir, options.seed, options.outDir, 
        options.outDir + "/hubs"))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
