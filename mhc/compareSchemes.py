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
    
    def __init__(self, fasta_list, true_maf, markov_model, seed, stats_dir,
        hub_root):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, and comparing against the given truth MAF, using the given
        Markov model to measure coding costs, using the specified RNG seed, and
        writing statistics to one directory, and assembly hubs for different
        pairs of genomes and schemes to another.
        
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
        
        # Save the random seed
        self.seed = seed
        
        # Save the stats directory to populate
        self.stats_dir = stats_dir
        
        # And the place to put the assembly hubs
        self.hub_root = hub_root
        
        self.logToMaster(
            "Creating SchemeAssessmentTarget with seed {}".format(seed))
        
   
    def generateSchemes(self):
        """
        Whatever schemes I want to test at the moment
        
        """
        
        # Plan out all the schemes as mismatch, credit, min_context,
        # add_context, mult_context, min_coding_cost
        scheme_plan = [
                        (False, False, 100, 0, 0, 0),
                        # Specified Min
                        (True, True, 60, 0, 0, 0),
                        (True, True, 80, 0, 0, 0),
                        (True, True, 100, 0, 0, 0),
                        (True, True, 120, 0, 0, 0),
                        # Specified add
                        (True, True, 0, 25, 0, 0),
                        (True, True, 0, 50, 0, 0),
                        (True, True, 0, 75, 0, 0),
                        (True, True, 0, 100, 0, 0),
                        # Coding cost context
                        (True, True, 0, 0, 0, 20),
                        (True, True, 0, 0, 0, 80),
                        (True, True, 0, 0, 0, 140),
                        (True, True, 0, 0, 0, 200),
                        (True, True, 0, 0, 0, 260),
                        # MultContext
                        (True, True, 0, 0, 1.25, 0),
                        (True, True, 0, 0, 1.5, 0),
                        (True, True, 0, 0, 1.75, 0),
                        (True, True, 0, 0, 2.0, 0),
                        (True, True, 0, 0, 2.25, 0),
                        (True, True, 0, 0, 2.50, 0),
                        (True, True, 0, 0, 2.75, 0),
                        (True, True, 0, 0, 3.0, 0)
                    ]
        
        for mismatch, credit, min_context, add_context, mult_context, \
            min_coding_cost in scheme_plan:
            
            # Start out with the context args
            extra_args = ["--context", str(min_context), "--addContext",
                str(add_context), "--multContext", str(mult_context), 
                "--minCodingCost", str(min_coding_cost)]
                
            # Make sure we have a Markov model
            if self.markov_model is None:
                raise Exception("We need a Markov model for coding cost!")
                
            # Send it along to the merger
            extra_args.append("--markovModel")
            extra_args.append(self.markov_model)
    
            # And give it a name to stick on our output files
            scheme_name = "Exact"
            if mismatch:
                # Add the args and scheme name component for mismatch
                extra_args.append("--mismatch")
                extra_args.append("--mismatches")
                extra_args.append("1")
                scheme_name = "Inexact"
            if credit:
                # Add the args and scheme name component for credit
                extra_args.append("--credit")
                scheme_name += "Credit"
                
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
                # Similarly for the min coding cost, replace the decimal.
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
        
        # This holds pairs of c2h and FASTA files
        c2h_fasta_pairs = []
        
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
            stats_filenames = [sonLib.bioio.getTempFile(suffix=".coverage",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's MAF alignment
            maf_filenames = [sonLib.bioio.getTempFile(suffix=".maf",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's hal alignment. We still need
            # individual hals because we need to have individual mafs for
            # comparison.
            hal_filenames = [sonLib.bioio.getTempFile(suffix=".hal",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's c2h alignment
            c2h_filenames = [sonLib.bioio.getTempFile(suffix=".c2h",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            # And another one to hold each child's fasta output for the c2h
            fasta_filenames = [sonLib.bioio.getTempFile(suffix=".fa",
                rootDir=self.getGlobalTempDir()) for i in xrange(num_children)]
                
            for hal in hal_filenames:
                # Make sure the HALs don't exist yet.
                os.unlink(hal)
                
            # Save these so we can compare them to other schemes later.
            alignments_by_scheme[scheme] = maf_filenames
            hals_by_scheme[scheme] = hal_filenames
                
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
                indel_filename = indel_filenames[i]
                tandem_filename = tandem_filenames[i]
                c2h_filename = c2h_filenames[i]
                fasta_filename = fasta_filenames[i]
                
                # Make a child to produce those, giving it a seed. Make sure to
                # give it only two FASTAs, reference first, so that when it
                # shuffles the non-reference ones it doesn't do anything. Also
                # make sure to tell it to use the spectrum output file.
                self.addChildTarget(ReferenceStructureTarget(
                    [reference_fasta, other_fasta], random.getrandbits(256), 
                    stats_filename, maf_filename, hal_filename=hal_filename,
                    spectrum_filename=spectrum_filename, 
                    indel_filename=indel_filename, 
                    tandem_filename=tandem_filename, c2h_filename=c2h_filename, 
                    fasta_filename=fasta_filename, extra_args=extra_args))
        
        
                
            # Make a follow-on job to merge all the child coverage outputs and
            # produce our coverage output file for this scheme.
            followOns.append(ConcatenateTarget(stats_filenames, self.stats_dir +
                "/coverage." + scheme))
                
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
            
        
        # What genome pairs did we run, in order? Make sure to strip extensions.
        genome_pairs = [(os.path.splitext(self.fasta_list[0])[0], 
            os.path.splitext(other)[0]) 
            for other in self.fasta_list[1:]]
                
        # Now we have all the followons for concatenating our stats and
        # comparing against the truth, we need a followon for comparing
        # corresponding MAFs from the same FASTA pair with different schemes,
        # for all combinations of schemes.
        agreement_target = AlignmentSchemeAgreementTarget(alignments_by_scheme,
            genome_pairs, random.getrandbits(256), 
            self.stats_dir + "/agreement", bed_root)
            
        # After that though, we need to take the BED files and the HAL files and
        # make assembly hubs in hub_root.
        hubs_target = AssemblyHubsTarget(hals_by_scheme, genome_pairs, bed_root,
            self.hub_root)
            
        # Do those two things in order.
        # TODO: These no longer work with unary-tree HALs.
        #followOns.append(SequenceTarget([agreement_target, hubs_target]))
        
        # We also want to combine all our c2h files into one massive hub.
        
        # What HAL will we use?
        merged_hal = sonLib.bioio.getTempFile(suffix=".hal",
                rootDir=self.getGlobalTempDir())
        os.unlink(merged_hal)
        
        # What merged c2h and fasta files will we use?
        merged_c2h = sonLib.bioio.getTempFile(suffix=".c2h",
                rootDir=self.getGlobalTempDir())
        merged_fasta = sonLib.bioio.getTempFile(suffix=".fa",
                rootDir=self.getGlobalTempDir())
        
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
            
        # And a target to make a hub from that
        merged_hub_target = AssemblyHubOnlyTarget(merged_hal, self.hub_root + 
            "/all")
            
        # Do those last two in order.
        followOns.append(SequenceTarget([merged_c2h_target, merged_hal_target,
            merged_hub_target]))
            
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
        
        elif len(self.c2h_fasta_pairs) == 1:
            # Merge the singleton against an empty file pair to add the prefix.
            
            # Make the empty pair
            c2h_empty = sonLib.bioio.getTempFile(suffix=".c2h",
                    rootDir=self.getGlobalTempDir())
            fasta_empty = sonLib.bioio.getTempFile(suffix=".fa",
                    rootDir=self.getGlobalTempDir())
            
            # Request the merge.
            self.addChildTarget(C2hMergeTarget(self.c2h_fasta_pairs + 
                [(c2h_empty, fasta_empty)], self.suffixes + [""], 
                self.merged_c2h, self.merged_fasta))
                
            # TODO: Do this whenever odd instead of at the leaves.
            
        elif len(self.c2h_fasta_pairs) == 2:
            
            print("Merging two pairs")
                
            # Merge the first two pairs
            check_call(self, ["../createIndex/cactusMerge", 
                self.c2h_fasta_pairs[0][0], self.c2h_fasta_pairs[0][1], 
                self.c2h_fasta_pairs[1][0], self.c2h_fasta_pairs[1][1],
                self.merged_c2h, self.merged_fasta, 
                "--suffix1", self.suffixes[0], "--suffix2", self.suffixes[1]])
        else:
            # Delegate to two children
            
            # What merged c2h and fasta files will we use for each?
            c2h_files = [sonLib.bioio.getTempFile(suffix=".c2h",
                    rootDir=self.getGlobalTempDir()) for i in xrange(2)]
            fasta_files = [sonLib.bioio.getTempFile(suffix=".fa",
                    rootDir=self.getGlobalTempDir()) for i in xrange(2)]
            
            # What should the first child get?
            split = len(self.c2h_fasta_pairs) / 2
            
            print("Splitting {} pairs at {}".format(
                len(self.c2h_fasta_pairs), split))
            
            # Tell it to go
            self.addChildTarget(C2hMergeTarget(self.c2h_fasta_pairs[0:split], 
                self.suffixes[0:split], c2h_files[0], fasta_files[0]))
                
            # And the second one
            self.addChildTarget(C2hMergeTarget(self.c2h_fasta_pairs[split:], 
                self.suffixes[split:], c2h_files[1], fasta_files[1]))
                
            # And merge their results when done, giving our result.
            self.setFollowOnTarget(C2hMergeTarget(zip(c2h_files, fasta_files), 
                ["", ""], self.merged_c2h, self.merged_fasta))
          
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
        
        # Make the base Target. Ask for 8gb of memory since this is hard.
        super(HalTarget, self).__init__(memory=8589934592)
        
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
                
        # Make the HAL.
        check_call(self, ["halAppendCactusSubtree", self.c2h_file, 
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
                
                if not os.path.exists(hub_dir):
                    # Don't let something else try making these in parallel.
                    os.makedirs(hub_dir)
                
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
        
        # Make the base Target. Ask for 8gb of memory since this is hard. And
        # ask for 32 CPUs so we can split up LOD.
        super(AssemblyHubTarget, self).__init__(memory=8589934592, cpu=32)
        
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
        
        self.logToMaster("Creating AssemblyHubTarget")
        
        
    def run(self):
        """
        Make the assembly hubs.
        """
        
        self.logToMaster("Starting AssemblyHubTarget")
        
                
        # Make a temporary directory for the jobTree tree. Make sure it's
        # absolute.
        tree_dir = os.path.abspath(sonLib.bioio.getTempDirectory(
            rootDir=self.getLocalTempDir()))
            
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
        bed_dirs = [self.bed_root + "/{}-{}-{}-{}".format(scheme1, scheme2,
            self.genome_pair[0], self.genome_pair[1]) 
            for (scheme1, scheme2) in possible_bed_pairs]
            
        # Keep the ones that exist and make a string of them. Make sure they are
        # absolute paths.
        bed_dirs_string = ",".join([os.path.abspath(directory) 
            for directory in bed_dirs if os.path.exists(directory)])
            
        if bed_dirs_string == "":
            bed_args = []
        else:
            # Only include the bedDirs option if there are beds.
            bed_args = ["--bedDirs", bed_dirs_string]
        
        # We want to make an assembly hub like so: 
                       
        # hal2assemblyHub.py data/alignment1.hal data/alignment1.hub
        # --jobTree data/tree --bedDirs data/beddir --hub=nocredit
        # --shortLabel="No Credit" --lod --cpHalFileToOut --noUcscNames
        
        # We need to do it in its own directory since it makes temp files in the
        # current directory.
        working_directory = sonLib.bioio.getTempDirectory(
            rootDir=self.getLocalTempDir())
            
        # Where should we come back to when done?
        original_directory = os.getcwd()
            
        # Turn our arguments into absolute paths (beddirs is already done).
        # May or may not be necessary.
        hal_abspath = os.path.abspath(self.hal)
        hub_abspath = os.path.abspath(self.hub)
        
        # Go in the temp directory and run the script
        os.chdir(working_directory)
        check_call(self, ["hal2assemblyHub.py", 
            hal_abspath, hub_abspath, "--jobTree", 
            tree_dir] + bed_args + ["--shortLabel", 
            self.scheme, "--lod", "--lodInMemory", "--lodNumProc=32",
            "--cpHalFileToOut", "--noUcscNames"])
        
        # Go back to the original directory
        os.chdir(original_directory)
            
        
        self.logToMaster("AssemblyHubTarget Finished")
        
class AssemblyHubOnlyTarget(jobTree.scriptTree.target.Target):
    """
    A target that makes a single assembly hub showing just a HAL.
    
    """
    
    def __init__(self, hal, hub):
        """
        Takes a hal alignment and makes a hub.
        
        """
        
        # Make the base Target. Ask for 8gb of memory since this is hard. And
        # ask for 32 CPUs so we can speed up LOD.
        super(AssemblyHubOnlyTarget, self).__init__(memory=8589934592, cpu=32)
        
        # Save the alignment file
        self.hal = hal
        
        # And the hub directory
        self.hub = hub
        
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
            
        # Turn our arguments into absolute paths (beddirs is already done).
        # May or may not be necessary.
        hal_abspath = os.path.abspath(self.hal)
        hub_abspath = os.path.abspath(self.hub)
        
        # Go in the temp directory and run the script
        os.chdir(working_directory)
        check_call(self, ["hal2assemblyHub.py", 
            hal_abspath, hub_abspath, "--jobTree", 
            tree_dir, "--lod", "--lodInMemory", "--lodNumProc=32",
            "--cpHalFileToOut", "--noUcscNames"])
        
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
        from compareSchemes import SchemeAssessmentTarget, \
            AlignmentSchemeAgreementTarget, BedSplitTarget, \
            AssemblyHubsTarget, AssemblyHubTarget
        
    # Make a stack of jobs to run
    stack = jobTree.scriptTree.stack.Stack(SchemeAssessmentTarget(
        options.fastas, options.trueMaf, options.markovModel, options.seed,
        options.outDir, options.outDir + "/hubs"))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
