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
    
    def __init__(self, fasta_list, true_maf, gene_bed_dir, seed, stats_dir,
    hub_root):
        """
        Make a new Target for building a several reference structures from the
        given FASTAs, and comparing against the given truth MAF, using the
        specified RNG seed, and writing statistics to one directory, and
        assembly hubs for different pairs of genomes and schemes to another.
        
        The BED files in gene_bed_dir, organized like
        gene_bed_dir/<genome>/<genome>.bed, will be used to evaluate alignment
        quality in light of gene annotations.
        
        The coverage statistics are just alignment coverage for each pair of
        genomes, in <genome number>\t<coverage fraction> TSVs per scheme.
        
        The alignment statistics are just <precision>\t<recall> TSVs per scheme.
        
        true_maf  may be None.
        
        Runs each subsequent FASTA against the first in each scheme.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(SchemeAssessmentTarget, self).__init__(memory=2147483648)
        
        # Save the FASTAs
        self.fasta_list = fasta_list
        
        # Save the filename of a MAF to compare all our MAFs against (or None)
        self.true_maf = true_maf
        
        # Save the gene bed directory to search
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
        Return a collection of tuples describing schemes.
        
        """
        
        # Plan out all the schemes as mismatch, credit, min_context,
        # add_context, mult_context, ignore_below, hamming_bound, hamming_max,
        # map_type
        return set([
            # Exact credit (tolerating 1 mismatch) with Hamming bound and
            # ignoring super short things.
            # Current best natural thing.
            (True, True, None, None, None, 10, 1, None, "natural"),
            (True, True, None, None, None, 10, 2, None, "natural"),
            (True, True, None, None, None, 10, 3, None, "natural"),
            (True, True, None, None, None, 10, 4, None, "natural"),
            (True, True, None, None, None, 10, 5, None, "natural"),
            (True, True, None, None, None, 10, 6, None, "natural"),
            (True, True, None, None, None, 10, 7, None, "natural"),
            (True, True, None, None, None, 10, 8, None, "natural"),
            # Exact credit (tolerating 1 mismatch) with Hamming bound and
            # Hamming distance allowance (i.e. mismatches again). Scheme under
            # test.
            (True, True, None, None, None, None, 2, 0, "natural"),
            (True, True, None, None, None, None, 2, 1, "natural"),
            (True, True, None, None, None, None, 3, 0, "natural"),
            (True, True, None, None, None, None, 3, 1, "natural"),
            (True, True, None, None, None, None, 3, 2, "natural"),
            (True, True, None, None, None, None, 4, 0, "natural"),
            (True, True, None, None, None, None, 4, 1, "natural"),
            (True, True, None, None, None, None, 4, 2, "natural"),
            (True, True, None, None, None, None, 4, 3, "natural")
        ])

    def generateSchemes(self):
        """
        Whatever schemes I want to test at the moment
        
        """
        
        # Get a big set of tuples succinctly describing all the schemes we want
        # to run.
        scheme_plan = self.getSchemePlan()
        
        for mismatch, credit, min_context, add_context, mult_context, \
            ignore_below, hamming_bound hamming_max, map_type in scheme_plan:
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
                
            # Yield the name with the args.
            yield (scheme_name, extra_args)
        
        
        
    def run(self):
        """
        Send off all the child targets to do comparisons.
        """
        
        self.logToMaster("Starting SchemeAssessmentTarget")
        
        # Unpack the genome names from the FASTA names. The first one is the
        # reference and the rest are queries.
        genome_names = [os.path.splitext(fasta)[0] for fasta in self.fasta_list]
        
        # Grab the reference and the queries
        reference_genome = genome_names[0]
        query_genomes = genome_names[1:]
        
        
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
            gene_beds = list(glob.glob(self.gene_bed_dir + "/*/*.bed"))
        else:
            # No gene BED files will be used.
            gene_beds = []
        
        # We need a few different follow-on jobs.
        followOns = []
        
        # We need a list of track targets that make tracks we want on the final
        # hub.
        track_targets = []
        
        # We need a directory to save out tree of BED files under for making the
        # assembly hubs.
        bed_root = sonLib.bioio.getTempDirectory(rootDir=self.getGlobalTempDir())
        
        # We need to keep track of all the bed dirs added under that. This is a
        # set so we can register a directory each time we add a genome's track.
        bed_dirs = set()
        
        # We have a wiggle root, but we make sure to keep that in our output so
        # the wiggles can be interrogated for context lengths.
        wiggle_root = self.stats_dir + "/wiggles"
        
        # And all the wiggle dirs. This is a set so we can register a directory
        # each time we add a genome's track, if we want to.
        wiggle_dirs = set()
        
        # We'll keep track of our random state so we can seed without messing up
        # temp filenames.
        random_state = None
        
        # This will hold lists of MAF alignments, in FASTA pair order, by scheme
        # name.
        alignments_by_scheme = {}
        
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
            
            if random_state is not None:
                # Pull out our unseeded random state so filenames on successive
                # iterations don't conflict.
                random.setstate(random_state)
            
            # Make a temp file for each of the children with this scheme to
            # write coverage stats to.
            coverage_filenames = [sonLib.bioio.getTempFile(suffix=".coverage",
                rootDir=self.getGlobalTempDir()) for _ in query_genomes]
                
            # And another one to hold each child's MAF alignment
            maf_filenames = [sonLib.bioio.getTempFile(suffix=".maf",
                rootDir=self.getGlobalTempDir()) for _ in query_genomes]
            
            # And another one to hold each child's c2h alignment
            c2h_filenames = [sonLib.bioio.getTempFile(suffix=".c2h",
                rootDir=self.getGlobalTempDir()) for _ in query_genomes]
                
            # And another one to hold each child's fasta output for the c2h
            fasta_filenames = [sonLib.bioio.getTempFile(suffix=".fa",
                rootDir=self.getGlobalTempDir()) for _ in query_genomes]
                
            # Save these so we can compare them to other schemes later.
            alignments_by_scheme[scheme] = maf_filenames
                
            # And another one to hold each child's adjacency component size
            # spectrum
            spectrum_filenames = [sonLib.bioio.getTempFile(suffix=".spectrum",
                rootDir=self.getGlobalTempDir()) for _ in query_genomes]
                
            # And another one to hold each child's indel lengths
            indel_filenames = [sonLib.bioio.getTempFile(suffix=".indels",
                rootDir=self.getGlobalTempDir()) for _ in query_genomes]
                
            # And another one to hold the tandem duplication counts.
            # TODO: these are probably just wrong at the moment.
            tandem_filenames = [sonLib.bioio.getTempFile(suffix=".tandem",
                rootDir=self.getGlobalTempDir()) for _ in query_genomes]
                
            # We also need files for the checkGenes category count outputs.
            checkgenes_count_filenames = [sonLib.bioio.getTempFile(
                suffix=".tsv", rootDir=self.getGlobalTempDir()) 
                for _ in query_genomes]
            
            # And for the checkGenes gene set output
            checkgenes_gene_filenames = [sonLib.bioio.getTempFile(
                suffix=".tsv", rootDir=self.getGlobalTempDir()) 
                for _ in query_genomes]
                
            # And we want a list of dicts of these things by mapping class.
            checkgenes_bed_filename_dicts = [
                {
                    classification: sonLib.bioio.getTempFile(
                        suffix=".bed", rootDir=self.getGlobalTempDir()) 
                    for classification in
                        # Skip the non2non mapping class because it's boring.
                        ["gene2gene", "gene2wrong", "gene2non", "non2gene"]
                }
                for _ in query_genomes
            ]
            
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
                    extra_args=extra_args))
                    
            
            # We need to check every MAF against the gene annotations. Out of
            # this we need to get mapping classification counts, the set of
            # query genes with any mappings in each category, and BED files for
            # each classification.
            
            # This list will hold the checkGenes targets that we need to run to
            # make all that.
            checkgenes_targets = []
            
            for maf, counts, genes, bed_dict in itertools.izip(maf_filenames, 
                checkgenes_count_filenames, checkgenes_gene_filenames,
                checkgenes_bed_filename_dicts):
                # For each set of parameters we need to run checkGenes on...
                
                # Make the target, and pass our global set of input gene BEDs.
                checkgenes_targets.append(MafGeneCheckerTarget(maf, gene_beds, 
                    counts, genes, bed_dict)) 
                    
            # Then for each BED that comes out of one of these checkGenes
            # targets, we need to BedSplitTarget it to rename the genome to
            # include the scheme, and to organize it in a proper bedDir for
            # passing to hal2assemblyHub.
            
            # This holds all the targets to do that
            bed_split_targets = []
            
            for query_genome, bed_dict in itertools.izip(query_genomes, 
                checkgenes_bed_filename_dicts):
                # For each genome and its corresponding set of BED files
                
                for classification, bed_file in bed_dict.iteritems():
                    # For each classification, we need a BedSplitTarget to split
                    # its BED and put it under the right bedDir for that
                    # classifications (so all the classifications become their
                    # own tracks). Make sure to re-name the query genome
                    # according to the scheme. TODO: Unify with the main suffix
                    # logic below.
                    bed_split_targets.append(BedSplitTarget(bed_file, 
                        [query_genome], bed_root + "/" + classification, 
                        genome_names=[query_genome + scheme]))
                        
                    # Make sure that the directory we are putting this bed in gets used
                    bed_dirs.add(bed_root + "/" + classification)
                
            
            # Check each MAF against all the genes, and then concatenate the
            # answers for this scheme.
            track_targets.append(SequenceTarget([
                # Run all the checkGenes runs
                RunTarget(checkgenes_targets),
                RunTarget([
                    # Concatenate the mapping classification counts 
                    ConcatenateTarget(checkgenes_count_filenames, 
                        self.stats_dir + "/checkgenes." + scheme),
                    # Concatenate the gene-with-mappings-in-category sets
                    ConcatenateTarget(checkgenes_gene_filenames, 
                        self.stats_dir + "/checkgenes-genesets." + scheme),
                    # Do all the BED splitting so we can have BED tracks on our
                    # hub
                    RunTarget(bed_split_targets)
                ])
            ]))
                
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
            pair_genomes += query_genomes
            # And record that each came from this scheme. Needs to keep the same
            # lenght as pair_genomes.
            pair_schemes += [scheme for _ in query_genomes]
            
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
                
        # What suffixes should we put on genomes?
        suffixes = [scheme for scheme in pair_schemes]
        
        # And what are the genomes with suffixes?
        genomes_with_suffixes = [genome + suffix for genome, suffix in 
            itertools.izip(pair_genomes, suffixes)]
            
        
        if self.gene_bed_dir is not None:
            # We want to go through this directory and duplicate all the genes
            # for each suffixed genome.
            
            # Make a bedDir for the genes
            genes_dir = bed_root + "/genes"
            os.mkdir(genes_dir)
            
            # If there are reference genes, where are they?
            reference_source = "{0}/{1}/{1}.bed".format(self.gene_bed_dir,
                reference_genome)
            
            if os.path.exists(reference_source):
                # Copy them over
                reference_destination = "{0}/{1}/{1}.bed".format(genes_dir,
                    reference_genome)
                    
                os.makedirs("{}/{}".format(genes_dir, reference_genome))
                    
                shutil.copyfile(reference_source, reference_destination)
            
            # Make sure we use these genes
            bed_dirs.add(genes_dir)
            
            for genome, suffix in itertools.izip(pair_genomes, suffixes):
                # For each genome and its suffix
                
                # Work out what it needs to be called.
                suffixed = genome + suffix
                
                # Where would the genes come from?
                genome_source = "{0}/{1}/{1}.bed".format(self.gene_bed_dir,
                    genome)
                    
                if os.path.exists(genome_source):
                    # They exist.
                    
                    # Where should the gene results go?
                    gene_out_dir = "{}/{}".format(genes_dir, suffixed)
                    if not os.path.exists(gene_out_dir):
                        # Make the directory now
                        os.makedirs(gene_out_dir) 
                    
                    # We need to change the contig name on every gene, which we
                    # can do with a BedSplitTarget. It will automatically decide
                    # to use a BED named after the new genome name.
                    track_targets.append(BedSplitTarget(
                        genome_source, [genome], genes_dir, 
                        genome_names=[suffixed]))
                        
                    # Report what we're up to.
                    self.logToMaster("Taking genes from {} for {}".format(
                        genome, suffixed))
        
        # Make a target to make the merged c2h/fasta files
        merged_c2h_target = C2hMergeTarget(c2h_fasta_pairs, suffixes, 
            merged_c2h, merged_fasta)
            
        # And a target to make the hal from them
        merged_hal_target = HalTarget(merged_c2h, merged_fasta, 
            reference_genome, genomes_with_suffixes, merged_hal)
            
        # And a target to make a hub from that. Make sure to pass all the track
        # directories.
        merged_hub_target = AssemblyHubOnlyTarget(merged_hal, self.hub_root + 
            "/all",
                # TODO: Disabling liftovers of these things for faster hub
                # generation. Make this a flag, preferably by type of track
                # (genes vs. gene alignment evaluation vs. debug context length
                # tracks vs. what have you).
                # bed_dirs=list(bed_dirs), 
                # wiggle_dirs=list(wiggle_dirs)
            )
            
        # Do those last ones in order. TODO: express actual order restrictions
        # with a few more intermediate targets.
        followOns.append(SequenceTarget([
            RunTarget([
                # We can do all these tracks in parallel with the HAL.
                RunTarget(track_targets),
                SequenceTarget([
                    # We need to merge the c2h files before the HAL.
                    merged_c2h_target, 
                    merged_hal_target
                ])
            ]), 
            # When the tracks are done, we can make the hub
            merged_hub_target
        ]))
            
        # So we need to run all of the follow-ons in parallel after this job
        self.setFollowOnTarget(RunTarget(followOns))
                
        self.logToMaster("SchemeAssessmentTarget Finished")

class MafGeneCheckerTarget(jobTree.scriptTree.target.Target):
    """
    A target that checks a MAF alignment against a set of BED files, and reports
    the number of mappings in various categories of correctness.
    
    """
    
    def __init__(self, maf_file, gene_bed_files, class_count_file, 
        gene_set_file, class_beds):
        """
        Reads the maf_file and the genes from the gene_bed_files list.
        Classifies every mapping in the alignment. Saves class counts to
        class_count_file ans a <class>\t<count> TSV. Saves genes with any
        mappings of each class in them to gene_set_file as a <class>\t<gene
        name> TSV.
        
        For each classification to BED filename mapping in class_beds, saves a
        BED file of mapped locations that are placed in that class.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(MafGeneCheckerTarget, self).__init__(memory=2147483648)
        
        # Save the arguments
        self.maf_file = maf_file
        self.gene_bed_files = gene_bed_files
        self.class_count_file = class_count_file
        self.gene_set_file = gene_set_file
        self.class_beds = class_beds
        
        self.logToMaster("Creating MafGeneCheckerTarget")
        
        
    def run(self):
        """
        Run the checks.
        """
        
        self.logToMaster("Starting MafGeneCheckerTarget")
        
        # Prepare arguments
        args = (["./checkGenes.py", "--maf", self.maf_file, "--beds"] + 
            self.gene_bed_files + ["--classCounts", self.class_count_file,
            "--geneSets", self.gene_set_file])
            
        for classification, bed_file in self.class_beds.iteritems():
            # We need to pass each mapping in this dict along. Since the options
            # are predictably named, we generate them.
            args.append("--" + classification + "Bed")
            args.append(bed_file)
        
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
        
        # Make the base Target. Ask for 32gb of memory since the intermediate
        # graphs can get big.
        super(C2hMergeTarget, self).__init__(memory=34359738368)
        
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
    A target which splits a BED up by genome. Also fixes them up a bit by
    collapsing adjacent identical features.
    
    """

    def __init__(self, bed_file, genomes, out_dir, feature_name=None,
        genome_names=None):
        """
        Split the given bed file per genome into <genome>/<genome>.bed in the
        given output directory, using the given list of genomes.
        
        If out_dir does not exist, it will be created.
        
        If a feature_name is specified, it is used to rename all the BED
        feratures.
        
        If genome_names is specified, each genome in genomes is renamed to the
        corresponding entry in new_names.
        
        Does not leave a trailing newline.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(BedSplitTarget, self).__init__(memory=2147483648)
        
        # Save the arguments
        self.bed_file = bed_file
        self.genomes = genomes
        self.out_dir = out_dir
        self.feature_name = feature_name
        
        if genome_names is not None:
            # Build a dict with which we can efficiently rename genomes.
            self.genome_map = {old_name: new_name 
                for old_name, new_name in itertools.izip(genomes, genome_names)}
        else:
            self.genome_map = None
    
    def get_name(self, genome_name):
        """
        Given the name of a genome, return the name we should use for it,
        subject to any renaming we are supposed to do.
        
        """
        
        if self.genome_map is not None:
            # Apply the renaming map
            return self.genome_map[genome_name]
        else:
            # Don't rename since there is no renaming map.
            return genome_name
            
    def run(self):
        """
        Run this target and do the splitting.
        
        """
    
        self.logToMaster("Starting BedSplitTarget")
        
        if not os.path.exists(self.out_dir):
            try:
                # Make sure out_dir exists
                os.makedirs(self.out_dir)
            except OSError:
                self.logToMaster("OSError trying to make {}. Did it get " +
                    "created after we looked for it?".format(self.out_dir))
        
        # Work out output file names
        out_filenames = [self.out_dir + "/{0}/{0}.bed".format(
            self.get_name(genome)) for genome in self.genomes]
            
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
        
        # Count lines processed
        lines_processed = 0
            
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
                
            # Say we're doing this line
            lines_processed += 1
            
            for genome, temp_file in itertools.izip(self.genomes, temp_files):
                if parts[0] == genome:
                    # Send the line to the file corresponding to the matching
                    # genome. TODO: use a dict or something instead of scanning
                    # for matching names.
                    
                    # Rename the genome if needed
                    parts[0] = self.get_name(genome)
                    
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
                
                if not os.path.exists(self.out_dir + "/" + 
                    self.get_name(genome)):
                    # Make the directory for the final BED file.
                    os.mkdir(self.out_dir + "/" + self.get_name(genome))
            
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
        
        
        self.logToMaster("BedSplitTarget Finished ({} lines)".format(
            lines_processed))
        
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
            try:
                # Make a directory for each genome with its suffix
                os.makedirs(self.out_dir + "/" + suffixed)
            except:
                # TODO: make sure this failed because the directory existed,
                # without getting an error if its parent doesn't exist.
                pass
        
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
    
    def __init__(self, hal, hub, bed_dirs=None, wiggle_dirs=None):
        """
        Takes a hal alignment and makes a hub.
        
        Can also optionally take a directory of BEDs and a directory of wiggles.
        Each is named after the name of the track, and arranged like:
        <track name>/<genome>/<genome>.<extension>
        
        """
        
        # Make the base Target. Ask for 8gb of memory and several CPUs since
        # this is hard.
        super(AssemblyHubOnlyTarget, self).__init__(memory=8589934592, cpu=32)
        
        # Save the alignment file
        self.hal = hal
        
        # And the hub directory
        self.hub = hub
        
        # And the BED directory list
        self.bed_dirs = bed_dirs        
        
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
        
        if self.bed_dirs is not None:
            # Make bed directories absolute paths if needed.
            self.bed_dirs = [os.path.abspath(path) 
                for path in self.bed_dirs]
        
        if self.wiggle_dirs is not None:
            # Make wiggle directories absolute paths if needed.
            self.wiggle_dirs = [os.path.abspath(path) 
                for path in self.wiggle_dirs]
        
        # Assemble the command
        command = ["hal2assemblyHub.py", 
            hal_abspath, hub_abspath, "--jobTree", 
            tree_dir, "--maxThreads", "32", # No LOD
            "--cpHalFileToOut", "--noUcscNames"]
            
        if self.bed_dirs is not None:
            # Add the option to use these bed dirs. TODO: these directories
            # can't have commas in their names
            command.append("--bedDirs")
            command.append(",".join(self.bed_dirs))
            
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
    parser.add_argument("--geneBedDir", default=None,
        help="hal2assemblyHub.py --bedDirs entry for evaluating alignments")
        
    
    
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
        options.fastas, options.trueMaf, options.geneBedDir, options.seed,
        options.outDir, options.outDir + "/hubs"))
    
    print "Starting stack"
    
    # Run it and see how many jobs fail
    failed_jobs = stack.startJobTree(options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print "All jobs completed successfully"
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
