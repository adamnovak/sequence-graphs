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

class SchemeAssessmentTarget(SchemeUsingTarget):
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
        
        true_maf may be None.
        
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
        
        # Make sure we have an RNG
        self.rng = random.Random()
        
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
        
        # Plan out all the schemes. They all start with map_type. For "natural"
        # schemes, we have map_type, min context, credit, mismatches for credit,
        # min edit distance bound, max edit distance. For "zip" schemes we have
        # map_type, min context, max range count.
        return set([
            # Natural with flat min thresholds
            ("natural", 20, False, False, None, None),
            ("natural", 50, False, False, None, None),
            ("natural", 100, False, False, None, None),
            ("natural", 150, False, False, None, None),
            ("natural", 200, False, False, None, None),
            # Do a 5-4-natural scheme with credit. 
            ("natural", None, True, True, 5, 4),
            # Do zip with flat mins
            ("zip", 20, 100),
            ("zip", 50, 100),
            ("zip", 100, 100),
            ("zip", 150, 100),
            ("zip", 200, 100)
        ])

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
            self.rng.seed(self.seed)
            
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
                    [reference_fasta, other_fasta], self.rng.getrandbits(256),
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
                maf_filenames, self.rng.getrandbits(256), 
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
