#!/usr/bin/env python2.7
"""
targets.py: a library of useful jobTree targets for working with createIndex and
reference structures.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import datetime
from xml.etree import ElementTree
import tsv

import jobTree.scriptTree.target
import jobTree.scriptTree.stack
import sonLib.bioio

def check_call(target, args):
    """
    Make a subprocess call, announcing that we are doing it.
    """
    print("Calling: {}".format(" ".join(args)))
    sys.stdout.flush()
    target.logToMaster("Calling: {}".format(" ".join(args)))
    return subprocess.check_call(args)

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
        
class SequenceTarget(jobTree.scriptTree.target.Target):
    """
    A target that runs a list of targets in succession. This is kind of a hack
    around the way JobTree wants you to think, but I need it to keep my targets
    from becoming dependent on what needs to use their results.
    
    """
    
    def __init__(self, targets):
        """
        Maker a new target which runs the targets in the list in the order they
        appear in.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(SequenceTarget, self).__init__(memory=2147483648)
        
        if len(targets) > 0:
            # If we have a first target, split that out.
            self.first = targets[0]
        else:
            self.first = None
        
        if len(targets) > 1:
            # If we have more than that, wrap them up in a linked list of
            # SequenceTargets.
            self.next = SequenceTarget(targets[1:])
        else:
            self.next = None
            
    def run(self):
        """
        Run this target and all subsequent ones.
        
        """
    
        self.logToMaster("Starting SequenceTarget")
    
        if self.first is not None:
            # Run the head of the linked list.
            self.addChildTarget(self.first)
            
        if self.next is not None:
            # Continue on with the rest of it.
            self.setFollowOnTarget(self.next)
        
        self.logToMaster("SequenceTarget Finished")
        
class SchemeUsingTarget(jobTree.scriptTree.target.Target):
    """
    A base class for targets that do something with multiple mapping/merging
    schemes. Contains utility functions to turn "scheme plans", which are lists
    of parameter tuples, into scheme names and extra_args values for use with a
    ReferenceStructureTarget
    
    """
    
    def __init__(self, *args, **kwargs):
        """
        Forward all the init arguments on to the base jobTree constructor.
        
        """
        # See <http://stackoverflow.com/a/6535962/402891>
        super(SchemeUsingTarget, self).__init__(*args, **kwargs)
        
    def getSchemePlan(self):
        """
        Return a collection of tuples describing schemes.
        
        They all start with map_type. For "natural" schemes, we have map_type,
        min context, credit, mismatches for credit, min edit distance bound, max
        edit distance. For "zip" schemes we have map_type, min context, max
        range count.
        
        This should be overridden in subclasses.
        """
        
        return []
    
    def generateSchemes(self):
        """
        Yield tuples of scheme names and the extra arguments createIndex needs
        to use those mapping schemes.
        
        This will automatically flesh out the schemes defined in getSchemePlan.
        """
        
        # Get a big set of tuples succinctly describing all the schemes we want
        # to run.
        scheme_plan = self.getSchemePlan()
        
        for plan_item in scheme_plan:
            # Pull out each planned scheme and see what map type it is using,
            # which determines its parameters.
            map_type = plan_item[0]
            
            # Start out with no configuration arguments
            extra_args = []
            
            # Give the scheme a name to stick on our output files
            scheme_name = "E"
            
            if map_type == "natural":
                # Unpack the parameters.
                map_type, min_context, credit, mismatch, hamming_bound, \
                    hamming_max = plan_item
                    
                
                
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
                    
            elif map_type == "zip":
                # Unpack the parameters
                map_type, min_context, range_count = plan_item
                
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
                
            if map_type == "natural":
            
                if hamming_bound is not None:
                    # For the natural mapping scheme, don't map on a maximum
                    # unique match run unless we can get a lower bound on its
                    # Hamming distance from all other reference locations that
                    # is at least this high.
                    extra_args.append("--minHammingBound")
                    extra_args.append(str(hamming_bound))
                    
                    # Mention it in the scheme name
                    scheme_name += "Ham{}".format(hamming_bound)
                    
                if hamming_max is not None:
                    # For the natural mapping scheme, allow this many mismatches
                    # in maximal unique match runs.
                    extra_args.append("--maxHammingDistance")
                    extra_args.append(str(hamming_max))
                    
                    # Mention it in the scheme name
                    scheme_name += "Mis{}".format(hamming_max)
                    
            if map_type == "zip":
                
                if range_count is not None:
                    # Handle a limit on the max number of ranges to visit
                    extra_args.append("--maxRangeCount")
                    extra_args.append(str(range_count))
                    
                    # Mention it in the scheme name
                    scheme_name += "Range{}".format(hamming_max)
        
            # Yield the name with the args.
            yield (scheme_name, extra_args)
    
    

class ReferenceStructureTarget(jobTree.scriptTree.target.Target):
    """
    A target that builds a reference structure.
    
    """
    
    def __init__(self, fasta_list, seed, coverage_filename, alignment_filename,
        hal_filename=None, spectrum_filename=None, indel_filename=None,
        tandem_filename=None, c2h_filename=None, fasta_filename=None, 
        extra_args=[]):
        """
        Make a new Target for building a reference structure from the given
        FASTAs, using the specified RNG seed, and writing coverage statistics to
        the specified file, and the alignment MAF to the other specified file.
        
        Each FASTA file must contain exactly one contig, named after the FASTA
        without the extension.
        
        Those coverage statistics are, specifically, alignment coverage of a
        genome vs. the order number at which the genome is added, and they are
        saved in a <genome number>\t<coverage fraction> TSV.
        
        If hal_filename is specified, the HAL file will be written there instead
        of a temporary file.
        
        If spectrum_filename is specified, saves the adjacency component size
        spectrum to the given file, as a TSV of <size>\t<count> lines.
        
        If indel_filename is specified, save the lengths of all simple (between
        two degree-two blocks) indels to that file.
        
        If tandem_filename is specified, save the count of tandem duplications
        detected (4-end rearrangements that involve two connected ends of the
        same block) to that file.
        
        If c2h_filename is specified, the c2h intermediate alignment file is
        saved to that location.
        
        If fasta_filename is specified, the FASTA that goes with the c2h file is
        saved there.
        
        If extra_args is specified, it should be a list of additional command-
        line arguments to createIndex, for specifying things like inexact
        matching.
        
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
        
        # And the (MAF) alignemnt filename to use
        self.alignment_filename = alignment_filename
        
        # And the HAL filename to use, if any
        self.hal_filename = hal_filename
        
        # And the spectrum filename to use, if any
        self.spectrum_filename = spectrum_filename
        
        # And the indel lengths filename, if any
        self.indel_filename = indel_filename
        
        # And the tandem duplication count filename, if any
        self.tandem_filename = tandem_filename
        
        # And the cactus2hal file for the alignment to be written to
        self.c2h_filename = c2h_filename
        
        # And the FASTA to go with it
        self.fasta_filename = fasta_filename
        
        # And the extra args
        self.extra_args = extra_args
        
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
        # Make sure the first one stays first.
        fasta_first = self.fasta_list[0]
        fasta_rest = list(self.fasta_list[1:])
        random.shuffle(fasta_rest)
        self.fasta_list = [fasta_first] + fasta_rest
        
        # Make a temp directory for the index
        index_dir = sonLib.bioio.getTempFile(rootDir=self.getGlobalTempDir())
        
        # Make a file for the cactus2hal, if not specified
        c2h_filename = self.c2h_filename or sonLib.bioio.getTempFile(
            rootDir=self.getLocalTempDir())
        
        # Make a file for the FASTA that we need to use the cactus2hal, if not
        # specified
        fasta_filename = self.fasta_filename or sonLib.bioio.getTempFile(
            rootDir=self.getLocalTempDir())
        
        # Put together the arguments to invoke
        args = ["../createIndex/createIndex", "--scheme", "greedy", 
            "--alignment", c2h_filename, "--alignmentFasta", fasta_filename, 
            index_dir] + self.fasta_list + self.extra_args
            
        if self.spectrum_filename is not None:
            # We want to keep the adjacency spectrum.
            args.append("--spectrum")
            args.append(self.spectrum_filename)
            
        if self.indel_filename is not None:
            # We want to keep the indel lengths.
            args.append("--indelLengths")
            args.append(self.indel_filename)
            
        if self.tandem_filename is not None:
            # We want to keep the tandem duplication count.
            args.append("--tandemDuplications")
            args.append(self.tandem_filename)
            
        # Announce our command we're going to run
        self.logToMaster("Invoking {}".format(" ".join(args)))
        print("Invoking {}".format(" ".join(args)))
        
        # Make an index with the FASTAs in that order, so we can read all the
        # logging output. Make sure to send any errors through stdout as well.
        process = subprocess.Popen(args,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        # Make a writer for the statistics
        writer = tsv.TsvWriter(open(self.coverage_filename, "w"))
    
        # Track how many lines of indexer output we process
        lines = 0
    
        # Tell our jobtree output we're reading from the process starting now,
        # in case something goes wrong.
        print("[{}] Reading process output...".format(datetime.datetime.now()))
    
        for line in process.stdout:            
            # Collect and parse the log output, and get the coverage vs. genome
            # number data.
            
            if "DEBUG:" not in line and "TRACE:" not in line:
                # Log things that came out at high logging priorities (TODO: or
                # are false positives) to our own log.
                print(line.strip())
            
            # Record that we processed a line
            lines += 1
            
            if "Coverage from alignment of genome" in line:
                # Grab the genome number
                genome = int(line.split(":")[-2].split("genome")[1])
                # Grab the coverage fraction
                coverage = float(line.split("=")[1])
            
                # Save the coverage vs. genome number data as a reasonable TSV.
                writer.line(genome, coverage)
                
        print("[{}] Process output done".format(datetime.datetime.now()))
                
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
        
        # What genomes do we have? Make a list of them in the order they were
        # added. TODO: This tightly couples us to both the genome name choosing
        # logic and the particular test case we are working on. Fix that.
        genomes = [os.path.splitext(os.path.split(fasta)[1])[0] 
            for fasta in self.fasta_list]
        
        # What tree should we use? Assume the genome names are the same as the
        # FASTA names without their extensions. Set up for a star tree rooted at
        # "rootSeq".
        tree = "(" + ",".join(genomes) + ")rootSeq;"
        
        if self.hal_filename is None:
            # Where should we save it? ("" isn't a filename so this or is OK)
            hal_filename = sonLib.bioio.getTempFile(
                rootDir=self.getLocalTempDir())
        else:
            hal_filename = self.hal_filename
        
        self.logToMaster("Creating HAL {} with tree: {}".format(hal_filename,
            tree))
        
        if(os.path.exists(hal_filename)):
            # Make sure HAL doesn't exist when we try to make it
            os.unlink(hal_filename)
        
        # Turn the c2h into a HAL with the given tree.
        check_call(self, ["halAppendCactusSubtree", c2h_filename, 
            fasta_filename, tree, hal_filename])
            
        self.logToMaster("Creating MAF")
            
        # Now we need to make that HAL into a MAF, saving to the file where
        # we're supposed to send out output. Use the first genome that we added
        # as the reference (since MAF can't represent a duplication in the
        # reference, and we will never call a duplication relative to this first
        # genome), and don't include rootSeq (since it is just an artifact of
        # the HAL format and not a real genome).
        check_call(self, ["hal2maf", "--refGenome", genomes[0], 
            "--targetGenomes", ",".join(genomes[1:]), hal_filename,
            self.alignment_filename])
            
        self.logToMaster("ReferenceStructureTarget Finished")
        
class AlignmentComparisonTarget(jobTree.scriptTree.target.Target):
    """
    A target that compares two MAF alignments.
    
    """
    
    def __init__(self, maf_a, maf_b, seed, output_filename, is_correct=False,
        bed=None):
        """
        Compare the two MAFs referred to by the given input filenames, and write
        a digest of mafComparator results to the given output filename. Uses a
        random seed to generate the mafComparator seed.
        
        If is_correct is true, the first alignment is considered ground truth,
        and a two column <precision>\t<recall> output file is produced.
        
        If bed is specified, mafComnparator's --dumpBed option is used to dump a
        BED of mismatching positions.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AlignmentComparisonTarget, self).__init__(memory=2147483648)
        
        # Save the parameters
        self.maf_a = maf_a
        self.maf_b = maf_b
        self.seed = seed
        self.output_filename = output_filename
        self.is_correct = is_correct
        self.bed = bed
        
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

        # Set up the mafComparator arguments. Make sure to ask for an absurd
        # number of samples so we check everything.
        args = ["mafComparator", "--maf1", self.maf_a, "--maf2",
            self.maf_b, "--out", xml_filename, "--seed", str(seed), "--samples",
            "100000000"]
        
        if(self.bed is not None):
            args.append("--dumpMismatches")
            args.append(self.bed)
        
        # Run mafComparator and generate the output XML
        check_call(self, args)
        
        # Now parse the XML output down to the statistics we actually want.
        tree = ElementTree.parse(xml_filename)
        
        # Grab the nodes with the aggregate results for each direction
        stats = tree.findall(
            "./homologyTests/aggregateResults/all")
            
        # Grab and parse the averages. There should be two.
        averages = [float(stat.attrib["average"]) for stat in stats]
        
        # Open the output file to save them to.
        writer = tsv.TsvWriter(open(self.output_filename, "w"))
        
        if self.is_correct:
            # Save precision and recall. The first average from mafComparator is
            # for recall (homologies in the first file that appear in the
            # second), and the second average from mafComparator is for
            # precision (homologies in the second file that appear in the
            # first). We save them as (precision, recall).
            writer.line(averages[1], averages[0])
        else:
            # Calculate the F score and save that, since it's symmetric and
            # doesn't care which of the two inputs is "truth".
            writer.line(2 * averages[0] * averages[1] / 
                (averages[0] + averages[1]))
            
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
        
class AlignmentTruthComparisonTarget(jobTree.scriptTree.target.Target):
    """
    A target that compares a set of MAF alignments against a single MAF
    alignment truth.
    
    """
    
    def __init__(self, truth_maf, maf_filenames, seed, output_filename):
        """
        Compare the given truth MAF filename to the MAFs referred to by the
        given input filenames, and write collated results to the given output
        filename. Uses a random seed to generate mafComparator seeds.
        
        """
        
        # Make the base Target. Ask for 2gb of memory since this is easy.
        super(AlignmentTruthComparisonTarget, self).__init__(memory=2147483648)
        
        # Save the parameters
        self.truth_maf = truth_maf
        self.maf_filenames = maf_filenames
        self.seed = seed
        self.output_filename = output_filename
        
        self.logToMaster("Creating AlignmentTruthComparisonTarget")
        
        
    def run(self):
        """
        Send off all the child targets to do alignment comparisons, and
        integrate their results.
        
        """
        
        self.logToMaster("Starting AlignmentTruthComparisonTarget")
        
        # Make a temp file for each of the children to write stats to
        comparison_filenames = [sonLib.bioio.getTempFile(
            rootDir=self.getGlobalTempDir()) for other in self.maf_filenames]
        
        # Seed the RNG after sonLib does whatever it wants with temp file names
        random.seed(self.seed)
        
        for maf_filename, comparison_filename in zip(self.maf_filenames,
            comparison_filenames):
            
            # Make a child to compare this MAF against the truth, and produce
            # this output file, giving it a 256-bit seed.
            self.addChildTarget(AlignmentComparisonTarget(self.truth_maf, 
                maf_filename, random.getrandbits(256), comparison_filename, 
                is_correct=True))
            
        # When we're done, concatenate all those comparison files together into
        # our output file.
        self.setFollowOnTarget(ConcatenateTarget(comparison_filenames,
            self.output_filename))
        
        self.logToMaster("AlignmentTruthComparisonTarget Finished")
        
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
    each is a binary tree of the given reference and some other genome, against
    a root. The corresponding suffix from the list of suffixes is assigned to
    each genome. Saves a merged c2h/fasta pair.
    
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
        # TODO: Make genome to merge on specifiable!
        args = ["../createIndex/cactusMerge", "refmhc", self.merged_c2h,
            self.merged_fasta]
        
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
        
        # Compose the tree. We know we are using the fake root sequence
        # alignment method.
        tree = "(" + ",".join(self.genomes + [self.reference]) + ")rootSeq;"
                
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
        
