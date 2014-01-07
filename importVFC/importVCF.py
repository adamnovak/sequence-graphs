#!/usr/bin/env python2.7
"""
importVCF.py: convert a VCF file to an Avro-based sequence graph format.

This script reads a VCF file containing simple SNP and indels (but no
rearrangements or structural variants) for a single sample. It produces a
sequence graph in two Avro-format files of Avro records: one holding
AlleleGroups, and one holding Adjacencies between AlleleGroups. Phasing
information from the VCF is preserved (at least between adjacent sites; if the
first chromosome in the VCF is always maternal and the second always paternal,
that information will be lost).

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, collections, math, itertools, logging

import avro.schema, avro.datafile, avro.io
import vcf

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
    parser.add_argument("vcf", type=argparse.FileType("r"),
        help="VCF file to read")
    parser.add_argument("allele_group_file", type=argparse.FileType("w"),
        help="output file for allele groups")
    parser.add_argument("adjacency_file", type=argparse.FileType("w"),
        help="output file for adjacencies")
    parser.add_argument("--allele_group_schema", 
        default="schemas/AlleleGroup.avsc", 
        help="Avro schema for AlleleGroups")
    parser.add_argument("--adjacency_schema", 
        default="schemas/Adjacency.avsc", 
        help="Avro schema for Adjacencies")
        
    # Logging options
    parser.add_argument("--loglevel", default="DEBUG", choices=["DEBUG", "INFO",
        "WARNING", "ERROR", "CRITICAL"],
        help="logging level to use")
        
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def set_loglevel(loglevel):
    """
    Given a string log level name, set the logging module's log level to that
    level. Raise an exception if the log level is invalid.
    
    """
    
    # Borrows heavily from (actually is a straight copy of) the recipe at
    # <http://docs.python.org/2/howto/logging.html>
    
    # What's the log level number from the module (or None if not found)?
    numeric_level = getattr(logging, loglevel.upper(), None)
    
    if not isinstance(numeric_level, int):
        # Turns out not to be a valid log level
        raise ValueError("Invalid log level: {}".format(loglevel))
    
    # Set the log level to the correct numeric value
    logging.basicConfig(format="%(levelname)s:%(message)s", level=numeric_level)

class AlleleGroup(dict):
    """
    Make a new AlleleGroup, which is a dict with some structure to it.
    
    An AlleleGroup must have the following keys:
        
        "id": a unique string ID tof the AlleleGroup.
        
        "fivePrime": a unique string ID for the 5' end of the AlleleGroup.
        
        "threePrime": a unique string ID for the 3' end of the AlleleGroup.
        
        "contig": a string naming the reference (or auxilliary) contig this
        AlleleGroup is located on.
        
        "start": an integer expressing wherte the AlleleGroup starts on the
        contig (0-based, inclusive).
        
        "end": an integer expressing wherte the AlleleGroup ends on the contig
        (0-based, exclusive).
        
    An AlleleGroup may have the following optional keys:
        
        "sequence": a string representing the DNA sequence for this AlleleGroup,
        if it differs from that of the reference. May be None. May also hold the
        reference sequence.
        
        "ploidy": an integer describing the number of copies of this AlleleGroup
        that are present. May be None.
        
        "sample": a string identifying what sample or individual this
        AlleleGroup belongs to. May be None.
              
    AlleleGroups can be hapily serialized to Avro records. The Avro records will
    deserialize back to normal dicts, though, so don't rely on this class having
    any attributes or methods; it's just here to provide a fancy constructor.
    
    """
    
    # Class-level things for generating IDs
    next_id = 0
    @classmethod
    def get_id(cls):
        """
        Get the next available unique string ID for AlleleGroups.
        
        """
        
        # Grab the next ID string
        to_return = "ag-{}".format(cls.next_id)
        
        # Advance the next ID
        cls.next_id += 1
        
        # Return the generated ID
        return to_return
    
    def __init__(self, contig, start, end):
        """
        Make a new AlleleGroup representing the region defined by the given
        reference coordinates.
        
        """
        
        # Set up unique IDs
        self["id"] = AlleleGroup.get_id()
        self["fivePrime"] = "{}-5'".format(self["id"])
        self["threePrime"] = "{}-3'".format(self["id"])
        
        # Save the reference coordinates
        self["contig"] = contig
        self["start"] = start
        self["end"] = end
        
class Adjacency(dict):
    """
    Make a new Adjacency, which is a dict with some structure to it.
    
    An Adjacency must have the following keys:
        
        "id": a unique string ID tof the Adjacency.
        
        "first": a reference to a "threePrime" or "fivePrime" ID from an
        AlleleGroup, from which the adjacency comes.
        
        "second": a reference to a "threePrime" or "fivePrime" ID from an
        AlleleGroup, to which the adjacency goes.
        
    An AlleleGroup may have the following optional keys:
        
        "ploidy": an integer describing the number of copies of this Adjacency
        that are present. May be None.
        
        "sample": a string identifying what sample or individual this Adjacency
        belongs to. May be None.
              
    Adjacencies can be hapily serialized to Avro records. The Avro records will
    deserialize back to normal dicts, though, so don't rely on this class having
    any attributes or methods; it's just here to provide a fancy constructor.
    
    """
    
    # Class-level things for generating IDs
    next_id = 0
    @classmethod
    def get_id(cls):
        """
        Get the next available unique string ID for AlleleGroups.
        
        """
        
        # Grab the next ID string
        to_return = "adj-{}".format(cls.next_id)
        
        # Advance the next ID
        cls.next_id += 1
        
        # Return the generated ID
        return to_return
    
    def __init__(self, first, second):
        """
        Make a new Adjacency representing an adjacency between the two endpoint
        IDs in the given order.
        
        """
        
        # Set up unique IDs
        self["id"] = Adjacency.get_id()
        
        # Save the AlleleGroup ends we tie together
        self["first"] = first
        self["second"] = second
        
def import_sample(vcf_records, sample_name, allele_group_writer,
    adjacency_writer):
    """
    Given a list of VCF records sorted in contig/coordinate order, a sample
    name, an Avro allele group writer, and an Avro adjacency writer, create Avro
    sequence graph records for the sample with the given name, describing its
    variation, and write those records to the writers.
    
    """
    
    logging.info("Importing sample {}".format(sample_name))
    
    # We would like to do this without all sorts of special cases for phasing.
    
    # We keep track of the last 3' ID on each of the two chromosomes. They are
    # both None when we start a chromosome, and both the same when we're going
    # through an area with no phasing.
    last_ends = [None, None]
    
    # We have some functions to manipulate these.
    def add_allele_group(allele_group, chromosome):
        """
        Take the given AlleleGroup, and create adjacencies to append it to the
        given chromosome (0 or 1). Updates the last end of that chromosome.
        
        """
        
        logging.debug("Appending {} onto chromosome {}".format(
            allele_group["id"], chromosome))
        
        if last_ends[chromosome] is not None:
            # We have something to append to
        
            # Make the adjacency tying the last 3' end to the new thing's 5'
            # end.
            adjacency = Adjacency(last_ends[chromosome],
                allele_group["fivePrime"])
                            
            # Set its ploidy to 1. If we append the same thing on both
            # chromosomes, and both chromosomes were unphased to begin with, we
            # end up with 2 ploidy-1 edges.
            adjacency["ploidy"] = 1
            
            # Tag with the sample's ID
            adjacency["sample"] = sample_name
            
            # Save the adjacency
            adjacency_writer.append(adjacency)
        
        # Update the chromosome end
        last_ends[chromosome] = allele_group["threePrime"]
        
    def new_contig():
        """
        Start a new contig. Clear the saved last endpoints.
        
        """
        
        last_ends = [None, None]
        
    # We keep track of the last call (sample's entry for a record), which lets
    # us know whether we need to have a phased or unphased intermediate region.
    last_call = None
    
    for record in vcf_records:
        # Read through the VCF a record at a time
        
        # Work out the reference position for this variant.
        # What contig is it on?
        reference_contig = record.CHROM
        # Where does it start?
        reference_start = record.POS
        # Where does it end? Depends on the length of the reference allele.
        reference_end = reference_start + len(record.REF)
        
        # Pull out the call we want. See
        # <http://pyvcf.readthedocs.org/en/latest/INTRO.html>. 
        # TODO: what if there's no call for this sample?
        call = record.genotype(sample_name)
        
        if last_call is not None:
            if last_call.site.CHROM == reference_contig:
                # It's on the same contig. Make an AlleleGroup for the region
                # between the last call and this one (the "intermediate"
                # region).
            
                # Where should it start? At the after-the-end position for the
                # last variant.
                intermediate_start = (last_call.site.POS +
                    len(last_call.site.REF))
                
                # Where should it end? At the start position for this variant.
                intermediate_end = reference_start
                
                # And it should obviously be on our contig.
            
                if last_call.phased and call.phased:
                    # It should be phased. We need to make two
                    
                    logging.debug("Intermediate sequence is phased.")
                    
                    for chromosome_number in [0, 1]:
                        # Make an AlleleGroup for each chromosome
                        allele_group = AlleleGroup(reference_contig,
                            intermediate_start, intermediate_end)
                
                        # Set the ploidy to 1 (since it's phased)
                        allele_group["ploidy"] = 1
                        
                        # Tag with the sample's ID
                        allele_group["sample"] = sample_name
                        
                        # Attach it to the correct chromosome.
                        add_allele_group(allele_group, chromosome_number)
                        
                        # Write the AlleleGroup to the file
                        allele_group_writer.append(allele_group)
                    
                else:
                    # It should be unphased.
                    
                    logging.debug("Intermediate sequence is unphased.")
                    
                    # Make one AlleleGroup for both chromosomes
                    allele_group = AlleleGroup(reference_contig,
                        intermediate_start, intermediate_end)
            
                    # Set the ploidy to 2 (since it's unphased)
                    allele_group["ploidy"] = 2
                    
                    # Tag with the sample's ID
                    allele_group["sample"] = sample_name
                    
                    for chromosome_number in [0, 1]:
                        # Attach it to both chromosomes
                        add_allele_group(allele_group, chromosome_number)
                        
                    # Write the AlleleGroup to the file
                    allele_group_writer.append(allele_group)
            else:
                # The last call was on a different contig than this one. Reset
                # our chromosome end tracker.
                logging.debug("New contig: {}".format(reference_contig))
                new_contig()
        
        
        # Now we can do this next variant.
        
        if call.phased:
            # We need to split the genotype sequences with | because they are
            # phased
            separator = "|"
            
            logging.debug("Call at {}:{}-{} is phased".format(reference_contig,
                reference_start, reference_end))
        else:
            # We need to split with / because genotypes are unphased.
            separator = "/"
            
            logging.debug("Call at {}:{}-{} is unphased".format(
                reference_contig, reference_start, reference_end))
        
        for chromosome_number, allele_sequence in enumerate(
            call.gt_bases.split(separator)):
            
            # For each allele base string in the first sample's genotype at this
            # location (along with the number of the chromosome it goes on, 0 or
            # 1)
            
            # Make an AlleleGroup (really just a dict that can generate its own
            # unique IDs.
            allele_group = AlleleGroup(reference_contig, reference_start, 
                reference_end)
                
            # Set the sequence
            allele_group["sequence"] = allele_sequence
            
            # Set the ploidy
            allele_group["ploidy"] = 1
            
            # Tag with the sample's ID
            allele_group["sample"] = sample_name
            
            # Add the AlleleGroup to the appropriate chromosome. If the sample
            # is homozygous we still end up with two AlleleGroups, they just say
            # the same allele.
            add_allele_group(allele_group, chromosome_number)
            
            # Write the AlleleGroup to the file
            allele_group_writer.append(allele_group)
            
        # Remember this call as the last one for the next time we go through
        last_call = call

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object

    # Set up logging
    set_loglevel(options.loglevel)
    
    # Load the Avro schemas
    
    try:
        # This holds the schema for AlleleGroups
        allele_group_schema = avro.schema.parse(
            open(options.allele_group_schema).read())
    except IOError:
        logging.critical("Could not load Avro AlleleGroup schema {}".format(
            options.allele_group_schema))
        logging.critical("Check your --allele_group_schema option")
        sys.exit(1)
        
    try:
        # This holds the schema for Adjacencies
        adjacency_schema = avro.schema.parse(
            open(options.adjacency_schema).read())
    except IOError:
        logging.critical("Could not load Avro Adjacency schema {}".format(
            options.allelegroup_schema))
        logging.critical("Check your --adjacency_schema option")
        sys.exit(1)
        
    # Make Avro-format output file writers. This one is for allele groups.
    allele_group_writer = avro.datafile.DataFileWriter(
        options.allele_group_file, avro.io.DatumWriter(), allele_group_schema)
    # This one is for adjacencies
    adjacency_writer = avro.datafile.DataFileWriter(options.adjacency_file,
        avro.io.DatumWriter(), adjacency_schema)
    
    # Make a VCF reader to read the input VCF
    vcf_reader = vcf.Reader(options.vcf)
    
    # Load all the VCF records into memory. TODO: implement streaming with state
    # for each sample.
    records = list(vcf_reader)
    
    for sample in vcf_reader.samples:
        # Process each sample one at a time
        import_sample(records, sample, allele_group_writer, adjacency_writer)
            
    # Close up the files
    allele_group_writer.close()
    adjacency_writer.close()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
