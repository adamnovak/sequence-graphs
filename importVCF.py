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
    
    # This list holds the two AlleleGroups from the last variant, if any
    last_allele_groups = []
    
    # This holds the last VCF record processed, if any
    last_record = None
    
    for record in vcf.Reader(options.vcf):
        # Read through the VCF a record at a time
        
        # Work out the reference position for this variant.
        # What contig is it on?
        reference_contig = record.CHROM
        # Where does it start?
        reference_start = record.POS
        # Where does it end? Depends on the length of the reference allele.
        reference_end = reference_start + len(record.REF)
        
        # This holds a list of reference-matching AlleleGroups that the next
        # pair of (phased) variant AlleleGroiups should come after, in order. If
        # phasing isn't known, this will just be the same AlleleGroup twice. If
        # there is no reference-matching AlleleGroup that the next pair of
        # variant AlleleGroups should come after, this will be two Nones.
        last_constant_segments = []
        
        if (len(last_allele_groups) > 0 and 
            last_allele_groups[0]["contig"] == reference_contig):
            # We had a previous variant on this chromosome. We need adjacencies
            # tying it to a AlleleGroup(s) for the intervening reference DNA,
            # and then later we'll need some adjacencies tying the intervening
            # reference DNA to our new AlleleGroups for this current variant.
            
            if last_record.samples[0].phased and record.samples[0].phased:
                # We need to preserve phasing
                for previous_variant in last_allele_groups:
                    # Make an AlleleGroup for the intervening reference DNA,
                    # running from the end of the last variant's AlleleGroup
                    # (which we are iterating over) to the start of the variant
                    # we are about to do.
                    constant_segment = AlleleGroup(reference_contig, 
                        previous_variant["end"], reference_start)
                        
                    # Set its ploidy to 1, since it represents only one
                    # chromosome
                    constant_segment["ploidy"] = 1
                        
                    # Write the AlleleGroup to the file
                    allele_group_writer.append(constant_segment)
                 
                    # Make an Adjacency from the end of the variant's
                    # AlleleGroup to the start of the constant one.
                    adjacency = Adjacency(previous_variant["threePrime"],
                        constant_segment["fivePrime"])
                        
                    # Set its ploidy to 1 since we're keeping the phasing.
                    adjacency["ploidy"] = 1
                    
                    # Save it
                    adjacency_writer.append(adjacency)
                    
                    # Keep the constant segment in the list so we can attach to
                    # it later.
                    last_constant_segments.append(constant_segment)
            else:
                # We need to discard phasing
            
                # Make an AlleleGroup for the intervening reference DNA, running
                # from the end of the last variant to the start of this one.
                constant_segment = AlleleGroup(reference_contig, 
                    last_allele_groups[0]["end"], reference_start)
                    
                # Set its ploidy to 2, since it represents both chromosomes
                constant_segment["ploidy"] = 2
                    
                # Write the AlleleGroup to the file
                allele_group_writer.append(constant_segment)
                
                for previous_variant in last_allele_groups:
                    # Make an Adjacency from the end of the variant's
                    # AlleleGroup to the start of the constant one.
                    adjacency = Adjacency(previous_variant["threePrime"],
                        constant_segment["fivePrime"])
                        
                    # Set its ploidy to 1 since we're splitting into phased
                    # chromosomes.
                    adjacency["ploidy"] = 1
                    
                    # Save it
                    adjacency_writer.append(adjacency)
                    
                    # TODO: this code is a bit repetitive. Can we condense it
                    # with the pahse case above?
                
                # Now the last constant segments are this one constant segment
                # we just made, twice.
                last_constant_segments = [constant_segment, constant_segment]
                
        else:
            # We had no previous variant, so there are no intervening constant
            # segments.
            last_constant_segments = [None, None]
            
        # Throw out the last pair of AlleleGroups created
        last_allele_groups = []
        
        if record.samples[0].phased:
            # We need to split the genotype sequences with | because they are
            # phased
            separator = "|"
        else:
            # We need to split with / because genotypes are unphased.
            separator = "/"
        
        for constant_segment, allele_sequence in itertools.izip(
            last_constant_segments, 
            record.samples[0].gt_bases.split(separator)):
            # For each allele base string in the first sample's genotype at this
            # location
            
            # Make an AlleleGroup (really just a dict that can generate its own
            # unique IDs.
            allele_group = AlleleGroup(reference_contig, reference_start, 
                reference_end)
                
            # Set the sequence
            allele_group["sequence"] = allele_sequence
            
            # Keep the AlleleGroup for us to attach our next between-variants
            # segment to.
            last_allele_groups.append(allele_group)
            
            # Write the AlleleGroup to the file
            allele_group_writer.append(allele_group)
            
            if constant_segment is not None:
                # We need to tie this AlleleGroup to the previous constant
                # segment AlleleGroup with an Adjacency.
                
                # Make an adjacency from the previous constant segment to the
                # variant.
                adjacency = Adjacency(constant_segment["threePrime"],
                    allele_group["fivePrime"])
                    
                # Set its ploidy to 1 since at least the variant end is phased.
                adjacency["ploidy"] = 1
                
                # Save it
                adjacency_writer.append(adjacency)
                
        # Remember that we just did this record, so for the next record we can
        # look and see if we need to keep phasing.
        last_record = record
            
    # Close up the files
    allele_group_writer.close()
    adjacency_writer.close()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
