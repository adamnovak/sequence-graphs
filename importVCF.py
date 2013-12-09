#!/usr/bin/env python2.7
"""
importVCF.py: convert a VCF file to an Avro-based sequence graph format.

This script reads a VCF file containing simple SNP and indels (but no
rearrangements or structural variants) for a single sample. It produces an
unphased sequence graph in two Avro-format files of Avro records: one holding
AlleleGroups, and one holding Adjacencies between AlleleGroups.

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
    parser.add_argument("vcf", type=argparse.FileType("r")
        help="VCF file to read")
    parser.add_argument("allelegroup_file", type=argparse.FileType("w"),
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

def AlleleGroup(dict):
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
    @staticmethod
    def get_id():
        """
        Get the next available unique string ID for AlleleGroups.
        
        """
        
        # Grab the next ID string
        to_return = "ag-{}".format(next_id)
        
        # Advance the next ID
        next_id += 1
        
        # Return the generated ID
        return to_return
    
    def __init__(self, contig, start, end):
        """
        Make a new AlleleGroup representing the region defined by the given
        reference coordinates.
        
        """
        
        # Set up unique IDs
        self["id"] = get_id()
        self["fivePrime"] = "{}-5'".format(self["id"])
        self["threePrime"] = "{}-3'".format(self["id"])
        
        # Save the reference coordinates
        self["contig"] = contig
        self["start"] = start
        self["end"] = end

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
        allelegroup_schema = avro.schema.parse(
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
    
    for record in vcf.Reader(options.vcf):
        # Read through the VCF a record at a time
        
        # Work out the reference position for this variant.
        # What contig is it on?
        reference_contig = record.CHROM
        # Where does it start?
        reference_start = record.POS
        # Where does it end? Depends on the length of the reference allele.
        reference_end = reference_start + len(record.REF)
        
        # This holds the reference-matching AlleleGroup between this variant and
        # the previous one, or None if there isn't any.
        last_constant_segment = None
        
        if (len(last_allele_groups) > 0 and 
            last_alle_groups[0]["contig"] == reference_contig):
            # We had a previous variant on this chromosome. We need adjacencies
            # tying it to an AlleleGroup for the intervening reference DNA, and
            # then later we'll need some adjacencies tying the intervening
            # reference DNA to our new AlleleGroups for this current variant.
            
            # TODO: implement adjacencies
            pass
        
        # Throw out the last pair of AlleleGroups created
        last_allele_groups = []
        
        for allele_sequence in record.samples[0].gt_bases.split("|/"):
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
            
        if last_constant_segment is not None:
            # Add Adjacencies between the previous non-variant AlleleGroup and
            # the two variant AlleleGroups just created.
            
            # TODO: implement this
            pass
            
    
    # Close up the files
    allele_group_writer.close()
    adjacency_writer.close()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
