#!/usr/bin/env python2.7
"""
checkGenes.py: check to see if a MAF aligns things properly according to a set
of gene annotations.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections

import tsv

from Bio import AlignIO, SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    parser.add_argument("--maf", type=argparse.FileType("r"), required=True, 
        help="MAF file of two genomes to read")
    parser.add_argument("--beds", nargs="+", required=True,
        help=".bed file(s) of genes on the genomes in the MAF")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def parse_bed(stream):
    """
    Parse a (tab-separated) BED from the given stream. Yield (contig, start,
    end, name) tuples.
    
    """
    
    # We make a TSV reader to do most of the work. TODO: assumes BED is tab-
    # separated
    reader = tsv.TsvReader(stream)
    
    for line in reader:
        # Break it out into a consistent tuple, and turn the numbers into ints.
        yield (line[0], int(line[1]), int(line[2]), line[3])

def get_matchings(maf_stream):
    """
    Given a stream of MAF data, yield individual base matchings from the
    alignment.
    
    Each matching is a tuple of (contig, base, other contig, other base, 
    orientation).
    
    """
    
    for alignment in AlignIO.parse(maf_stream, "maf"):
        records = list(alignment)
        
        # This will hold alignment records arranged by the contig they belong
        # to.
        records_by_contig = collections.defaultdict(list)
        
        for record in records:
            # Save it to the right list
            records_by_contig[record.id].append(record)
            
        for (contig1, records1), (contig2, records2) in itertools.combinations(
            records_by_contig.iteritems(), 2):
            # For each pair of contigs, make all the matchings between them...
            
            for record1, record2 in itertools.product(records1, records2):
                # For each pair of records which may induce matchings...
                
                # Where are we on the first sequence?
                index1 = record1.annotations["start"]
                # What direction do we go in?
                delta1 = record1.annotations["strand"]
                
                # Where are we on the second sequence?
                index2 = record2.annotations["start"]
                # What direction do we go in?
                delta2 = record2.annotations["strand"]
                
                for char1, char2 in itertools.izip(record1, record2):
                    if char1 != "-" and char2 != "-":
                        # This is an aligned character. Yield a base matching.
                        yield (contig1, index1, contig2, index2, 
                            (delta1 == delta2))
                            
                    if char1 != "-":
                        # Advance in record 1
                        index1 += delta1
                    if char2 != "-":
                        # Advance in record 2
                        index2 += delta2
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Load all the BEDs
    bedRegions = [list(parse_bed(open(bed))) for bed in options.beds]
    
    
    for matching in get_matchings(options.maf):
        print matching
    
    return
    
    for alignment in AlignIO.parse(options.maf, "maf"):
        records = list(alignment)
        for record in records:
            # Grab the original source sequence name and coordinates
            source_name = record.id
            source_start = record.annotations["start"]
            source_length = record.annotations["size"]
            source_strand = record.annotations["strand"]
            source_end = (source_start + source_length if source_strand == "+1"
                else source_start - source_length)
            
            print record
            print record.annotations
            print("{}:{}-{}".format(source_name, source_start, source_end))
    
    
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
