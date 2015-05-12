#!/usr/bin/env python2.7
"""
diffGeneSets.py: Compares two files of gene/category multisets.

Each file is a list of <category>\t<gene> lines, which may be repeated if the
same gene ended up in the same category coming from several genomes.

We will work out how many copies of each gene enter or leave each category
between the two files, and report that.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections

import tsv

from Bio import AlignIO, SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# We need yet another comprehensive bio library in Python since BioPython hasn't
# got interval trees.
from bx.intervals import IntervalTree



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
    parser.add_argument("set1", type=argparse.FileType("r"),
        help="TSV file for first gene set")
    parser.add_argument("set2", type=argparse.FileType("r"),
        help="TSV file for second gene set")
    parser.add_argument("--mode", choices=["set", "multiset"],
        default="multiset",
        help="Mode in which to count gene differences")
    parser.add_argument("--outFile", type=argparse.FileType("w"),
        default=sys.stdout,
        help="TSV file for first gene set")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def parse_gene_set(stream):
    """
    Read the stream, parse the gene/category multiset described, and return it.
    
    Returns a defaultdict of counters by category, where each counter has the
    count of genes in that category.
    
    """
    
    # Make the structure we're going to return, counts by category and then
    # gene.
    to_return = collections.defaultdict(collections.Counter)
    
    for parts in tsv.TsvReader(stream):
        # Count this gene in this category
        to_return[parts[0]][parts[1]] += 1
        
    return to_return
    
def diff_gene_sets(set1, set2, mode="multiset"):
    """
    Make a dict by category and then gene of the count differences between these
    two gene sets.
    
    Mode can be "set" or "multiset".
    
    If it is "set", gives +1 for completely new genes in a category, -1 for
    completely removed genes in a category, and 0 for anything else.
    
    If it is "multiset", gives actual count differences in each category.
    
    """
    
    # Make a structure to hold the differences.
    differences = collections.defaultdict(collections.Counter)
    
    for category in set(set1.iterkeys()) | set(set2.iterkeys()):
        # For each category we have anything in
        for gene in (set(set1[category].iterkeys()) | 
            set(set2[category].iterkeys())):
            
            if mode == "multiset":
            
                # For each gene we have anything for, count up the number of
                # copies gained/lost in this category.
                difference = set2[category][gene] - set1[category][gene]
                
            elif mode == "set":
                # Just flag added/removed genes.
                
                # This is our flag. 0 for no change, 1 for added, -1 for
                # removed.
                difference = 0
                
                if set1[category][gene] == 0 and set2[category][gene] > 0:
                    # Gene was added
                    difference = 1
                elif set1[category][gene] > 0 and set2[category][gene] == 0:
                    # Gene was removed
                    difference = -1
            else:
                raise Exception("Mode {} not valid".format(mode))
            
            # Save it
            differences[category][gene] = difference
    
    return differences
                
def describe_difference(differences, stream, numbers=True):
    """
    Given a dict of gene set count differences by category and then gene,
    describe said differences to the given output stream.
    
    """
    
    for category, gene_dict in differences.iteritems():
        # What's the net change for this category
        net = sum(gene_dict.itervalues())
    
        # Each category gets a header.
        stream.write("{}".format(category))
        if numbers:
            # And maybe a net copies added/removed count.
            stream.write(" ({:+})".format(net))
        stream.write(":\n")
        
        for gene, difference in gene_dict.iteritems():
            if difference > 0:
                # Say this category gained copies.
                stream.write("\t+++ {: <16}".format(gene))
                if numbers:
                    # Say how many
                    stream.write("\t{:+}".format(difference))
                stream.write("\n")
            if difference < 0:
                # Say this category lost copies.
                stream.write("\t--- {: <16}".format(gene))
                if numbers:
                    # Say how many
                    stream.write("\t{:+}".format(difference))
                stream.write("\n")
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Parse the two gene sets, diff them, and describe the results.
    describe_difference(diff_gene_sets(parse_gene_set(options.set1),
        parse_gene_set(options.set2), options.mode), options.outFile,
        options.mode == "multiset")
    
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
