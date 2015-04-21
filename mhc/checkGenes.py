#!/usr/bin/env python2.7
"""
checkGenes.py: check to see if a MAF aligns things properly according to a set
of gene annotations.

Writes a file of <category>\t<count> TSV lines categorizing all the mappings.

Mapping categories:

"ortholog": when two positions in genes with the same name are mapped together.

"paralog": when two positions in two genes with different names are mapped
together.

"unannotated": when a position that is in a gene and a position that is not are
mapped together.

"background": when two positions, neither of which is in a gene, are mapped
together.

All the pairwise mappings described by the MAF are counted and classified, so a
MAF column with n non-gap characters will produce n choose 2 mappings.

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
    parser.add_argument("maf", type=argparse.FileType("r"), 
        help="MAF file of two genomes to read")
    parser.add_argument("--beds", nargs="+", required=True,
        help=".bed file(s) of genes on the genomes in the MAF")
    parser.add_argument("--classCounts", type=argparse.FileType("w"),
        default=sys.stdout,
        help="output file to save mapping class counts to")
    parser.add_argument("--geneSets", type=argparse.FileType("w"),
        default=sys.stdout,
        help="output file to save sets of genes with mappings in each class to")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def parse_bed(stream):
    """
    Parse a (tab-separated) BED from the given stream. Yield (contig, start,
    end, name, strand) tuples, where strand is 1 or -1. Start is inclusive and
    end is exclusive.
    
    """
    
    # We make a TSV reader to do most of the work. TODO: assumes BED is tab-
    # separated
    reader = tsv.TsvReader(stream)
    
    for line in reader:
        # Break it out into a consistent tuple, and turn the numbers into ints.
        yield (line[0], int(line[1]), int(line[2]), line[3], 
            1 if line[5] == "+" else -1)

def get_columns_from_maf(maf_stream):
    """
    Given a stream of MAF data between two or more contigs, yield individual
    columns from the alignment (for columns that specify alignment).
    
    Each column is a list of tuples of (reference contig, reference base,
    is reverse).
    """
    
    for alignment in AlignIO.parse(maf_stream, "maf"):
        if len(alignment) < 2:
            # We don't have two sequences aligned here.
            continue
            
        for column in xrange(len(alignment[0])):
            # For every column, pull out the bases
            column_bases = alignment[:, column]
            
            # What per-base tuples will we have?
            column_tuples = []
            
            for base, record in itertools.izip(column_bases, alignment):
                # For each record in the alignment
                if base == "-":
                    # Skip stuff that isn't actually aligned here
                    continue
                    
                if record.id == "rootSeq":
                    # This record is fake. TODO: decide how to configure this
                    # sanely.
                    continue
                    
                start_index = record.annotations["start"]
                delta = record.annotations["strand"]
                
                # Put this base in the column, with true for the reverse strand
                # and false for the forward strand as an orientation.
                column_tuples.append((record.id, start_index + column * delta,
                    delta == -1))
                    
            if len(column_tuples) > 1:
                # Send out the column data if there's a real alignment
                yield column_tuples                    

def classify_mappings(columns, genes):
    """
    Given a source of lists of (contig, base, orientation) columns, and a soucre
    of (contig, start, end, gene name, strand number) genes, yield a (class,
    gene, other gene) tuple for each class that each column belongs to.
    
    """
    
    # First we need to build an interval tree for each contig so we can look up
    # what genes overlap each end of each mapping.
    geneTrees = collections.defaultdict(IntervalTree)
    
    for contig, start, end, gene, strand in genes:
        # Put the gene in under the given interval, with strands 1 and -1.
        geneTrees[contig].insert(start, end, (gene, strand))
        
    for column in columns:
        # Look at each column
        
        # Keep track of how many time each gene appears
        gene_counts = collections.Counter()
        
        # Get the gene, orientation pairs for each record in the column
        gene_sets = []
        
        for contig, base, strand in column:
            # For each position in the column, see what genes in what
            # orientations are there.
            genes_found = geneTrees[contig].find(base, base)
            
            # Put them in the opposite direction if we need to based on the
            # strand of the position we mapped to.
            genes_found = {(name, orientation if not strand else -orientation)
                for name, orientation in genes_found}
        
            # Save the set of gene, strand pairs.
            gene_sets.append(genes_found)
                    
        for gene_set in gene_sets:
            # Then total over all the sets of gene, orientation pairs
            for gene in gene_set:
                # Count each gene, orientation pair
                gene_counts[gene] += 1
        
        
        if len(gene_counts) == 0:
            # No genes found. Try the next column
            yield ("background", None, None)
            continue
            
        # Look at the most common gene/orientation pair and its count
        top_gene, top_count = gene_counts.most_common(1)[0]
        
        if top_count > 1:
            # Two or more things have this gene in this orientation, so we have
            # an ortholog mapping.
            yield ("ortholog", top_gene[0], top_gene[0])
            
        if top_count == len(column):
            # This gene in this orientation appears appear in every aligned
            # record.
            continue
            
        # Otherwise there's some record it doesn't appear in. Either a paralog
        # mapping or a mapping to no gene (or a mapping to the same gene in
        # reverse, which we treat as a paralog mapping).
            
        for gene_set in gene_sets:
            # So let's look for a record it isn't in
            
            if top_gene not in gene_set and len(gene_set) > 0:
                # We found a gene set without the most common gene/orientation
                # pair but with some other gene/orientation pair.
                
                # Say there's a paralog mapping between the most common gene and
                # this gene. If their names are the same, their orientations
                # must differ.
                yield ("paralog", top_gene[0], list(gene_set)[0][0])
                break
                
        
        for gene_set in gene_sets:
            # Let's go again but this time find the first empty gene set.
            
            if len(gene_set) == 0:
                # We found a gene set with no genes at all. Say we align
                # this most common gene to a place with no genes.
                yield ("unannotated", top_gene[0], None)
                break
            
def check_genes(maf_filename, genes_filenames):
    """
    Given a MAF filename, and a list of BED filenames holding genes, return a
    Counter from class name to count of mappings in that class, a defaultdict
    from class to a set of gene names involved in each, and a set of gene name
    pairs that align together.
    
    """
    
    # Load all the BEDs, and concatenate them
    genes = itertools.chain.from_iterable([parse_bed(open(bed)) 
        for bed in genes_filenames])
    
    # Get the column iterator, which yields lists of (contig, base, orientation)
    # tuples.
    columns = get_columns_from_maf(maf_filename)
    
    # Classify each column in light of the genes, spitting out classes and gene
    # name pairs.
    classified = classify_mappings(columns, genes)
    
    # This will count up the instances of each class
    class_counts = collections.Counter()
    
    # This will contain sets of genes with any columns in each class
    gene_sets = collections.defaultdict(set)
    
    # This will contain sets of gene pairs aligned together
    gene_pairs = set()
        
    for classification, gene1, gene2 in classified:
        # For each classification
        
        # Count it
        class_counts[classification] += 1
        
        # Record the involved genes
        if gene1 is not None:
            gene_sets[classification].add(gene1)
        if gene2 is not None:
            gene_sets[classification].add(gene2)
            
        if gene1 > gene2:
            # Make sure genes are sorted.
            gene1, gene2 = gene2, gene1
            
        if (gene1, gene2) not in gene_pairs:
            gene_pairs.add((gene1, gene2))
            
        
    # Return the counts and the gene sets and the confusion pairs
    return class_counts, gene_sets, gene_pairs
    
            

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # This will count up the instances of each class, and the genes involved in
    # each class
    class_counts, gene_sets = check_genes(options.maf, options.beds)
    
    for classification, count in class_counts.iteritems():
        # For each class
        
        # Dump a TSV of base counts (over all gene pairs) by classification
        options.classCounts.write("{}\t{}\n".format(classification, count))
        
    for classification, genes in gene_sets.iteritems():
        # For each class
        for gene in genes:
            # For each gene, say it has some mappings in the class.
            options.geneSets.write("{}\t{}\n".format(classification, gene))
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
