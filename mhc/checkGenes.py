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

def get_mappings_from_maf(maf_stream):
    """
    Given a stream of MAF data between two or more contigs, yield individual
    base mappings from the alignment.
    
    Each mappings is a tuple of (reference contig, reference base, query contig,
    query base, is backwards).
    """
    
    for alignment in AlignIO.parse(maf_stream, "maf"):
        records = list(alignment)
        
        if len(records) < 2:
            # We don't have two sequences aligned here.
            continue
        
        for record1, record2 in itertools.combinations(records, 2):
            # For each pair of records which may induce mappings...
            
            if record1.id == record2.id:
                # Skip self alignments
                continue
            
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
                    # This is an aligned character. Yield a base mapping. We
                    # only yield each mapping in one direction.
                    yield (record1.id, index1, record2.id, index2,
                        (delta1 != delta2))
                        
                if char1 != "-":
                    # Advance in record 1
                    index1 += delta1
                if char2 != "-":
                    # Advance in record 2
                    index2 += delta2

def classify_mappings(mappings, genes):
    """
    Given a source of (contig, base, other contig, other base, orientation)
    mappings, and a soucre of (contig, start, end, gene name, strand number)
    genes, yield a (class, gene, other gene, mapping) tuple for each mapping.
    
    """
    
    # First we need to build an interval tree for each contig so we can look up
    # what genes overlap each end of each mapping.
    geneTrees = collections.defaultdict(IntervalTree)
    
    for contig, start, end, gene, strand in genes:
        # Put the gene in under the given interval
        geneTrees[contig].insert(start, end, (gene, strand))
        
    for mapping in mappings:
        # Look at each mapping
        
        if len(mapping) == 5:
        
            # Unpack the mapping
            contig1, base1, contig2, base2, orientation = mapping
            
            # Pull the name, strand tuples for both ends
            genes1 = set(geneTrees[contig1].find(base1, base1))
            genes2 = set(geneTrees[contig2].find(base2, base2))
            
            if orientation:
                # This is a backwads mapping, so flip orientations on one of the
                # sets
                genes2 = {(name, -strand) for name, strand in genes2}
            
            # What gene names should we use to describe this mapping, if any?
            
            # TODO: How many times do genes with different names actually
            # overlap? I bet never.
            
            # Make a set of all the contig1 gene names.
            gene_names1 = {name for name, _ in genes1}
            # Make a set of all the contig2 gene names.
            gene_names2 = {name for name, _ in genes2}
            
            if(len(genes1 & genes2)) > 0:
                # At least one shared gene exists in the right orientation to
                # explain this mapping.
                
                # Grab a plausible gene name.
                gene_name = next(iter(genes1 & genes2))[0]
                
                yield "ortholog", gene_name, gene_name, mapping
            elif len(genes1) == 0 and len(genes2) > 0:
                # We mapped a gene to somewhere where there is no gene
                
                # Grab a plausible gene name.
                gene_name = next(iter(gene_names2))
                
                yield "unannotated", gene_name, None, mapping
            
            elif len(genes1) > 0 and len(genes2) > 0:
                # We mapped a gene to a place where there are genes, but not the
                # right gene in the right orientation
                
                # Grab a plausible gene name from.
                gene_from = next(iter(gene_names2))
                
                # And where it could be mapped to
                gene_to = next(iter(gene_names1))
                
                yield "paralog", gene_from, gene_to, mapping
            elif len(genes1) > 0 and len(genes2) == 0:
                # We mapped a non-gene to somewhere where there are genes
                
                # What's a plausible gene name we could be going to?
                gene_to = next(iter(gene_names1))
                
                yield "unannotated", None, gene_to, mapping
            else:
                # We mapped a non-gene to somewhere where there are no genes.
                yield "background", None, None, mapping
                
            # OK that's all the cases for mapped bases.
            
def check_genes(maf_filename, genes_filenames):
    """
    Given a MAF filename, and a list of BED filenames holding genes, return a
    Counter from class name to count of mappings in that class, and a
    defaultdict from class to a set of gene names involved in each.
    
    """
    
    # Load all the BEDs, and concatenate them
    genes = itertools.chain.from_iterable([parse_bed(open(bed)) 
        for bed in genes_filenames])
    
    # Get all the mappings. Mappings are tuples in the form (reference, 
    # offset, query, offset, query read name, is backwards). In a maf the
    # query read name is just the query name.
    mappings = get_mappings_from_maf(maf_filename)
    
    # Classify each mapping in light of the genes, tagging it with a
    # classiication, and the genes it is from and to (which may be None).
    classified = classify_mappings(mappings, genes)
    
    # This will count up the instances of each class
    class_counts = collections.Counter()
    
    gene_sets = collections.defaultdict(set)
        
    for classification, gene1, gene2, mapping in classified:
        # For each classified mapping
        
        # We have two positions and an orientatin in the mapping.
        contig1, base1, contig2, base2, orientation = mapping
        
        # Record it in its class for its gene pair (which may be Nones)
        class_counts[classification] += 1
        
        # Record the involved genes
        if gene1 is not None:
            gene_sets[classification].add(gene1)
        if gene2 is not None:
            gene_sets[classification].add(gene2)
            
        
    # Return the counts and the gene sets
    return class_counts, gene_sets
    
            

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
