#!/usr/bin/env python2.7
"""
checkGenes.py: check to see if a MAF aligns things properly according to a set
of gene annotations.

Writes a file of <category>\t<count> TSV lines categorizing all the mappings.

TODO: document categories.

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
    parser.add_argument("--maf", type=argparse.FileType("r"), required=True, 
        help="MAF file of two genomes to read")
    parser.add_argument("--beds", nargs="+", required=True,
        help=".bed file(s) of genes on the genomes in the MAF")
    parser.add_argument("--gene2wrongBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file in which to save location of gene2wrong mappings")
    parser.add_argument("--gene2geneBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file in which to save location of gene2gene mappings")
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

def get_mappings(maf_stream):
    """
    Given a stream of MAF data between a single reference contig and a single
    query contig, yield individual base mappings from the alignment.
    
    Each mappings is a tuple of (reference contig, reference base, query contig,
    query base, orientation).
    
    """
    
    for alignment in AlignIO.parse(maf_stream, "maf"):
        records = list(alignment)
        
        if len(records) < 2:
            # We don't have two sequences aligned here.
            continue
        
        # Guess the reference and query
        reference = records[0].id
        query = None
        for record in records:
            if record.id != reference:
                query = record.id
                break
        if query is None:
            # We have multiple records but no query. This might happen if we
            # turned things around on purpose, but for now we should fail.
            raise Exception("Could not find query in {} records".format(
                len(records)))
            
        # This will hold alignment records arranged by the contig they belong
        # to.
        records_by_contig = collections.defaultdict(list)
        
        for record in records:
            # Save it to the right list
            records_by_contig[record.id].append(record)
            
            if record.id not in (reference, query):
                # Complain if we get something not from the reference or query.
                raise Exception("Extra record: {}".format(record.id))
            
        for record1, record2 in itertools.product(records_by_contig[reference],
            records_by_contig[query]):
            # For each pair of records which may induce mappings...
            
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
                    # This is an aligned character. Yield a base mapping.
                    yield (reference, index1, query, index2, 
                        (delta1 == delta2))
                        
                if char1 != "-":
                    # Advance in record 1
                    index1 += delta1
                if char2 != "-":
                    # Advance in record 2
                    index2 += delta2

def classify_mappings(mappings, genes):
    """
    Given a source of (contig, base, other contig, other base, orientation)
    mappings (in reference, query order), and a soucre of (contig, start, end,
    gene name, strand number) genes, yield a (class, gene, mapping) tuple for
    each mapping.
    
    """
    
    # First we need to build an interval tree for each contig so we can look up
    # what genes overlap each end of each mapping.
    geneTrees = collections.defaultdict(IntervalTree)
    
    for contig, start, end, gene, strand in genes:
        # Put the gene in under the given interval
        geneTrees[contig].insert(start, end, (gene, strand))
        
    for mapping in mappings:
        # Look at each mapping
        
        # Unpack the mapping
        contig1, base1, contig2, base2, orientation = mapping
        
        # Pull the name, strand tuples for both ends
        genes1 = set(geneTrees[contig1].find(base1, base1))
        genes2 = set(geneTrees[contig2].find(base2, base2))
        
        if orientation == False:
            # This is a backwads mapping, so flip orientations on one of the
            # sets
            genes2 = {(name, -strand) for name, strand in genes2}
        
        # Reference was first, so we want to see if the genes we were supposed
        # to have in the query were recapitulated in the right orientation in
        # the reference.
        
        # What gene names from the query should we use to describe this mapping,
        # if any?
        
        # TODO: How many times do genes with different names actually overlap? I
        # bet never.
        
        # Make a set of all the reference gene names.
        gene_names1 = {name for name, base in genes1}
        # Make a set of all the query gene names.
        gene_names2 = {name for name, base in genes2}
        
        if(len(genes1 & genes2)) > 0:
            # At least one shared gene exists in the right orientation to
            # explain this mapping.
            
            # Grab a plausible gene name.
            gene_name = next(iter(genes1 & genes2))[0]
            
            yield "gene2gene", gene_name, mapping
        elif len(genes1) == 0 and len(genes2) > 0:
            # We mapped a gene to somewhere where there is no gene
            
            # Grab a plausible gene name.
            gene_name = next(iter(gene_names2))
            
            yield "gene2non", gene_name, mapping
        
        elif len(genes1) > 0 and len(genes2) > 0:
            # We mapped a gene to a place where there are genes, but not the
            # right gene in the right orientation
            
            # Grab a plausible gene name.
            gene_name = next(iter(gene_names2))
            
            yield "gene2wrong", gene_name, mapping
        elif len(genes1) > 0 and len(genes2) == 0:
            # We mapped a non-gene to somewhere where there are genes
            yield "non2gene", None, mapping
        else:
            # We mapped a non-gene to somewhere where there are no genes.
            yield "non2non", None, mapping

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Make a dict from class to BED output stream, or None
    classBeds = collections.defaultdict(lambda: None)
    
    # Populate it with the streams from the options
    classBeds["gene2wrong"] = options.gene2wrongBed
    classBeds["gene2gene"] = options.gene2geneBed
    
    
    # Load all the BEDs, and concatenate them
    genes = itertools.chain.from_iterable([parse_bed(open(bed)) 
        for bed in options.beds])
    
    # Get all the mappings
    mappings = get_mappings(options.maf)
    
    # Classify each mapping in light of the genes
    classified = classify_mappings(mappings, genes)
    
    # This will count up the instances of each class
    classCounts = collections.Counter()
    
    # This will hold sets of genes that have any bases in each category.
    geneSets = collections.defaultdict(set)
    
    for classification, gene, mapping in classified:
        # For each classified mapping
        
        # Record it in its class
        classCounts[classification] += 1
        
        if classBeds[classification] is not None:
            # We want to write a BED of this class.
            
            # Unpack the mapping
            contig1, base1, contig2, base2, orientation = mapping
            
            # Dump a BED record in query coordinates
            classBeds[classification].write("{}\t{}\t{}\t{}".format(contig2, 
                base2, base2 + 1, classification))
                
        if gene is not None:
            # Note that this gene has some bases in this mapping class
            geneSets[classification].add(gene)
        
        
    
    for classification, count in classCounts.iteritems():
        # For each class
        
        # Dump a TSV of bases by classification
        options.classCounts.write("{}\t{}\n".format(classification, count))
        
    for classification, genes in geneSets.iteritems():
        # For each set of query genes with mappings in a class
        
        for gene in genes:
            # Save each observed gene, one per line
            options.geneSets.write("{}\t{}\n".format(classification, gene))
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
