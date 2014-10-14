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
    parser.add_argument("--maf", type=argparse.FileType("r"), 
        help="MAF file of two genomes to read")
    parser.add_argument("--tsv", type=argparse.FileType("r"),
        help="TSV file of mappings of reads to load")
    parser.add_argument("--beds", nargs="+", required=True,
        help=".bed file(s) of genes on the genomes in the MAF")
    parser.add_argument("--gene2wrongBed", type=argparse.FileType("w"),
        default=None,
        help=".bed for mappings to the wrong gene")
    parser.add_argument("--gene2geneBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file for mappings to the right gene")
    parser.add_argument("--non2geneBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file for mappings of non-genes to genes")
    parser.add_argument("--gene2nonBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file for mappings of genes to non-genes")
    parser.add_argument("--non2nonBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file for mappings of non-genes to non-genes")
    parser.add_argument("--gene2unmappedBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file for unmapped positions in genes")
    parser.add_argument("--non2unmappedBed", type=argparse.FileType("w"),
        default=None,
        help=".bed file for unmapped positions in non-genes")
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
    Given a stream of MAF data between a single reference contig and a single
    query contig, yield individual base mappings from the alignment.
    
    Each mappings is a tuple of (reference contig, reference base, query contig,
    query base, is backwards).
    
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
                        (delta1 != delta2))
                        
                if char1 != "-":
                    # Advance in record 1
                    index1 += delta1
                if char2 != "-":
                    # Advance in record 2
                    index2 += delta2
                    
def get_mappings_from_tsv(tsv_stream):
    """
    Given a TSV stream of mappings like
    <reference>\t<position>\t<query>\t<position>\t<isBackwards> (or
    <query>\t<position> for unmapped bases), where the query name is of the form
    <actual query contig>:<start>-<end>, parse out and yield mappings between
    the actual full query sequence and the reference.
    
    """
    
    for parts in tsv.TsvReader(tsv_stream):
        if len(parts) == 5:
            # This base is mapped
    
            # Parse out everything
            # What is the reference we're aligned to?
            reference_name = parts[0]
            # Where on it are we?
            reference_position = int(parts[1])
            # What is the query we're aligning?
            query = parts[2]
            # Where on it are we?
            query_position = int(parts[3])
            # Are we aligned to the reverse strand?
            is_backwards = bool(int(parts[4]))
            
            # Parse the query name more
            query_parts = query.split(":")
            # What is the name of the orifginal un-split query sequence?
            query_name = query_parts[0]
            query_parts = query_parts[1].split("-")
            # Where on it did this sequence start?
            query_offset = int(query_parts[0])
            
            # Assemble a mapping of the original query sequence and yield it.
            yield (reference_name, reference_position, query_name, 
                query_position + query_offset, is_backwards)
        elif len(parts) == 2:
            # This base is unmapped
            # Parse out everything
            # What query holds the base?
            query = parts[0]
            # Where on it are we?
            query_position = int(parts[1])
            
            # Parse the query name more
            query_parts = query.split(":")
            # What is the name of the orifginal un-split query sequence?
            query_name = query_parts[0]
            query_parts = query_parts[1].split("-")
            # Where on it did this sequence start?
            query_offset = int(query_parts[0])
            
            yield (query_name, query_position + query_offset)
            

def classify_mappings(mappings, genes):
    """
    Given a source of (contig, base, other contig, other base, orientation) or
    (other contig, other base) mappings (in reference, query order), and a
    soucre of (contig, start, end, gene name, strand number) genes, yield a
    (class, fromGene, toGene, mapping) tuple for each mapping.
    
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
            
            # Reference was first, so we want to see if the genes we were
            # supposed to have in the query were recapitulated in the right
            # orientation in the reference.
            
            # What gene names from the query should we use to describe this
            # mapping, if any?
            
            # TODO: How many times do genes with different names actually
            # overlap? I bet never.
            
            # Make a set of all the reference gene names.
            gene_names1 = {name for name, _ in genes1}
            # Make a set of all the query gene names.
            gene_names2 = {name for name, _ in genes2}
            
            if(len(genes1 & genes2)) > 0:
                # At least one shared gene exists in the right orientation to
                # explain this mapping.
                
                # Grab a plausible gene name.
                gene_name = next(iter(genes1 & genes2))[0]
                
                yield "gene2gene", gene_name, gene_name, mapping
            elif len(genes1) == 0 and len(genes2) > 0:
                # We mapped a gene to somewhere where there is no gene
                
                # Grab a plausible gene name.
                gene_name = next(iter(gene_names2))
                
                yield "gene2non", gene_name, None, mapping
            
            elif len(genes1) > 0 and len(genes2) > 0:
                # We mapped a gene to a place where there are genes, but not the
                # right gene in the right orientation
                
                # Grab a plausible gene name from.
                gene_from = next(iter(gene_names2))
                
                # And where it could be mapped to
                gene_to = next(iter(gene_names1))
                
                yield "gene2wrong", gene_from, gene_to, mapping
            elif len(genes1) > 0 and len(genes2) == 0:
                # We mapped a non-gene to somewhere where there are genes
                
                # What's a plausible gene name we could be going to?
                gene_to = next(iter(gene_names1))
                
                yield "non2gene", None, gene_to, mapping
            else:
                # We mapped a non-gene to somewhere where there are no genes.
                yield "non2non", None, None, mapping
                
            # OK that's all the cases for mapped bases.
            
        elif len(mapping) == 2:
            # Unpack the mapping (numbered 2 for consistency with above)
            contig2, base2 = mapping
            
            # Pull a set of name, strand pairs for this position
            genes2 = set(geneTrees[contig2].find(base2, base2))
            
            # Pull out just the names
            gene_names = [name for name, _ in genes2]
            
            # Grab a plausible gene name from.
            gene_from = gene_names[0] if len(gene_names) > 0 else None
            
            if gene_from is None:
                # No gene was here.
                yield "non2unmapped", None, None, mapping
            else:
                # There was a gene here. Report it.
                yield "gene2unmapped", gene_name, None, mapping
            

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Make a dict from class to BED output stream, or None
    class_beds = collections.defaultdict(lambda: None)
    
    # Populate it with TSV writers for the streams from the options
    # TODO: this is awkward and needs a map-over-option like Scala.
    class_beds["gene2wrong"] = tsv.TsvWriter(options.gene2wrongBed) \
        if options.gene2wrongBed is not None else None
    class_beds["gene2gene"] = tsv.TsvWriter(options.gene2geneBed) \
        if options.gene2geneBed is not None else None
    class_beds["gene2non"] = tsv.TsvWriter(options.gene2nonBed) \
        if options.gene2nonBed is not None else None
    class_beds["gene2unmapped"] = tsv.TsvWriter(options.gene2unmappedBed) \
        if options.gene2unmappedBed is not None else None
    class_beds["non2gene"] = tsv.TsvWriter(options.non2geneBed) \
        if options.non2geneBed is not None else None
    class_beds["non2non"] = tsv.TsvWriter(options.non2nonBed) \
        if options.non2nonBed is not None else None
    class_beds["non2unmapped"] = tsv.TsvWriter(options.non2unmappedBed) \
        if options.non2unmappedBed is not None else None
    
    # Make a writer for gene set output
    gene_set_writer = tsv.TsvWriter(options.geneSets)
    
    # Load all the BEDs, and concatenate them
    genes = itertools.chain.from_iterable([parse_bed(open(bed)) 
        for bed in options.beds])
    
    if options.maf is not None:
        # Get all the mappings
        mappings = get_mappings_from_maf(options.maf)
    elif options.tsv is not None:
        # Get all the mappings by re-assembling reads mapped in a TSV
        mappings = get_mappings_from_tsv(options.tsv)
    else:
        # TODO make argparse check this
        raise Exception("No mappings provided")
        
    
    # Classify each mapping in light of the genes
    classified = classify_mappings(mappings, genes)
    
    # This will count up the instances of each class
    class_counts = collections.Counter()
    
    # This maps from classification, source gene, and gene mapped to to genome set that has that mapping.
    gene_sets = collections.defaultdict(lambda: collections.defaultdict(
        lambda: collections.defaultdict(set)))
    
    for classification, gene_from, gene_to, mapping in classified:
        # For each classified mapping
        
        # Unpack the mapping
        if len(mapping) == 5:
            # We have query and reference positions
            _, _, contig2, base2, _ = mapping
        elif len(mapping) == 2:
            # We only have the wuery position
            contig2, base2 = mapping
        
        # Record it in its class
        class_counts[classification] += 1
        
        if class_beds[classification] is not None:
            # We want to write a BED of this class.
            
            # Dump a BED record in query coordinates, named after the gene we
            # mapped to.
            class_beds[classification].line(contig2, base2, base2 + 1,
                gene_to)
                
        # Note that this mapped genome has some bases mapping from this one gene
        # to this other gene (either of which may be None).
        gene_sets[classification][gene_from][gene_to].add(contig2)
        
        
    
    for classification, count in class_counts.iteritems():
        # For each class
        
        # Dump a TSV of bases by classification
        options.classCounts.write("{}\t{}\n".format(classification, count))
        
    for classification, gene_mappings in gene_sets.iteritems():
        # For each set of query genes with mappings in a class
        
        for gene, other_genes in gene_mappings.iteritems():
            # For each gene, what other genes was it observed mapping to?
            
            for other_gene, genome_set in other_genes.iteritems():
                # For each mapping of this gene, what genomes did it appear in?
                gene_set_writer.line(classification, gene, other_gene,
                    ",".join(sorted(genome_set)))
                

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
