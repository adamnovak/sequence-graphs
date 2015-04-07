#!/usr/bin/env python2.7
"""
fetchRegion.py: Fetch the sequence data, GRC alignments, and gene sets for a GRC
region (like "LRC_KIR" or "MHC") by name.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections, urllib2

import tsv

from Bio import AlignIO, SeqIO, Align, Entrez
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
    parser.add_argument("region",
        help="name of the region to download, and the output directory")
    parser.add_argument("--assembly_url", 
        default=("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
            "GCA_000001405.17_GRCh38.p2/"
            "GCA_000001405.17_GRCh38.p2_assembly_structure"),
        help="URL for the assembly, containing genomic_region_definitions.txt")
    parser.add_argument("--email", default="anovak@soe.ucsc.edu",
        help="E-mail address to report to Ensembl")

    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def getRegionInfo(region_name, assembly_root):
    """
    Go download the genomic_region_definitions.txt from the specified assembly,
    and return the (contig, start, end) of the named region.
    
    """
    
    # Open the region definitions
    for parts in tsv.TsvReader(urllib2.urlopen(assembly_root +
        "/genomic_region_definitions.txt")):
        
        # For every region in the list
        
        if parts[0] == region_name:
            # We found it. Parse out contig, start, end.
            # Contig is like "CM000663.2" and not all that useful...
            return (parts[1], int(parts[2]), int(parts[3]))
            
    # If we get here, there's no region by that name.
    raise RuntimeError("No region named " + region_name)
    
def getRegionSequences(region_name, assembly_root):
    """
    Given the name of a region and the root URL for the assembly, yield all the
    alt locus sequence names (Genbank IDs) that are in the region.
    
    """
    
    # Open the alt locus placement file
    for parts in tsv.TsvReader(urllib2.urlopen(assembly_root +
        "/all_alt_scaffold_placement.txt")):
        # For every alt locus...
        
        if parts[7] == region_name:
            # We found one in the correct region (which happens to be column 7)
            
            # Give its sequence ID/accession with version (column 3)
            yield parts[3]
    
    
    
def getGINumber(grc_id):
    """
    Given a GRC-style (genbank) ID with version, like "CM000663.2" or
    "GL383549.1" or "KI270832.1", get the GI number associated with that ID from
    Entrez.
    
    """
    
    # First just search the ID as a search term.
    search_results = Entrez.read(Entrez.esearch("nuccore", term="GL000209.2"))
    # Grab the handle thingy for the first search result
    first_handle = Entrez.read(Entrez.epost("nuccore",
        id=search_results["IdList"][0]))
    # Actually download that record
    record = Entrez.read(Entrez.esummary(db="nuccore",
        webenv=first_handle["WebEnv"], query_key=first_handle["QueryKey"]))
        
    # Return the GI number. TODO: should this be the ID instead because of how
    # we use it next? Are they ever different?
    return record["Gi"]
    
def getSequence(gi_id, start=None, end=None):
    """
    Get a sequence by numerical GI number, optionally with start and end
    parameters (in 0-based (?) coordinates from the left). If start is
    specified, end must also be specified.
    
    """
    
    if start is None:
        # Go fetch the whole record
        fetch_handle = Entrez.efetch(db="nucleotide", id=gi_id, rettype="fasta")
    else:
        # Just fetch part of it
        fetch_handle = Entrez.efetch(db="nucleotide", id=gi_id, rettype="fasta",
            seq_start=start, seq_end=end)
    
    # Load up the FASTA record
    record = SeqIO.read(fetch_handle, "fasta")
    
    # Change the record FASTA ID to just GIwhatever
    record.id = "GI{}".format(gi_id)
    
    # Return the fixed-up record
    return record   
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Set Ensembl e-mail
    Entrez.email = options.email
    
    # Go get the region of the reference we're talking about
    ref_acc, ref_start, ref_end = getRegionInfo(options.region,
        options.assembly_url)
        
    # Make our output directory
    os.makedirs(options.region)
    
    # Get the reference's GI
    ref_gi = getGINumber(ref_acc)
    
    print("Reference for {} is GI{}:{}-{}".format(options.region, ref_gi,
        ref_start, ref_end))
    
    # Grab the reference sequence
    ref_seq = getSequence(ref_gi, ref_start, ref_end)
    
    # Change it to be just called "ref"
    ref_seq.id = "ref"
    
    # Write it to <region>/ref.fa
    SeqIO.write([ref_seq], open("{}/ref.fa".format(options.region), "w"),
        "fasta")
    
    for alt_acc in getRegionSequences(options.region, options.assembly_url):
        # For every alt in the region
        
        # Get its GI number
        alt_gi = getGINumber(alt_acc)
        
        print("Downloading alt GI{}".format(alt_gi))
        
        # Grab the sequence data
        alt_seq = getSequence(alt_gi)
        
        # Write it to <region>/GI<number>.fa
        SeqIO.write([alt_seq], open("{}/GI{}.fa".format(options.region, alt_gi),
            "w"), "fasta")
        
            
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
