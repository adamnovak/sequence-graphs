#!/usr/bin/env python2.7
"""
fetchRegion.py: Fetch the sequence data, GRC alignments, and gene sets for a GRC
region (like "LRC_KIR" or "MHC") by name.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections, urllib2, shutil, subprocess, glob

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
        help="E-mail address to report to Entrez")

    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def url_open_tsv(url):
    """
    Open a TSV URL and loop through the lines as lists.
    
    """
    
    try:
        reader = tsv.TsvReader(urllib2.urlopen(url))
    except urllib2.URLError as err:
        print("Could not open " + url)
        raise err
        
    return reader

def get_region_info(region_name, assembly_root):
    """
    Go download the genomic_region_definitions.txt from the specified assembly,
    and return the (contig, start, end) of the named region.
    
    """
    
    # Open the region definitions
    for parts in url_open_tsv(assembly_root + "/genomic_regions_definitions.txt"):
        # For every region in the list
        
        if parts[0] == region_name:
            # We found it. Parse out contig, start, end.
            # Contig is like "CM000663.2" and not all that useful...
            return (parts[1], int(parts[2]), int(parts[3]))
            
    # If we get here, there's no region by that name.
    raise RuntimeError("No region named " + region_name)
    
def get_region_sequences(region_name, assembly_root):
    """
    Given the name of a region and the root URL for the assembly, yield all the
    alt locus sequence names (Genbank IDs) and assembly unit names
    (ALT_REF_LOCI_###) that are in the region.
    
    """
    
    # Open the alt locus placement file
    for parts in url_open_tsv(assembly_root + "/all_alt_scaffold_placement.txt"):
        # For every alt locus...
        
        if parts[7] == region_name:
            # We found one in the correct region (which happens to be column 7)
            
            # Give its sequence ID/accession with version (column 3) and the
            # assembly unit name
            yield parts[3], parts[0]

def get_record_by_grc(grc_id):
    """
    Given a GRC ID, return the Entrez nucleotide DocSum record.
    
    """

    # First just search the ID as a search term.
    search_results = Entrez.read(Entrez.esearch("nucleotide", term=grc_id))
    
    if len(search_results["IdList"]) > 1:
        # We should only get one result. If we have many, we might be looking at
        # the wrong one.
        print(search_results)
        raise RuntimeError("Too many results!")
    
    # Grab the handle thingy for the first search result
    first_handle = Entrez.read(Entrez.epost("nucleotide",
        id=search_results["IdList"][0]))
        
    # Actually download that record
    record = Entrez.read(Entrez.esummary(db="nucleotide",
        webenv=first_handle["WebEnv"], query_key=first_handle["QueryKey"]))[0]
        
    # Return it
    return record
    
def get_ucsc_name(grc_id, alt_parent_grc_id=None):
    """
    Given a GRC-style (genbank) ID with version, like "CM000663.2" or
    "GL383549.1" or "KI270832.1", get the UCSC name for that sequence, like
    "chr6_GL000252v2_alt".
    
    If the sequence is an alt, the GRC id of its parent chromosome must be
    specified.
    
    """
    
    if alt_parent_grc_id is None:
        # Simple case; it's a primary chromosome.
        
        # Fetch the record
        record = get_record_by_grc(grc_id)
        
        # Parse out all the "extra" fields
        extra_parts = record["Extra"].split("|")
        
        # Find the "gnl" key
        gnl_index = extra_parts.index("gnl")
        
        # The chromosome number/letter is two fields later.
        chromosome_character = extra_parts[gnl_index + 2]
        
        # Make it chrThat
        ucsc_name = "chr{}".format(chromosome_character)
        
    else:
        # We do have a parent. Get its UCSC name.
        parent_name = get_ucsc_name(alt_parent_grc_id)
        
        # Convert from .2 or whatever to v2 or whatever
        name_middle = grc_id.replace(".", "v")
        
        # Put them in the name pattern template to generate the name
        ucsc_name = "{}_{}_alt".format(parent_name, name_middle)
        
    # Report the result
    print("{} is {} at UCSC".format(grc_id, ucsc_name))
    return ucsc_name
        
        
    
def get_gi_number(grc_id):
    """
    Given a GRC-style (genbank) ID with version, like "CM000663.2" or
    "GL383549.1" or "KI270832.1", get the GI number associated with that ID from
    Entrez.
    
    """
    
    # Go fetch the record
    record = get_record_by_grc(grc_id)
        
    print("{} = {}".format(grc_id, record["Gi"]))
        
    # Return the GI number. TODO: should this be the ID instead because of how
    # we use it next? Are they ever different?
    return record["Gi"]
    
def get_length(gi_id):
    """
    Get the length of a sequence given its numerical GI number.
    
    """
    
    # Grab the handle thingy for the record with this ID
    handle = Entrez.read(Entrez.epost("nucleotide", id=str(gi_id)))
        
    # Actually download that record
    record = Entrez.read(Entrez.esummary(db="nucleotide",
        webenv=handle["WebEnv"], query_key=handle["QueryKey"]))[0]
        
    # Return the length of the sequence
    return record["Length"]
    
    
    
def get_sequence(gi_id, start=None, end=None):
    """
    Get a sequence by numerical GI number, optionally with start and end
    parameters (in 1-based coordinates from the left). If start is
    specified, end must also be specified.
    
    """
    
    if start is None:
        # Go fetch the whole record. We need to make the ID a str or the API
        # client freaks out.
        fetch_handle = Entrez.efetch(db="nucleotide", id=str(gi_id),
            rettype="fasta")
    else:
        # Just fetch part of it
        fetch_handle = Entrez.efetch(db="nucleotide", id=str(gi_id),
            rettype="fasta", seq_start=start, seq_end=end)
    
    # Load up the FASTA record
    record = SeqIO.read(fetch_handle, "fasta")
    
    # Change the record FASTA ID to just GIwhatever
    record.id = "GI{}".format(gi_id)
    
    # Return the fixed-up record
    return record
    
def download_gff3(ref_acc, alt_acc, alt_unit, assembly_root, out_filename):
    """
    Download the GFF3 alignment between the given reference accession and the
    given alt accession (in the given assembly unit), from the given assembly
    root URL, and save it to the given output filename.
    
    """
    
    # Figure out what we want to download
    gff3_url = "{}/{}/alt_scaffolds/alignments/{}_{}.gff".format(assembly_root,
        alt_unit, alt_acc, ref_acc)

    # Open the URL to read
    in_stream = urllib2.urlopen(gff3_url)
    with open(out_filename, "w") as out_stream:
        # Copy everything to the output file as in
        # <http://stackoverflow.com/a/5397438/402891>
        shutil.copyfileobj(in_stream, out_stream)
        
def get_genes(grc_id, out_name, start=1, end=None, alt_parent_grc_id=None,
    db="hg38"):
    """
    Given a GRC ID (like "CM000663.2"), the name of the contig on which to
    report the genes, optional start and end coordinates (1-based) and the GRC
    ID of the parent chromosome if it is an alt, yield BED lines for all the
    genes in the specified region.
    
    If start is specified, coordinates will be given relative to that position.
    
    Assumes "hgsql" is installed and configured and available on the PATH.
    
    Uses the hg38 database unless told otherwise.
    
    All inputs must be trusted and not permitted to contain SQL injection.
    
    """
    
    # Convert to 0-based not-end-inclusive coordinates.
    start -= 1
    
    # Get the name to look up in the database.
    query_contig = get_ucsc_name(grc_id, alt_parent_grc_id)
    
    # Spec out the query. TODO: Can I not say the database name constantly?
    query_parts = ["SELECT \"", out_name, "\", ", db, ".knownGene.txStart - ", 
        start, ", ", db, ".knownGene.txEnd - ", start, ", ", db, 
        ".kgXref.geneSymbol, 0, ", db, ".knownGene.strand FROM ", db, 
        ".knownGene LEFT OUTER JOIN ", db, ".kgXref ON ", db, 
        ".knownGene.name = ", db, ".kgXref.kgID WHERE ", db, 
        ".knownGene.txStart != ", db, ".knownGene.cdsStart AND ", db, 
        ".knownGene.chrom = \"", query_contig, "\" AND ", db, 
        ".knownGene.txStart >= ", start]
    
    if end is not None:
        # Require the end criterion to be met too.
        query_parts += [" AND ", db, ".knownGene.txEnd < ", end]
    
    # Finish off the query.
    query_parts.append(";")
    
    # Put together the whole query.
    query = "".join([str(part) for part in query_parts])
    
    # Build the hgsql command line
    args = ["hgsql", "-e",  query]
    
    # Open the process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
        
    for line in itertools.islice(process.stdout, 1, None):
        # For all lines except the first, yield them because they are BED lines.
        yield line
            
    if process.wait() != 0:
        raise RuntimeError("hgsql")
            
    # We are done with this process.
    process.stdout.close()
    
def open_gene_bed(region, sequence_id):
    """
    Given the region name and the sequence ID ("ref" or "GI<whatever>") for a
    sequence, give back an output file object to which a BED of the genes in
    that sequence may be written.
    
    """
    
    # Each bed goes in a folder named after its sequence, so hal2assemblyHub can
    # use it.
    bed_dir = "{}/genes/{}".format(region, sequence_id)
    
    if not os.path.exists(bed_dir):
        # Make sure we have a place to put the genes
        os.makedirs(bed_dir)
    
    # Open a bed file in there for writing and return it.
    return open(bed_dir + "/genes.bed", "w")
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Set Entrez e-mail
    Entrez.email = options.email
    
    # Go get the region of the reference we're talking about. Starts and ends
    # are 1-based.
    ref_acc, ref_start, ref_end = get_region_info(options.region,
        options.assembly_url)
        
    # Make our output directory
    if not os.path.exists(options.region):
        os.makedirs(options.region)
        
    # We're going to write a chrom.sizes file with accessions (not the GI
    # numbers) for the gff3->psl conversion step
    acc_chrom_sizes = tsv.TsvWriter(open(options.region + "/acc.chrom.sizes",
        "w"))
        
    # Get the reference's GI
    ref_gi = get_gi_number(ref_acc)
    
    print("Reference for {} is GI{}:{}-{} 1-based".format(options.region,
        ref_gi, ref_start, ref_end))
    
    # Grab the reference sequence
    ref_seq = get_sequence(ref_gi, ref_start, ref_end)
    
    print("Got {}bp for a {}bp reference".format(len(ref_seq),
        ref_end - ref_start + 1))
        
        
    if len(ref_seq) > ref_end - ref_start + 1:
        # Clip it down if it's too long. Assuming we have the correct sort of
        # coordinates, and that we got served the data starting at the correct
        # offset.
        ref_seq = ref_seq[0:ref_end - ref_start + 1]
    elif len(ref_seq) < ref_end - ref_start:
        raise RuntimeError("Didn't get enough sequence from the API!")
    
    # Change it to be just called "ref"
    ref_seq.id = "ref"
    
    # Write it to <region>/ref.fa
    SeqIO.write([ref_seq], open("{}/ref.fa".format(options.region), "w"),
        "fasta")
        
    # Write a chromosome size entry for the reference by its accession
    acc_chrom_sizes.line(ref_acc, get_length(ref_gi))
    
    print("Writing genes for ref")
    
    # Make a BED to put reference genes in
    ref_bed = open_gene_bed(options.region, "ref") 
    
    for line in get_genes(ref_acc, "ref", ref_start, ref_end):
        # Write all the BED lines for the appropriate region of the reference to
        # that file.
        ref_bed.write(line)
        
    ref_bed.close()
    
    for alt_acc, alt_unit in get_region_sequences(options.region,
        options.assembly_url):
        # For every alt in the region
        
        # Get its GI number
        alt_gi = get_gi_number(alt_acc)
        
        print("Downloading alt GI{}".format(alt_gi))
        
        # Grab the sequence data
        alt_seq = get_sequence(alt_gi)
        
        # Write it to <region>/GI<number>.fa
        SeqIO.write([alt_seq], open("{}/GI{}.fa".format(options.region, alt_gi),
            "w"), "fasta")
            
        # Add this alt to the chromosome-sizes-by-accession file
        acc_chrom_sizes.line(alt_acc, get_length(alt_gi))
        # Sneak into the TSV writer and flush, so the sizes file can now be
        # read.
        acc_chrom_sizes.stream.flush()
        
        # Where should we put the GFF alignment for this alt to the reference?
        alt_gff3 = "{}/GI{}.gff3".format(options.region, alt_gi)
        
        print("Downloading alignment")
        
        # Go download it
        download_gff3(ref_acc, alt_acc, alt_unit, options.assembly_url,
            alt_gff3)
        
        # And we need to convert that to PSL
        alt_psl = "{}/GI{}.psl".format(options.region, alt_gi)
        
        print("Converting to PSL")
        
        # Run the conversion with the bit of the sizes file we have so far. We
        # need to pass the chrom.sizes file twice now because gff3ToPsl has
        # changed its interface.
        subprocess.check_call(["gff3ToPsl", options.region + "/acc.chrom.sizes",
            options.region + "/acc.chrom.sizes", alt_gff3, alt_psl])
        
        # Edit the output to point to the GI instead of the accession
        subprocess.check_call(["sed", "-i", "s/{}/GI{}/g".format(alt_acc,
            alt_gi), alt_psl])
            
        print("Writing genes for GI{}".format(alt_gi))
    
        # Make a BED to put alt genes in
        alt_bed = open_gene_bed(options.region, "GI{}".format(alt_gi)) 
        
        for line in get_genes(alt_acc, "GI{}".format(alt_gi),
            alt_parent_grc_id=ref_acc):
            # Write all the BED lines for the alt to the file
            alt_bed.write(line)
            
        alt_bed.close()
       
    # Now we need to do psl2maf, complete with globbing.
            
    print("Creating GRC MAF")
            
    # Find the psl2maf.py script
    psl2maf = (os.path.dirname(os.path.realpath(__file__)) +
        "/../mhc/psl2maf.py")
        
    # Go call psl2maf, moving the reference stuff over to "ref" and shifting it
    # back so that the first base we clipped out of the reference is 0,
    # splitting apart mismatches, and making sure to use all the PSLs and MAFs
    # in our output directory. We make sure to add 1 to the reference start in
    # the offset, because some basedness-conversion needs to happen. TODO: Make
    # this a function or make this use an import or somehow de-uglify it.
    args = ([psl2maf, "--maf", 
        options.region + "/GRCAlignment.maf", "--referenceOffset", 
        str(-ref_start + 1), "--referenceSequence", "ref", "--noMismatch",
        "--psls"] + glob.glob(options.region + "/*.psl") + ["--fastas"] + 
        glob.glob(options.region + "/*.fa"))
        
    print("Calling: {}".format(" ".join(args)))
        
    subprocess.check_call(args)
        
        
            
        
            
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
