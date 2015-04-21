#!/usr/bin/env python2.7
"""
evaluateMappings.py: evaluate a HAL alignment with a collection of evaluations
and metrics.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections, csv, tempfile, xml.etree

import doctest, pprint

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
    parser.add_argument("hal",
        help="HAL file to evaluate")
    parser.add_argument("--truth",
        help=".maf file of a true alignment for precision and recall")
    parser.add_argument("--beds", nargs="*",
        help=".bed file(s) of genes on the genomes in the HAL")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def parse_halstats_stats(stream):
    """
    Parse the output of halStats on a HAL file, from the given iterable.
    
    Returns a list of dicts from halStats column header (GenomeName,
    NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments) to
    string or int value as appropriate.
    
    >>> lines = ["",
    ... "hal v2.1",
    ... "(GI568335879,ref)rootSeq;",
    ... "",
    ... "GenomeName, NumChildren, Length, NumSequences, NumTopSegments, "
    ... "NumBottomSegments",
    ... "rootSeq, 10, 5524595, 1, 0, 67141",
    ... "GI568335879, 0, 4672374, 1, 25067, 0",
    ... "ref, 0, 4970457, 1, 49220, 0",
    ... ""]
    >>> parsed = list(parse_halstats_stats(lines))
    >>> len(parsed)
    3
    >>> parsed[1]["GenomeName"]
    'GI568335879'
    >>> parsed[0]["Length"]
    5524595

    """
    
    # Get a reader and skip the first 4 lines
    reader = csv.reader(stream)
    for i in xrange(4):
       next(reader)
        
    # Read the table header with the field names, which is the next nonempty
    # line.
    header = next(reader)
    
    # Throw out whitespace
    header = [x.strip() for x in header]
    
    for data in reader:
        if len(data) == 0:
            # We finished the table
            break
        
        elif len(data) != len(header):
            raise RuntimeError("Got {} fields for {} columns".format(len(data),
                len(header)))
                
        # Drop all the whitespace
        data = [x.strip() for x in data]
            
        # We'll fill in a dict with all the type-auto-detected values for each
        # column. We convert to int if a field is all digits (which is fine
        # since nothing will be negative).
        record = dict(zip(header, [int(x) if x.isdigit() else x for x in data]))
            
        # Yield each record as we parse it
        yield record
        
def get_halstats_stats(hal_filename):
    """
    Go get the halStats results, and return a list of dicts of stats, one per
    genome in the HAL.
    
    The halStats program must be on the PATH.
    
    """
    
    # Build the halStats command line
    args = ["halStats", hal_filename]
    
    # Open the process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
        
    # Parse the results
    stats_list = list(parse_halstats_stats(process.stdout))
            
    if process.wait() != 0:
        raise RuntimeError("halStats failed")
            
    # We are done with this process.
    process.stdout.close()
    
    return stats_list
    
def parse_halstats_coverage(stream):
    """
    Parse the output of halStats --coverage on a HAL file, from the given
    iterable.
    
    Returns a dict from alignment partner name to total number of bases in that
    alignment partner covered (i.e. the sum of the at least once, at least
    twice, at least thrice, etc. numbers.)
    
    >>> lines = [
    ... "Genome, # of sites mapping at least once, twice, thrice, ...",
    ... "ref, 81189, 11, 5",
    ... "GI262359905, 81187, 22, 3, 1",
    ... "GI528476558, 81007, 1"]
    >>> pprint.pprint(parse_halstats_coverage(lines))
    {'GI262359905': 81213, 'GI528476558': 81008, 'ref': 81205}
    
    """
    
    # Get a reader and skip the first line
    reader = csv.reader(stream)
    for i in xrange(1):
       next(reader)
        
    # This is the dict we will populate
    to_return = {}
        
    for data in reader:
        if len(data) == 0:
            # We finished the table
            break
        
        
        # Drop all the whitespace
        data = [x.strip() for x in data]
            
        # Sum up all the numbers in the row and save them under the contents of
        # the first column. These numbers may be big enough to need exponents.
        to_return[data[0]] = sum([float(x) for x in data[1:]])
    
    # Return the populated dict
    return to_return
    
def get_halstats_coverage(hal_filename, genome):
    """
    Get the number of bases that each other genome has aligned to the given
    genome in the given HAL file.
    
    The halStats program must be on the PATH.
    
    """
    
    # First we need to get the coverage
    args = ["halStats", hal_filename, "--coverage", genome]

    # Open the process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
    
    # Keep the dict from name to total columns shared with the reference.
    coverage_dict = parse_halstats_coverage(process.stdout)
            
    if process.wait() != 0:
        raise RuntimeError("halStats --coverage failed")
        
    return coverage_dict
    
def parse_halstats_basecomp(stream):
    r"""
    Parse the output of halStats --baseComp on a HAL file and a genomne, from
    the given iterable.
    
    Returns the total fraction of the genome tested that is not Ns, clamped to
    1.
    
    >>> lines = ["0.290273\t0.220496\t0.208913\t0.280317"]
    >>> print(parse_halstats_basecomp(lines))
    0.999999
    
    >>> lines = ["0.290273\t0.220496\t0.208913\t0.280319"]
    >>> print(parse_halstats_basecomp(lines))
    1.0
    
    """
    
    # Grab the first (only) line
    line = next(iter(stream))
    
    # Split it by tabs, parse the numbers, sum up, and min with 1.
    return min(1.0, sum([float(x) for x in line.split("\t")]))
    
def get_halstats_basecomp(hal_filename, genome):
    """
    Get the portion of the given genome that is not N in the given HAL, clamped
    to 1.
    
    """
    
    # First we need to write the args.
    args = ["halStats", hal_filename, "--baseComp", "{},1".format(genome)]

    # Open the process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
    
    # Parse out the portion that is not N.
    non_n = parse_halstats_basecomp(process.stdout)
            
    if process.wait() != 0:
        raise RuntimeError("halStats --baseComp failed")
        
    return non_n

        
def metric_halstats(hal_filename, reference_id="ref"):
    """
    A metric that runs halStats on the given HAL, and returns the results.
    
    Results include coverage against the genome specified by reference_id, but
    ignoring bases in each other genome that are N and thus cannot be aligned.
    
    The halStats program must be on the PATH.
    
    """
    
    # Get the list of dicts of per-genome stats.
    status_list = get_halstats_stats(hal_filename)
    
    # Get the dict from genome name to total columns shared with the reference.
    coverage_dict = get_halstats_coverage(hal_filename, reference_id)
            
    for entry in status_list:
        # For each genome, we want the coverage against the reference.
        
        # Grab the genome name
        genome_name = entry["GenomeName"]
        
        if not coverage_dict.has_key(genome_name):
            # This is probably the root sequence and didn't get a covrage for
            # some reason. At any rate, the root sequence would be all Ns
            continue
        
        # Figure out how much of it is not Ns
        non_n = get_halstats_basecomp(hal_filename, genome_name)
        
        # How many bases are eligible?
        eligible = float(entry["Length"] * non_n)
        
        if eligible == 0:
            # No coverage is defined
            entry["Coverage"] = float("NaN")
            continue
    
        # Compute and save the coverage for each entry, by dividing bases
        # aligned by bases eligible.
        entry["Coverage"] = coverage_dict[genome_name] / eligible
            
    # Return the results
    return status_list
    
def hal2maf(hal_filename):
    """
    Convert the HAL at the given filename to a MAF in a temporary file. Returns
    the name of the temp file, which must be deleted by the caller.
    
    hal2maf must be installed and on the PATH.
    
    """
    
    # Make a temporary file
    handle, maf_filename = tempfile.mkstemp()
    os.close(handle)
    
    # Run the conversion
    subprocess.check_call(["hal2maf", hal_filename, maf_filename])
    
    # Return the filename of the MAF
    return maf_filename
    
def parse_mafcomparator_output(stream):
    """
    Given a stream of XML output from mafComparator, between a truth MAF and a
    MAF under test (in that order), parse the XML and produce a (precision,
    recall) tuple.
    
    """
    
    # Parse the XML output
    tree = xml.etree.ElementTree.parse(stream)
    
    # Grab the nodes with the aggregate results for each comparison direction
    stats = tree.findall(
        "./homologyTests/aggregateResults/all")
        
    # Grab and parse the averages. There should be two.
    averages = [float(stat.attrib["average"]) for stat in stats]
    
    # Return them in precision, recall order. First we put portion of MAF under
    # test represented in truth (precision), then portion of truth represented
    # in MAF under test (recall).
    return (averages[1], averages[0])    

def metric_mafcomparator(maf_filename, maf_truth):
    """
    Run mafComparator to compare the given MAF against the given truth MAF.
    
    Returns a precision, recall tuple.
    
    """

    # Make a temporary file
    # Make a temporary file
    handle, xml_filename = tempfile.mkstemp()
    os.close(handle)
    

    # Set up the mafComparator arguments. We sample only a subset of the
    # columns, but still enough that we will probably get a good number of
    # sig figs on coverage.
    
    # Use a hardcoded seed for consistency
    args = ["mafComparator", "--maf1", maf_truth, "--maf2",
        maf_filename, "--out", xml_filename, "--seed", "1", "--samples",
        "1000000"]
    
    # Run mafComparator and generate the output XML
    subprocess.check_call(args)
    
    # Read the output and get precision and recall
    stats = parse_mafcomparator_output(open(xml_filename))
    
    # Delete the temporary file
    os.unlink(xml_filename)
    
    # Return the precision and recall
    return stats
    
    
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Print the halStats output
    pprint.pprint(metric_halstats(options.hal))
    
    if options.beds is not None or options.truth is not None:
        # We need a MAF for checkGenes and for the precision/recall
        # calculations.
        maf_filename = hal2maf(options.hal)
    
    if options.beds is not None:
        
        # We're going to check the alignment agains the genes
        import checkGenes
        class_counts, gene_sets, gene_pairs = checkGenes.check_genes(
            maf_filename, options.beds)
            
        # Print the output
        pprint.pprint(class_counts)
        pprint.pprint(gene_sets)
        pprint.pprint(gene_pairs)
        
    if options.truth is not None:
    
        # We're going to get precision and recall against the truth.
        precision, recall = metric_mafcomparator(maf_filename, options.truth)

        # TODO: Output better        
        print(precision, recall)
        
    if options.beds is not None or options.truth is not None:
        # Clean up the MAF
        os.unlink(maf_filename)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
