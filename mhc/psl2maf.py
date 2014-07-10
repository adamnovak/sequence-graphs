#!/usr/bin/env python2.7
"""
psl2maf.py: Convert a series of PSLs to a single multiple alignment MAF, pulling
sequence from FASTA files.

"""

from Bio import SearchIO, AlignIO, SeqIO, Align

import argparse, sys, os, os.path, random, subprocess, shutil, itertools


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
    parser.add_argument("--psl", type=str, required=True, 
        action="append",
        help="PSL file(s) to convert")
    parser.add_argument("--maf", type=argparse.FileType("w"), 
        help="MAF file to save output to")
    parser.add_argument("--fasta", required=True, action="append",
        help=".fasta file(s) to obtain sequence from")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Load all the FASTAs, indexed
    fastaDicts = [SeqIO.index(fasta, "fasta") for fasta in options.fasta]
    
    def getSequence(name):
        """
        Get a sequence by ID from the first FASTA that has it.
        """
        
        for fastaDict in fastaDicts:
            if fastaDict.has_key(name):
                return fastaDict[name]
        
        raise ValueError("No sequence {} in any FASTA".format(name))
    
    # Make a multiple alignment that will align all the sequences.
    totalMSA = None
        
    for psl in options.psl:
        # For each PSL we want in out MAF
        for result in SearchIO.parse(psl, "blat-psl"):
            # Parse the PSL and go through the results
            
            print result
            
            for hit in result:
                # For every hit in the result
                
                # Grab the query that matched in this hit (ends up being the
                # thing we hit with the alignment somehow)
                queryID = hit.query_id
                querySeqRecord = getSequence(queryID)
                
                # Grab the hit ID, which is the thing it hit, and the
                # sequence for that.
                hitID = hit.id
                hitSeqRecord = getSequence(hitID)
                
                # This is going to hold an MSA that we append onto
                global hitMSA
                hitMSA = None
                # Where does the alignment end in the query sequence?
                global hitMSAQueryPos
                hitMSAQueryPos = 0
                # And in the hit sequence?
                global hitMSAHitPos
                hitMSAHitPos = 0

                def appendMSA(newMSA, hitStart, queryStart):
                    """
                    Adjoin a new MultipleSeqAlignment, which starts at the given
                    index in the hit and the other given index in the query, to
                    the right of the current one for the hit, padding with
                    gapped/unaligned sequence if necessary.
                    
                    The MSA must have the hit first and the query second. If the
                    MSA is None, the alignment is padded out to the given
                    lrngths in each sequence.
                    
                    Indexes are assumed to be on whatever strand would make them
                    go up as new MSAs are added.
                    
                    TODO: Handle reverse strand right. It's not really possible
                    in MultipleSeqAlignment.
                    
                    """
                    
                    # TODO: This needs to be in some sort of positioned-partial-
                    # alignment object.
                    global hitMSA
                    global hitMSAQueryPos
                    global hitMSAHitPos
                    
                    # Pad with some gapped sequence to the right position in
                    # each to accept the new alignment.
                    
                    # How much padding do we need in each sequence?
                    queryPaddingNeeded = queryStart - hitMSAQueryPos
                    hitPaddingNeeded = hitStart - hitMSAHitPos
                    
                    print("Padding {} in query and {} in hit".format(
                        queryPaddingNeeded, hitPaddingNeeded))
                    
                    # Grab the sequence from the query to go opposite the
                    # hit gap
                    queryActualSequence = querySeqRecord[hitMSAQueryPos:
                        queryStart]
                        
                    # Grab the sequence from the hit to go opposite the
                    # query gap
                    hitActualSequence = hitSeqRecord[hitMSAHitPos:hitStart]
                    
                    # Make the padding alignment from those sequences,
                    # adding in the gaps.
                    paddingAlignment = Align.MultipleSeqAlignment([
                        queryActualSequence + "-" * hitPaddingNeeded, 
                        "-" * queryPaddingNeeded + hitActualSequence])
                    
                    if hitMSA is None:
                        # We're just starting. Start with the padding alignment.
                        hitMSA = paddingAlignment
                    else:
                        # Add in the padding alignment
                        hitMSA += paddingAlignment
                        
                    if newMSA is None:
                        # Just pad out to the given lengths.
                        hitMSAHitPos = hitStart
                        hitMSAQueryPos = queryStart
                    else:
                        # Add in the new alignment we actually wanted to add.
                        hitMSA += newMSA
                        
                        # Update the indicators of where we're up to
                        hitMSAHitPos = hitStart + len(newMSA[0])
                        hitMSAQueryPos = queryStart + len(newMSA[1])
                    
                    print("Now at {} in hit and {} in query".format(
                        hitMSAHitPos, hitMSAQueryPos))
                        
                    
                        
                        
                        
                
                for hsp in hit:
                    # For every HSP (high-scoring pair) in the hit
                    
                    for fragment in hsp:
                        # For every HSP fragment in the HSP (actual alignment
                        # bit. Why is it so insanely nested?)
                        
                        # Fix up the fragment by going and fetching its hit
                        # sequence piece from the appropriate SeqRecord.
                        hitFragment = hitSeqRecord[fragment.hit_start:
                            fragment.hit_end]
                            
                        if fragment.hit_strand == "-":
                            # We meant to get the other strand.
                            # TODO: Properly support this!
                            hitFragment = hitFragment.reverse_complement()
                        
                        # Make sure we got the right number of bases.
                        assert(len(hitFragment) == fragment.hit_span)
                        
                        # Put it in.
                        fragment.hit = hitFragment
                        
                        # Now grab the bit of the query sequence involved in
                        # this fragment.
                        queryFragment = querySeqRecord[fragment.query_start:
                            fragment.query_end]
                            
                        # Make sure we got the right number of bases.
                        assert(len(queryFragment) == fragment.query_span)
                            
                        if fragment.query_strand == "-":
                            # We meant to get the other strand
                            # TODO: Properly support this!
                            queryFragment = queryFragment.reverse_complement()
                            
                        # Put it in
                        fragment.query = queryFragment
                        
                        # Get the MultipleSeqAlignment which the fragment then
                        # creates.
                        alignment = fragment.aln
                        
                        # Stick that onto the end of our growing MSA for this
                        # hit. Padding to align this alignment to the right
                        # place in each sequence will be automatically added.
                        appendMSA(alignment, fragment.hit_start,
                            fragment.query_start)
                            
                # Pad out the MSA with the unaligned sequence in each sequence
                appendMSA(None, len(hitSeqRecord), len(querySeqRecord))
                print(hitMSA)
                        
                # Now we have completed an MSA for just this hit
                # TODO: Merge it into the master MSA keying on the reference.
                # For now just take the last one
                totalMSA = hitMSA
                
    # Save the multiple alignment as a MAF.
    AlignIO.write(totalMSA, options.maf, "maf")
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
