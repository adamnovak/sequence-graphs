#!/usr/bin/env python2.7
"""
psl2maf.py: Convert a series of PSLs to a single multiple alignment MAF, pulling
sequence from FASTA files. The PSLs must all be against the same reference
sequence, and the alignment to that reference is used to induce the multiple
alignment (i.e. nothing in the other sequences will align except as mediated by
alignment to the reference).

TODO: snake_case this file.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import doctest, logging, pprint

from Bio import SearchIO, AlignIO, SeqIO, Align
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
    parser.add_argument("--psls", nargs="+", required=True, 
        help="PSL file(s) to convert")
    parser.add_argument("--maf", type=argparse.FileType("w"),
        default=sys.stdout,
        help="MAF file to save output to")
    parser.add_argument("--fastas", nargs="+", required=True,
        help=".fasta file(s) to obtain sequence from")
    parser.add_argument("--referenceOffset", type=int, default=0,
        help="offset all reference coordinates by the given amount")
    parser.add_argument("--referenceSequence", type=str, default=None,
        help="override reference sequence name with this one")
    parser.add_argument("--noMismatch", action="store_true",
        help="only align bases which match")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def gapMismatches(alignment):
    """
    Given an alignment (an MSA with just a reference and a query), replace any
    mismatches with gaps in each sequence.
    
    Return the processed alignment.
    """
    
    # Make lists of characters that we will join into the new reference and
    # query sequences.
    gappedReference = []
    gappedQuery = []
    
    # How many mismatches did we gap?
    mismatches_gapped = 0
    # How many aligned bases did we check?
    bases_checked = 0
    
    
    # Where are we in the alignment? 
    for column in xrange(len(alignment[0])):
        # Just go through all the columns in the alignment's reference.
        
        # Pull out the reference and query characters at this position.
        refChar = alignment[0, column]
        queryChar = alignment[1, column]
        
        bases_checked += 1
        
        if "-" in [refChar, queryChar] or refChar == queryChar:
            # We have a gap or a match. Pass it through to bioth sequences.
            gappedReference.append(refChar)
            gappedQuery.append(queryChar)
        else:
            # We have a mismatch. Gap one and then the other.
            gappedReference.append("-")
            gappedQuery.append(queryChar)
            
            gappedReference.append(refChar)
            gappedQuery.append("-")
            
            mismatches_gapped += 1
            
    # Now we need to manufacture the MultipleSeqAlignment to return from these
    # lists of characters.

    # What names do the sequences in this alignment have?
    seqNames = [record.id for record in alignment]
    
    # Make a SeqRecord for each list of properly gapped-out characters, with the
    # appropriate name.
    seqRecords = [SeqRecord(Seq("".join(alignedList)), name) 
        for alignedList, name in zip([gappedReference, gappedQuery], seqNames)]
    
    for i in xrange(len(seqRecords)):
        # Se tannotations on all the new records
        seqRecords[i].annotations = alignment[i].annotations
    
    if float(mismatches_gapped) / bases_checked > 0.5 and bases_checked > 100:    
        # If this gets too high, it means we have a bad offset somewhere. Yell
        # at the user.
        logging.warning("{}/{} bases gapped due to mismatch".format(
            mismatches_gapped, bases_checked))
        
    # Make the records into a proper MSA and return it.
    return Align.MultipleSeqAlignment(seqRecords)
    
def tree_reduce(items, operator, default_value=None):
    """
    Reduce the items in the given list using the given binary function. Unlike
    the normal Python reduce, which reduces by a fold from the left, this reduce
    reduces in a tree. This means that if the operation is, say, a list
    concatenation, the reduction will be more efficient than if it were done in
    a left to right order, as each input element will only be copied log(n)
    times, instead of n for the leftmost element in a normal reduce.
    
    If items is empty, the given default_value is returned.
    
    """
    
    if len(items) == 0:
        # Base case for an empty list
        return None
    elif len(items) == 1:
        # Base case for a one-item list
        return items[0]
    else:
        # We have two or more items (and will never have to deal with 0)
        
        # Where is the middle? Do integer division.
        middle = len(items) / 2
        
        # Figure out the answers on the left and the right
        left_answer = tree_reduce(items[0:middle], operator)
        right_answer = tree_reduce(items[middle:], operator)
        
        logging.info("Merging:")
        logging.info(left_answer)
        logging.info("And:")
        logging.info(right_answer)
        
        # Recurse down both sides, do the operator to combine them, and return.
        to_return = operator(left_answer, right_answer)
        
        logging.info("Yields:")
        logging.info(to_return)
        
        return to_return
            
        
            
def smart_adjoin(msa1, msa2, sequence_source):
    """
    Given two Multiple Sequence Alignments (MSAs) on the same source sequences,
    with correct annotations, concatenate them together, with the intervening
    sequences unaligned.
    
    Either MSA may be None, in which case the other is returned.
    
    Requires a function that, when passed a sequence ID, returns the SeqRecord
    for the full sequence.
    
    Requires that there be a valid way to attach the two sequences together
    (i.e. the same sequence doesn't run in different directions in the two
    blocks).
    
    Raises a RuntimeError if the two MSAs cannot be adjoined.
    
    """
    
    if msa1 is None:
        # Nothing plus something equals that thing.
        return msa2
        
    if msa2 is None:
        # Nothing plus something equals that thing.
        return msa1
        
    logging.debug("Adjoining {}bp and {}bp reference alignments".format(
        msa1[0].annotations["size"], msa2[0].annotations["size"]))
    
    for seq1, seq2 in itertools.izip(msa1, msa2):
        # Check all the sequences
        
        if seq1.annotations["strand"] != seq2.annotations["strand"]:
            # These alignments are to opposite reference strands and cannot be
            # adjoined.
            raise RuntimeError("Can't adjoin alignments on opposite strands")
            
    if msa2[0].annotations["start"] < msa1[0].annotations["start"]:
        # Whatever strand we're on for the first sequence, alignment 2 needs to
        # happen first.
        msa2, msa1 = msa1, msa2
        
    # We're going to get the sequence needed to go from the end of MSA1 to the
    # start of MSA2.
    intervening_sequences = []
    
    for seq1, seq2 in itertools.izip(msa1, msa2):
        # For each pair of sequence pieces, we need the sequence from #1 to #2,
        # on the appropriate strand.
        
        # Where does the intervening sequence start along the strand in
        # question? Remember MAF coordinates are 0-based.
        intervening_start = seq1.annotations["start"] + seq1.annotations["size"]
        
        # And where does it end? (1 past the end)
        intervening_end = seq2.annotations["start"]
        
        if intervening_end < intervening_start:
            # We're always going up in strand-local coordinates.
            raise RuntimeError("Sequence is trying to go backwards!")
        
        if seq1.annotations["strand"] == -1:
            # Convert to the correct strand.
            
            intervening_start = seq1.annotations["srcSize"] - intervening_start
            intervening_end = seq1.annotations["srcSize"] - intervening_end
            
            intervening_start, intervening_end = (intervening_end, 
                intervening_start)
            
        # Go get and clip out the intervening sequence.    
        intervening_sequence = sequence_source(seq1.id)[
            intervening_start:intervening_end]
            
        if seq1.annotations["strand"] == -1:
            # Make sure it is on the correct strand
            intervening_sequence = intervening_sequence.reverse_complement()
            
        # Put the clipped-out, correctly-oriented unaligned sequence in the
        # list.
        intervening_sequences.append(intervening_sequence)
        
    # We'll tack these additional alignments onto msa1
    to_return = msa1
        
    for i in xrange(len(intervening_sequences)):
        # Now for each intervening sequence, I need an MSA consisting of that
        # sequence in its correct row and gaps in all the other rows.
        
        # Make all the rows for this bit of unaligned sequence, as SeqRecords.
        alignment_rows = [SeqRecord(Seq("-" * len(intervening_sequences[i])))
            if j != i else intervening_sequences[i]
            for j in xrange(len(intervening_sequences))]
            
        # Make them into an alignment and stick it on
        to_return = to_return + Align.MultipleSeqAlignment(alignment_rows)

    # Now stick on msa2
    to_return = to_return + msa2
    
    for i in xrange(len(to_return)):
        # Do the annotations for each record in the alignment
        
        # Set the ID
        to_return[i].id = msa1[i].id
        
        # Start with the annotations from msa1, so start is correct
        to_return[i].annotations.update(msa1[i].annotations)
        
        # Compute the actual sequence length that outght to be used here.
        to_return[i].annotations["size"] = (msa2[i].annotations["start"] + 
            msa2[i].annotations["size"] - msa1[i].annotations["start"])
            
        # Make sure size is correct correct.
        assert(len(str(to_return[i].seq).replace("-","")) == 
            to_return[i].annotations["size"])
    
    # Give back the final adjoined alignment
    return to_return
        
def reverse_msa(msa):
    """
    Given a MultipleSeqAlignment with MAF annotations, reverse-complement it,
    correcting the annotations.
    
    >>> ref1 = SeqRecord(Seq("AT-ATATAT"), "first")
    >>> ref1.annotations = {"strand": 1, "start": 0, "size": 8, "srcSize": 18}
    >>> alt1 = SeqRecord(Seq("ATAATATAT"), "second")
    >>> alt1.annotations = {"strand": -1, "start": 0, "size": 9, "srcSize": 9}
    >>> msa1 = Align.MultipleSeqAlignment([ref1, alt1])
    >>> rev = reverse_msa(msa1)
    
    >>> print(msa1)
    Alphabet() alignment with 2 rows and 9 columns
    AT-ATATAT first
    ATAATATAT second
    >>> print(rev)
    Alphabet() alignment with 2 rows and 9 columns
    ATATAT-AT first
    ATATATTAT second
    >>> pprint.pprint(rev[0].annotations)
    {'size': 8, 'srcSize': 18, 'start': 10, 'strand': -1}
    >>> pprint.pprint(rev[1].annotations)
    {'size': 9, 'srcSize': 9, 'start': 0, 'strand': 1}
    
    >>> rev2 = reverse_msa(rev)
    >>> print(rev2)
    Alphabet() alignment with 2 rows and 9 columns
    AT-ATATAT first
    ATAATATAT second
    >>> pprint.pprint(rev2[0].annotations)
    {'size': 8, 'srcSize': 18, 'start': 0, 'strand': 1}
    >>> pprint.pprint(rev2[1].annotations)
    {'size': 9, 'srcSize': 9, 'start': 0, 'strand': -1}
    
    
    """
    
    logging.debug("Reversing {}bp reference MSA".format(
        msa[0].annotations["size"]))
    
    # Make an alignment with all the sequences reversed.
    to_return = Align.MultipleSeqAlignment((record.reverse_complement() 
        for record in msa))
        
    for i in xrange(len(to_return)):
        # Fix up the annotations on each sequence
        
        # Start with the original annotations
        to_return[i].annotations.update(msa[i].annotations)
        
        # We need to flip the strand
        to_return[i].annotations["strand"] = -msa[i].annotations["strand"]
        
        # And count the start from the other end.
        to_return[i].annotations["start"] = (msa[i].annotations["srcSize"] - 
            msa[i].annotations["start"] - msa[i].annotations["size"])
            
        # Set the id
        to_return[i].id = msa[i].id
        
    # We finished it.
    return to_return
        

    
def mergeMSAs(msa1, msa2, full_ref):
    """
    Given two MultipleSeqAlignment objects sharing a first (reference) sequence,
    merge them on the reference. Returns a MultipleSeqAlignment containing all
    the sequences from each alignment, in the alignment induced by the shared
    reference sequence.
    
    Also needs access to the full reference SeqRecord in case it needs bases to
    fill in a gap.
    
    The first sequence may actually be only a subrange in either MSA, and either
    MSA may be on either strand of it.
    
    Either MSA may be None, in which case the other MSA is returned.
    
    >>> ref = SeqRecord(Seq("ATATATATGCATATATAT"), "first")
    >>> ref.annotations = {"strand": 1, "start": 0, "size": 18, "srcSize": 18}
    >>> ref1 = SeqRecord(Seq("AT-ATATAT"), "first")
    >>> ref1.annotations = {"strand": 1, "start": 0, "size": 8, "srcSize": 18}
    >>> alt1 = SeqRecord(Seq("ATAATATAT"), "second")
    >>> alt1.annotations = {"strand": -1, "start": 0, "size": 9, "srcSize": 9}
    >>> ref2 = SeqRecord(Seq("ATATATAT--"), "first")
    >>> ref2.annotations = {"strand": -1, "start": 0, "size": 8, "srcSize": 18}
    >>> alt2 = SeqRecord(Seq("ATATGG--AT"), "third")
    >>> alt2.annotations = {"strand": 1, "start": 0, "size": 8, "srcSize": 8}
    
    >>> msa1 = Align.MultipleSeqAlignment([ref1, alt1])
    >>> msa2 = Align.MultipleSeqAlignment([ref2, alt2])
    >>> merged = mergeMSAs(msa1, msa2, ref)
    >>> print(merged)
    Alphabet() alignment with 3 rows and 21 columns
    AT-ATATATGC--ATATATAT first
    ATAATATAT------------ second
    -----------AT--CCATAT third
    >>> pprint.pprint(merged[0].annotations)
    {'size': 18, 'srcSize': 18, 'start': 0, 'strand': 1}
    >>> pprint.pprint(merged[1].annotations)
    {'size': 9, 'srcSize': 9, 'start': 0, 'strand': -1}
    >>> pprint.pprint(merged[2].annotations)
    {'size': 8, 'srcSize': 8, 'start': 0, 'strand': -1}
    
    
    >>> ref3 = SeqRecord(Seq("ATGCAT"), "first")
    >>> ref3.annotations = {"strand": 1, "start": 6, "size": 6, "srcSize": 18}
    >>> alt3 = SeqRecord(Seq("ATCCAT"), "fourth")
    >>> alt3.annotations = {"strand": 1, "start": 5, "size": 6, "srcSize": 15}
    >>> msa3 = Align.MultipleSeqAlignment([ref3, alt3])
    
    >>> merged2 = mergeMSAs(merged, msa3, ref)
    >>> print(merged2)
    Alphabet() alignment with 4 rows and 21 columns
    AT-ATATATGC--ATATATAT first
    ATAATATAT------------ second
    -----------AT--CCATAT third
    -------ATCC--AT------ fourth
    >>> pprint.pprint(merged2[0].annotations)
    {'size': 18, 'srcSize': 18, 'start': 0, 'strand': 1}
    >>> pprint.pprint(merged2[1].annotations)
    {'size': 9, 'srcSize': 9, 'start': 0, 'strand': -1}
    >>> pprint.pprint(merged2[2].annotations)
    {'size': 8, 'srcSize': 8, 'start': 0, 'strand': -1}
    >>> pprint.pprint(merged2[3].annotations)
    {'size': 6, 'srcSize': 15, 'start': 5, 'strand': 1}
    
    

    
    """
    
    if msa1 is None:
        # No merging to do.
        return msa2
    
    if msa2 is None:
        # No merging to do this way either.
        return msa1
        
    if msa1[0].annotations["strand"] == -1:
        # MSA 1 needs to be on the + strand of the reference
        msa1 = reverse_msa(msa1)
        
    if msa2[0].annotations["strand"] == -1:
        # MSA 2 also needs to be on the + strand of the reference
        msa2 = reverse_msa(msa2)
        
    if msa2[0].annotations["start"] < msa1[0].annotations["start"]:
        # msa2 starts before msa1. We want msa1 to start first, so we need to
        # flip them.
        msa1, msa2 = msa2, msa1
        
    logging.debug("Zipping {}bp/{} sequence and {}bp/{} sequence  reference "
        "alignments".format(msa1[0].annotations["size"], len(msa1),
        msa2[0].annotations["size"], len(msa2)))
        
    # Make sure we are joining on the right sequence.
    assert(msa1[0].id == msa2[0].id)
        
    logging.debug("Merging")
        
    logging.debug(msa1)
    logging.debug(msa1[0].annotations)
    logging.debug(msa2)
    logging.debug(msa2[0].annotations)
        
    # Compute the offset: number of extra reference columns that msa2 needs in
    # front of it. This will always be positive or 0.
    msa2_leading_offset = (msa2[0].annotations["start"] - 
        msa1[0].annotations["start"])
    
    logging.debug("{}bp between left and right alignment starts".format(
        msa2_leading_offset))
        
    # It would be nice if we could shortcut by adjoining compatible alignments,
    # but the IDs wouldn't match up at all.
    
    # Make lists for each sequence we are going to build: those in msa1, and
    # then those in msa2 (except the duplicate reference).
    merged = [list() for i in xrange(len(msa1) + len(msa2) - 1)]
    
    # Start at the beginning of both alignments.
    msa1Pos = 0
    msa2Pos = 0
    
    # How many reference characters have been used?
    refChars = 0
    
    while refChars < msa2_leading_offset and msa1Pos < len(msa1[0]):
        # Until we're to the point that MSA 2 might have anything to say, we
        # just copy MSA 1.
        
        for i, character in enumerate(msa1[:, msa1Pos]):
            # For each character in the first alignment in this column
            
            # Put that character as the character for the appropriate
            # sequence.
            merged[i].append(character)
            
        for i in xrange(len(msa1), len(msa1) + len(msa2) - 1):
            # For each of the alignment rows that come from msa2, put a gap.
            merged[i].append("-")
            
        
        if msa1[0, msa1Pos] != "-":
            # We consumed a reference character.
            refChars += 1
            
        # We used some of MSA1
        msa1Pos += 1
        
    logging.debug("Used {}/{} offset".format(refChars, msa2_leading_offset))
        
    while refChars < msa2_leading_offset:
        # We have a gap between the first MSA and the second, and we need to
        # fill it with reference sequence.
        
        # We know we are refChars after the beginning of the first reference, so
        # we use that to know what base to put here.
        merged[0].append(full_ref[msa1[0].annotations["start"] + refChars])
        
        for i in xrange(1, len(msa1) + len(msa2) - 1):
            # And gap out all the other sequences.
            merged[i].append("-")
            
        # We consumed (or made up) a reference character
        refChars += 1
    
    while msa1Pos < len(msa1[0]) and msa2Pos < len(msa2[0]):
        # Until we hit the end of both sequences
        
        if refChars % 10000 == 0:
            logging.debug("Now at {} in alignment 1, {} in alignment 2, {} in "
                "reference".format(msa1Pos, msa2Pos, refChars))
        
        if(msa1[0, msa1Pos] == "-"):
            # We have a gap in the first reference. Put this column from the
            # first alignment alongside a gap for every sequence in the second
            # alignment.
            for i, character in enumerate(msa1[:, msa1Pos]):
                # For each character in the first alignment in this column
                
                # Put that character as the character for the appropriate
                # sequence.
                merged[i].append(character)
                
            for i in xrange(len(msa1), len(msa1) + len(msa2) - 1):
                # For each of the alignment rows that come from msa2, put a gap.
                merged[i].append("-")
                
            # Advance in msa1. We'll keep doing this until it doesn't have a gap
            # in its reference.
            msa1Pos += 1
            
        elif(msa2[0, msa2Pos] == "-"):
            # We have a letter in the first reference but a gap in the second.
            # Gap out the merged reference and all the columns from alignment 1,
            # and take the non-reference characters from alignment 2.
            
            for i in xrange(len(msa1)):
                # For the reference and all the sequences in msa1, add gaps
                merged[i].append("-")
           
            for i, character in zip(xrange(len(msa1),
                len(msa1) + len(msa2) - 1), msa2[1:, msa2Pos]):
                
                # For each of the alignment rows that come from msa2, put the
                # character from that row.
                merged[i].append(character)
                
            # Advance in msa2. We'll keep doing this until both msa1 and msa2
            # have a non-gap character in their references. We make it an
            # invariant that this will always be the same character.
            msa2Pos += 1
            
        else:
            # Neither has a gap. They both have real characters.
            
            if(msa1[0, msa1Pos] != msa2[0, msa2Pos]):
                logging.error(msa1)
                logging.error(msa2)
                raise RuntimeError("{} in reference 1 does not match {} "
                    "in reference 2".format(msa1[0, msa1Pos], msa2[0, msa2Pos])) 
            
            for i, character in enumerate(msa1[:, msa1Pos]):
                # Copy all the characters from msa1's column
                merged[i].append(character)
                
            for character, i in zip(msa2[1:, msa2Pos], xrange(len(msa1), 
                len(msa1) + len(msa2) - 1)):
                # Copy all the characters from msa2's column, except its
                # reference
                merged[i].append(character)
                
            # Advance both alignments
            msa1Pos += 1
            msa2Pos += 1
            
            # Say we used a reference character
            refChars += 1
            
        for otherMerged in merged[1:]:
            # Make sure we aren't dropping characters anywhere.
            assert(len(otherMerged) == len(merged[0]))
           
    logging.debug("At {}/{} of msa2, {}/{} of msa1".format(msa2Pos,
        len(msa2[0]), msa1Pos, len(msa1[0])))
        
    # By here, we must have finished one of the MSAs. Only one can have anything
    # left.
    assert(msa1Pos == len(msa1[0]) or msa2Pos == len(msa2[0]))
            
    while msa1Pos < len(msa1[0]):
        # MSA2 finished first and now we have to finish up with the tail end of
        # MSA1
        
        for i, character in enumerate(msa1[:, msa1Pos]):
            # For each character in the first alignment in this column
            
            # Put that character as the character for the appropriate
            # sequence.
            merged[i].append(character)
            
        for i in xrange(len(msa1), len(msa1) + len(msa2) - 1):
            # For each of the alignment rows that come from msa2, put a gap.
            merged[i].append("-")
            
        # Advance in msa1, until we finish it.
        msa1Pos += 1
            
    while msa2Pos < len(msa2[0]):
        # MSA1 finished first and now we have to finish up with the tail end of
        # MSA2
        
        # For the reference, put whatever it has in MSA2
        merged[0].append(msa2[0][msa2Pos])
        
        for i in xrange(1, len(msa1)):
            # For all the sequences in msa1, add gaps
            merged[i].append("-")
       
        for i, character in zip(xrange(len(msa1),
            len(msa1) + len(msa2) - 1), msa2[1:, msa2Pos]):
            
            # For each of the alignment rows that come from msa2, put the
            # character from that row.
            merged[i].append(character)
            
        # Advance in msa2, until we finish it.
        msa2Pos += 1
        
    # Now we have finished populating these aligned lists. We need to make a
    # MultipleSeqAlignment from them.
    
    # What names do the sequences in this alignment have? All the ones from
    # msa1, and then all the ones from msa2 except the first (which is the
    # reference)
    seqNames = [record.id for record in msa1] + [record.id 
        for record in msa2[1:]]
    
    # Make a SeqRecord for each list of properly gapped-out characters, with the
    # appropriate name.
    seqRecords = [SeqRecord(Seq("".join(alignedList)), name) 
        for alignedList, name in zip(merged, seqNames)]
        
    # Make the records into a proper MSA
    merged = Align.MultipleSeqAlignment(seqRecords)
    
    # Do the annotations for the reference
    merged[0].annotations.update(msa1[0].annotations)
    # Calculate the total reference bases used. It will be the distance between
    # the rightmost alignment end and the start of msa1, along the reference.
    merged[0].annotations["size"] = (max(msa2[0].annotations["start"] + 
        msa2[0].annotations["size"], msa1[0].annotations["start"] +
        msa1[0].annotations["size"]) - msa1[0].annotations["start"])
    
    for i in xrange(1, len(msa1)):
        # Copy over annotations from MSA1
        merged[i].annotations.update(msa1[i].annotations)
        
    for i in xrange(len(msa1), len(msa1) + len(msa2) - 1):
        # Copy over annotations from MSA2, starting after the reference.
        merged[i].annotations.update(msa2[i - len(msa1) + 1].annotations)
        
    # The merged result reverence needs to be longer than the input references.
    #assert(len(merged[0]) >= len(msa1[0]))
    #assert(len(merged[0]) >= len(msa2[0]))
        
    # Give back the merged MSA
    return merged
                
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    
    if len(args) == 2 and args[1] == "--test":
        # Turn on debug logging
        logging.basicConfig(level=logging.DEBUG)
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
        
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Load all the FASTAs, indexed
    fastaDicts = [SeqIO.index(fasta, "fasta") for fasta in options.fastas]
    
    # Cache sequences from getSequence. Will waste memory but save disk IO.
    cache = {}
    
    def getSequence(name):
        """
        Get a sequence by ID from the first FASTA that has it.
        """
        
        if cache.has_key(name):
            return cache[name]
        
        for fastaDict in fastaDicts:
            if fastaDict.has_key(name):
                cache[name] = fastaDict[name]
                return cache[name]
        
        raise ValueError("No sequence {} in any FASTA".format(name))
    
    # Save them all here
    all_msas = []
        
    for psl in options.psls:
        logging.debug("Processing {}".format(psl))
    
        # For each PSL we want in out MAF
        for result in SearchIO.parse(psl, "blat-psl"):
            # Parse the PSL and go through the results
            
            for hit in result:
                # For every hit in the result
                
                if options.referenceSequence is not None:
                    # Override the reference (i.e. query) sequence first. 
                    
                    # TODO: why is reference actually query and query actually
                    # hit???
                    hit.query_id = options.referenceSequence
                
                # Grab the query that matched in this hit (ends up being the
                # thing we hit with the alignment somehow)
                queryID = hit.query_id
                querySeqRecord = getSequence(queryID)
                
                # Grab the hit ID, which is the thing it hit, and the
                # sequence for that.
                hitID = hit.id
                hitSeqRecord = getSequence(hitID)
                
                for hsp in hit:
                    # For every HSP (high-scoring pair) in the hit (which
                    # corresponds to a PSL line, and will be on a consistent
                    # strand)
                    
                    # Let's make an MSA describing the entire HSP.
                    hsp_msa = None
                    
                    logging.info("Starting a new HSP")
                    
                    for fragment in hsp:
                        # For every HSP fragment in the HSP (actual alignment
                        # bit. Why is it so insanely nested?)
                        
                        # Offset the reference (i.e. query) coordinates
                        fragment.query_start += options.referenceOffset
                        fragment.query_end += options.referenceOffset
                        
                        # Fix up the fragment by going and fetching its hit
                        # sequence piece from the appropriate SeqRecord.
                        hitFragment = hitSeqRecord[fragment.hit_start:
                            fragment.hit_end]
                            
                        if fragment.hit_strand == -1:
                            # We meant to get the other strand.
                            hitFragment = hitFragment.reverse_complement()
                        
                        # Make sure we got the right number of bases.
                        assert(len(hitFragment) == fragment.hit_span)
                        
                        # Work out the start on the appropriate strand
                        if fragment.hit_strand == 1:
                            # On the forward strand we can use the normal start
                            hit_start = fragment.hit_start
                        else:
                            # We have to calculate the start index on the
                            # reverse strand. Do it for 0-based coordinates.
                            hit_start = (len(hitSeqRecord) - 
                                fragment.hit_start - fragment.hit_span)
                        
                        # Annotate the hit (alt) with strand, start, size, and
                        # srcSize, as MafIO uses
                        hitFragment.annotations = {
                            "strand": fragment.hit_strand,
                            "size": fragment.hit_span,
                            "start": hit_start,
                            "srcSize": len(hitSeqRecord)
                        }
                        
                        # Put it in.
                        fragment.hit = hitFragment
                        
                        # Now grab the bit of the query sequence involved in
                        # this fragment.
                        queryFragment = querySeqRecord[fragment.query_start:
                            fragment.query_end]
                            
                        if fragment.query_strand == -1:
                            # We meant to get the other strand
                            queryFragment = queryFragment.reverse_complement()
                            
                        # Make sure we got the right number of bases.
                        if len(queryFragment) != fragment.query_span:
                            raise RuntimeError("Query fragment has {} bases "
                                "instead of {}".format(len(queryFragment),
                                fragment.query_span))
                            
                        # Work out the start on the appropriate strand
                        if fragment.query_strand == 1:
                            # On the forward strand we can use the normal start
                            query_start = fragment.query_start
                        else:
                            # We have to calculate the start index on the
                            # reverse strand. Do it for 0-based coordinates.
                            query_start = (len(querySeqRecord) - 
                                fragment.query_start - fragment.query_span)
                                
                        # Annotate the query (ref) with strand, start, size, and
                        # srcSize, as MafIO uses
                        queryFragment.annotations = {
                            "strand": fragment.query_strand,
                            "size": fragment.query_span,
                            "start": query_start,
                            "srcSize": len(querySeqRecord)
                        }
                            
                        # Put it in
                        fragment.query = queryFragment
                        
                        # Get the MultipleSeqAlignment which the fragment then
                        # creates. Query (ref) is first.
                        alignment = fragment.aln

                        if options.noMismatch:
                            # We only want to have match operations in our
                            # alignment. If two bases don't match, we need to
                            # gap them apart.
                            alignment = gapMismatches(alignment)
                            
                        # Add in the alignment for this hit to its MSA, doing
                        # whatever reordering is needed to make it work.
                        hsp_msa = smart_adjoin(hsp_msa, alignment, getSequence)
                        
                    # Save the HSP MSA
                    all_msas.append(hsp_msa)
                        
                    logging.info("Produced MSA:")
                    logging.info(hsp_msa)
                    
    logging.info("Merging MSAs...")
                    
    # Now merge all the MSAs together with a tree reduce.
    combined_msa = tree_reduce(all_msas, lambda a, b: mergeMSAs(a, b,
        getSequence(a[0].id)))
    
    logging.info("Writing output")
    
    # Save all the alignments in one MAF.
    AlignIO.write(combined_msa, options.maf, "maf")
            

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
