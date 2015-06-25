#!/usr/bin/env python2.7
"""
unitigStats.py: Read a (possibly gzipped) unitig .mag file, put it in a
database, and report some statistics about it.

MAG format is a variant of FASTQ format, with specially formatted header lines
and quality scores.

The header of a record contains four fields, separated by tab characters. The
first field gives the numerical identifiers for the two ends of the sequence,
separated by a colon. The second field gives the number of reads that went into
the record. The third field, if it is not ".", gives semicolon-terminated end-
id,offset pairs detailing the end with which the left end of this record
overlaps, and the number of bases of overlap. The fourth field is the same, but
for the right end of this sequence.

Example record header line:

    @3056104052:2484561266	1040	1855531145,98;	743607148,86;183408,88;

The qualities in the FASTQ record give the number of reads covering each base,
in chr(33 + value) format.

Note that the + line for qualities will not contain the same header information.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, zlib, time, collections, hashlib

import Bio.SeqIO
import Bio.SeqRecord
import sqlite3


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
    parser.add_argument("--in_file", type=argparse.FileType("r"),
        default=sys.stdin,
        help="MAG-format file of the unitig graph to read")
    parser.add_argument("database",
        help="Database filename to write")
    parser.add_argument("fasta",
        help="FASTA filename to write")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

class TransparentUnzip(object):
    """
    A class that represents a transparently-un-gzipping stream.
    
    Right now only supports readline
    """
    
    def __init__(self, stream):
        """
        Encapsulate the given stream in an unzipper.
        """
        
        # Keep the stream
        self.stream = stream
        
        # Make a decompressor with the accept a gzip header flag.
        # See <http://stackoverflow.com/a/22311297/402891>
        self.decompressor = zlib.decompressobj(zlib.MAX_WBITS + 16)
        
        # We need to do lines ourselves, so we need to keep a buffer
        self.line_buffer = ""
        
        # Track out throughput
        self.compressed_bytes = 0
        self.uncompressed_bytes = 0
        
    def readline(self):
        """
        Return the next line with trailing "/n", or "" if there is no next line.
        """
        
        # See if we have a line to spit out
        newline_index = self.line_buffer.find("\n")
        
        if newline_index == -1:
            # No line is in the buffer
            # Go get more data from the stream (in 16 k blocks)
            compressed = self.stream.read(16 * 2 ** 10)
            
            # Track the amount read
            self.compressed_bytes += len(compressed)
            
            if compressed == "":
                if self.decompressor is not None:
                    # No mor einput data; flush out the output data.
                    self.line_buffer += self.decompressor.flush()
                    self.decompressor = None
                    # Spit out a line from that
                    return self.readline()
                
                # We didn't find a newline, and there's no more data, and
                # nothing to flush. Take what we have (with no trailing \n)
                line = self.line_buffer
                # Clear the buffer so next time we return "" as desired.
                self.line_buffer = ""
                
                # Track the uncompressed bytes
                self.uncompressed_bytes += len(line)
                
                return line
                
            # Otherwise we found more data
            
            try:
                # Decompress each block if needed
                decompressed = self.decompressor.decompress(compressed)
            except zlib.error:
                # Just skip decompressing; it's probably not actually
                # compressed.
                decompressed = compressed
            
            # Stick it in the buffer
            self.line_buffer += decompressed
            
            # Try again
            return self.readline()
        
        else:
            # We have a line. Grab it.
            line = self.line_buffer[0:newline_index + 1]
            
            # Pop it off
            self.line_buffer = self.line_buffer[newline_index + 1:]
            
            # Track the uncompressed bytes
            self.uncompressed_bytes += len(line)
            
            # Return it
            return line
            
    def start_stats(self):
        """
        Reset byte stat tracking
        """
        
        self.compressed_bytes = 0
        self.uncompressed_bytes = 0
        
    def get_stats(self):
        """
        Return the number of compressed, uncompressed bytes processed.
        """
        
        return self.compressed_bytes, self.uncompressed_bytes        
        
def parse_mag(stream):
    """
    Yield MAG records from the given (possibly compressed) FASTQ stream.
    
    Yields records of the form [[left end id, right end id], read count,
    [[left end overlap end, left end overlap length], ...], 
    [[right end overlap end, right end overlap length], ...], sequence]
    
    """
    
    for fastq_record in Bio.SeqIO.parse(stream, "fastq"):
        # For each FASTQ record
        
        # Parse the header line
        parts = fastq_record.description.split("\t")
        
        # Parse the ends
        parts[0] = [int(end_id) for end_id in parts[0].split(":")]
        
        # Parse the read count
        parts[1] = int(parts[1])
        
        # Parse the overlaps
        for i in xrange(2,4):
            # We have left end and right end overlaps we want to treat the same
            # way
            if parts[i] == ".":
                # This means no overlaps
                parts[i] = []
            else:
                # There are some overlaps
                
                # Grab the ;-terminated overlaps
                overlaps = parts[i].split(";")[:-1]
                
                # Split each overlap on the comma, and int the parts
                parts[i] = [[int(x) for x in overlap.split(",")]
                    for overlap in overlaps]
                    
        # Stick in the sequence
        parts.append(fastq_record.seq)
        
        # We've filled in the fields here, so we can spit this out.
        yield parts

class Sequence(object):
    """
    Represents a piece of sequence in a side graph.
    
    """
    
    def __init__(self, sequence_id, length):
        """
        Make a new sequence with the given ID of the given length.
        """
        
        # Save the arguments
        self.id = sequence_id
        self.length = length
   
class Side(object):
    """
    Represents an oriented position on a Sequence.
    """
    
    def __init__(self, sequence_id, offset, is_right):
        """
        Make a new Position on the given Sequence, at the given base offset, and
        either on the left or right face of that base as specified.
        
        """
        
        # Save the arguments
        self.sequence_id = sequence_id
        self.offset = offset
        self.is_right = is_right
        
    def add_local_offset(self, offset):
        """
        Produce a new Side by moving this Side the specified number of bases
        outward along its strand.
        """
        
        if not self.is_right:
            # If we're a left side, outward is negative. 
            offset = -offset
        
        return Side(self.sequence_id, self.offset + offset, self.is_right)
        
    def flip(self):
        """
        Produce a new Side for the opposite side of this Side's base.
        """
        
        return Side(self.sequence_id, self.offset, not self.is_right)
        
    def __repr__(self):
        """
        Return a string representation of this Side.
        """
        
        return "Side({}, {}, {})".format(self.sequence_id, self.offset,
            self.is_right)
            
    def __lt__(self, other):
        """
        Return ture if this Side is less than the other one, and false
        otherwise.
        """
        
        if self.sequence_id < other.sequence_id:
            return True
        if self.sequence_id > other.sequence_id:
            return False
        
        # IDs must be equal    
        
        if self.offset < other.offset:
            return True
        if self.offset > other.offset:
            return False
            
        # Offsets must be equal
            
        if self.is_right < other.is_right:
            return True
        if self.is_right > other.is_right:
            return False
            
        # Must be equal
        return False
        
class Join(object):
    """
    Represents an adjacency between two Sides.
    """
    
    def __init__(self, side1, side2):
        """
        Make a new Join between the given Sides.
        
        """
        
        # Make a set of the sides.
        self.sides = set([side1, side2])
        
class SideGraph(object):
    """
    Represents a side graph composed of sequences and joins.
    """
    
    def __init__(self):
        """
        Make a new empty side graph.
        """
        
        # Hold sequences by ID
        self.sequences = {}
        
        # Hold joins indexed by each of their sequences
        self.joins_by_sequence = collections.defaultdict(list)
    
    def add_sequence(self, sequence_id, length, sequence):
        """
        Add a new sequence to the graph. Takes the ID, the length, and the
        bases (as a BioPython Seq object).
        
        """
        
        # Make and add the sequence
        self.sequences[sequence_id] = Sequence(sequence_id, length)
        
    def add_join(self, side1, side2):
        """
        Add a join between the two given sides.
        
        """
        
        # Make the join
        join = Join(side1, side2)
        
        print("Joining {} to {}".format(side1, side2))
        
        # Insert it in the list under the sequence IDs
        self.joins_by_sequence[side1.sequence_id].append(join)
        self.joins_by_sequence[side2.sequence_id].append(join)
        
class DatabaseSideGraph(object):
    """
    A side graph backed by a database.
    
    TODO: just pull in the GA4GH server code?
    """
    
    def __init__(self, database_filename, fasta_filename):
        """
        Make a new side graph, writing to the given sqlite database and FASTA
        file.
        
        TODO: reading.
        """
        
        # Connect to the database
        self.database = sqlite3.connect(database_filename)
        
        # Open the FASTA file for writing
        self.fasta = open(fasta_filename, "w")
        
        # Set up the database. Copied from Maciek's schema
        self.database.executescript("""
        CREATE TABLE FASTA (ID INTEGER PRIMARY KEY,
	        fastaURI TEXT NOT NULL);
        --
        --
        -- Sequences and joins - basis of graph topology
        --
        CREATE TABLE Sequence (ID INTEGER PRIMARY KEY,
	        fastaID INTEGER NOT NULL REFERENCES FASTA(ID), -- the FASTA file that contains this sequence's bases.
	        sequenceRecordName TEXT NOT NULL, -- access to the sequence bases in the FASTA file ONLY.
	        md5checksum TEXT NOT NULL, -- checksum of the base sequence as found in the FASTA record.
	        length INTEGER NOT NULL); -- length of the base sequence as found in the FASTA record.
        --
        --
        CREATE TABLE GraphJoin (ID INTEGER PRIMARY KEY,
	        --
	        -- by convention, side1 < side2 in the lexicographic ordering defined by (sequenceID, position, forward).
	        --
	        side1SequenceID INTEGER NOT NULL REFERENCES Sequence(ID),
	        side1Position INTEGER NOT NULL, -- 0 based indexing, counting from 5' end of sequence.
	        side1StrandIsForward BOOLEAN NOT NULL, -- true if this side joins to 5' end of the base
	        -- 
	        side2SequenceID INTEGER NOT NULL REFERENCES Sequence(ID),
	        side2Position INTEGER NOT NULL,
	        side2StrandIsForward BOOLEAN NOT NULL);
        """)
        
        # Add a FASTA file entry for our FASTA
        self.database.execute("INSERT INTO FASTA VALUES (0, ?);", (
            fasta_filename,))
        
    
    def add_sequence(self, sequence_id, length, sequence):
        """
        Add a new sequence to the graph. Takes the ID, the length, and the
        bases (as a BioPython Seq object).
        
        Sequence ID must be an int.
        
        """
        
        # Slap an ID on the sequence and write it to the FASTA
        Bio.SeqIO.write(Bio.SeqRecord.SeqRecord(sequence, str(sequence_id)),
            self.fasta, "fasta")
        
        # Hash the sequence
        hasher = hashlib.md5()
        hasher.update(str(sequence))
        md5sum = hasher.hexdigest()
        
        
        # Put it in the database, pointing to FASTA 0
        self.database.execute("INSERT INTO Sequence VALUES (?, 0, ?, ?, ?);",
            (sequence_id, str(sequence_id), md5sum, len(sequence))) 
        
        
    def add_join(self, side1, side2):
        """
        Add a join between the two given sides.
        
        """
        
        if side1 > side2:
            # Flip around
            side1, side2 = side2, side1
        
        # Stick it in the database. Don't put a primary key; we will get one for
        # free.
        self.database.execute("INSERT INTO GraphJoin(side1SequenceID, "
            "side1Position, side1StrandIsForward, side2SequenceID, "
            "side2Position, side2StrandIsForward) VALUES (?, ?, ?, ?, ?, ?);",
            (side1.sequence_id, side1.offset, not side1.is_right,
            side2.sequence_id, side2.offset, not side2.is_right))
            
        
def get_max_overlap(overlaps, sides_by_endpoint):
    """
    Get the max overlap form any of the (endpoint, overlap) pairs that
    corresponds to an existing sequence. Returns (endpoint, overlap). If there
    is no such overlap, returns (None, 0).
    """
    
    # The max is of course no overlap to start
    best = (None, 0)
    
    for endpoint, overlap in overlaps:
        # For every overlap instance, we have to see if it's to a thing we have
        # already.
        if sides_by_endpoint.has_key(endpoint):
            # This is to an existing sequence.
            
            if best[0] is None or best[1] < overlap:
                # We have a new best pair
                best = (endpoint, overlap)
            
    return best


def process_overlaps(overlaps, longest_overlap, novel_sequence_side,
    sides_by_endpoint, graph):
    """
    Given a set of overlap (endpoint, length) tuples, the tuple out of these
    that is the longest while still pointing to an endpoint registered in the
    endpoint-to-sides dict, the Side of the newly added sequence that joins
    directly to that maximally-overlapping endpoint, the endpoints-to-sides
    dict, and the graph we are building, make all the joins in the graph to
    express the overlaps to endpoints that have been created already.
    
    """
    
    for endpoint, overlap in overlaps:
        # For all the overlaps...
        # If there are none on this end we'll just skip the whole loop
    
        if not sides_by_endpoint.has_key(endpoint):
            # Skip all the overlaps with things that don't exist yet.
            continue 
            
        # How many bases do we have to go into the sequence we are joining
        # onto by the longest overlap?
        overlap_offset = longest_overlap[1] - overlap
        
        # Get the side we should be welding to at the distance into the newly
        # added record for this overlap.
        stand_in_side = get_merged_side(longest_overlap, novel_sequence_side,
            overlap, sides_by_endpoint)
        
        # Attach to it
        graph.add_join(sides_by_endpoint[endpoint], stand_in_side)
        
    # Now we have added joins to the end of the novel sequence for all the
    # longest overlaps, and to the designated merged-into sequence for shorter
    # overlaps.
    
def get_merged_side(longest_overlap, novel_sequence_side, in_from_endpoint,
    sides_by_endpoint):
    """
    Given the longest overlap (endpoint, overlap) tuple for the appropriate end
    of an added sequence, the end Side of the added novel sequence for that end,
    the distance in from the endpoint of the record that produced the novel
    sequence, and the endpoint-to-side dict, get the outward-faceing Side for
    the given distance in from the *endpoint* of the record that produced the
    newly added sequence (even if that Side is on a base that has been merged in
    elsewhere due to the longer overlap).
    
    """
    
    if longest_overlap[1] == 0:
        # No bases are overlapped
        
        # We're counting from the end of the novel sequence, and we need to wind
        # it inwards.
        return novel_sequence_side.add_local_offset(-in_from_endpoint)
    elif longest_overlap[1] <= in_from_endpoint:
        # We need to be in the newly added sequence, because we are further in
        # than the longest overlap.
        
        # We need to walk in the difference between the longest overlap and the
        # offset we wanted.
        return novel_sequence_side.add_local_offset(longest_overlap[1] -
            in_from_endpoint)
    else:
        # We want a distance in from the endpoint less than the longest overlap,
        # so the answer is going to be in the sequence that provided the longest
        # overlap.
    
        # How many bases do we have to go into the sequence we are joining
        # onto?
        overlap_offset = longest_overlap[1] - in_from_endpoint
        
        # We need to go one or more bases into the sequence we
        # overlapped onto.
        
        # This Side stands in for the Side that we wanted to attach to, which is
        # off the end of the new Sequence. First we flip the Side we officially
        # merged to, and then (having counted switching to it as a base) we walk
        # in the remainder of the offset.
        return sides_by_endpoint[longest_overlap[0]].flip().add_local_offset(
            overlap_offset - 1)
    
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
    
    # We want the average unitig length, so we need to track these totals.
    total_length = 0
    total_unitigs = 0
    
    # We need to build a side graph in a database.
    graph = DatabaseSideGraph(options.database, options.fasta)
    
    # Keep track of MAG endpoint numbers. This maps endpoint numbers to where
    # they fall in the graph. It's not necessary to keep more than one, because
    # when you come along overlapping an endpoint, we just need a Join from that
    # endpoint to the part of your sequence not doing any overlap.
    sides_by_endpoint = {}
    
    # We want to know how fast we read stuff
    last_time = time.clock()
    # How many records per reporting period?
    interval = 10000
    
    # Grab a handle to this so we can track bytes
    expander = TransparentUnzip(options.in_file)
    expander.start_stats()
    
    for ends, _, left_overlaps, right_overlaps, seq in parse_mag(expander):
        # Unpack each unitig
        
        # Add each unitig to the totals
        total_unitigs += 1
        total_length += len(seq)
        
        # Find the longest overlap on each side with existing sequence
        max_left_overlap = get_max_overlap(left_overlaps, sides_by_endpoint)
        max_right_overlap = get_max_overlap(right_overlaps, sides_by_endpoint)
        
        # Add in the novel sequence
        sequence_id = total_unitigs
        sequence_length = len(seq) - max_left_overlap[1] - max_right_overlap[1]
        if sequence_length == 0:
            # TODO: Figure out how to handle this
            raise Exception("Sequences can't be all overlap.")
        
        # There is any novel sequence here. Make a sequence with that sequence
        # data.
        graph.add_sequence(sequence_id, sequence_length, 
            seq[max_left_overlap[1]:len(seq) - max_right_overlap[1]])
        
        # Define the sequence left and right Sides
        sequence_left_side = Side(sequence_id, 0, False)
        sequence_right_side = Side(sequence_id, sequence_length - 1, True)
            
        # Do the joins for all the left overlaps
        process_overlaps(left_overlaps, max_left_overlap, sequence_left_side,
            sides_by_endpoint, graph)
        # And the right ones
        process_overlaps(right_overlaps, max_right_overlap, sequence_right_side,
            sides_by_endpoint, graph)
                       
        # Add the Sides for the endpoints of this sequence. They probably got
        # merged in, if we had an overlap on either side. Just look 0 in from 
        sides_by_endpoint[ends[0]] = get_merged_side(max_left_overlap,
            sequence_left_side, 0, sides_by_endpoint)
        sides_by_endpoint[ends[1]] = get_merged_side(max_right_overlap,
            sequence_right_side, 0, sides_by_endpoint)
        
        
        if total_unitigs % interval == 0:
            # We should print status
            
            # What time is it? How long have we been working?
            now_time = time.clock()
            elapsed = now_time - last_time
            last_time = now_time
            
            # How many bytes did we process
            compressed_bytes, uncompressed_bytes = expander.get_stats()
            expander.start_stats()
        
            print("Processed {} unitigs, {:.2f} per second, "
                "{:.2f} kbps/{:.2f} kbps".format(
                total_unitigs, interval / elapsed,
                compressed_bytes / elapsed / 1000,
                uncompressed_bytes / elapsed / 1000))
        
    print("Average unitig length: {}".format(total_length / 
        float(total_unitigs)))
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

