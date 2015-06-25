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
import networkx


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
    parser.add_argument("in_file", type=argparse.FileType("r"),
        default=sys.stdin,
        help="MAG-format file of the unitig graph to read")
    
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
    
    # Let's build an overlap graph between ends
    graph = networkx.Graph()
    
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
        
        # Add the sequence edge to the graph, annotating it with properties
        graph.add_edge(ends[0], ends[1], {"type": "sequence", 
            "length": len(seq)})
            
        # Add the overlap edges
        for destination, length in left_overlaps:
            # Each of the left overlaps connects to the left end.
            # This edge may already exist, from the other side.
            graph.add_edge(destination, ends[0], {"type": "overlap",
                "length": length})
                
        for destination, length in right_overlaps:
            # Each of the right overlaps connects to the right end.
            # This edge may already exist, from the other side.
            graph.add_edge(destination, ends[1], {"type": "overlap",
                "length": length})
        
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
        
    # Count all the overlaps
    total_overlaps = 0
    for (node1, node2, data) in graph.edges_iter():
        # Check each edge
        
        if data["type"] == "overlap":
            # And count the overlaps
            total_overlaps += 1
            
    print("{} overlaps among {} sequences".format(total_overlaps,
        total_unitigs))
        
    # Now do the connected components
    print("Connected components: {}".format(
        networkx.number_connected_components(graph)))
        
    # Throw out the sequences and ask that question again.
    # First we need to make a list (so we aren't iterating while removing)
    bad_edges = [(u, v) for (u, v, data) in
        graph.edges_iter() if data["type"] == "sequence"]
        
    graph.remove_edges_from(bad_edges)
    
    print("Overlap-connected components: {}".format(
        networkx.number_connected_components(graph)))
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

