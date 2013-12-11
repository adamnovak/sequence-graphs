#!/usr/bin/env python2.7
"""
visualize.py: draw a GraphViz graph of an Avro-format variant sequence graph.

This script reads two Avro-format files: a file of AlleleGroups and a file of
Adjacencies. It produces a visualization of the graph contained in the files,
using pygraphviz.

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, collections, math, itertools, logging

import avro.schema, avro.datafile, avro.io
import pygraphviz

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
    parser.add_argument("allele_group_file", type=argparse.FileType("r"),
        help="input file for allele groups")
    parser.add_argument("adjacency_file", type=argparse.FileType("r"),
        help="input file for adjacencies")
    parser.add_argument("graph_file",
        help="file to save graph picture in")
    parser.add_argument("--samples", nargs="+", type=set, default=None,
        help="list of samples to restrict to")
        
    # Logging options
    parser.add_argument("--loglevel", default="DEBUG", choices=["DEBUG", "INFO",
        "WARNING", "ERROR", "CRITICAL"],
        help="logging level to use")
        
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def set_loglevel(loglevel):
    """
    Given a string log level name, set the logging module's log level to that
    level. Raise an exception if the log level is invalid.
    
    """
    
    # Borrows heavily from (actually is a straight copy of) the recipe at
    # <http://docs.python.org/2/howto/logging.html>
    
    # What's the log level number from the module (or None if not found)?
    numeric_level = getattr(logging, loglevel.upper(), None)
    
    if not isinstance(numeric_level, int):
        # Turns out not to be a valid log level
        raise ValueError("Invalid log level: {}".format(loglevel))
    
    # Set the log level to the correct numeric value
    logging.basicConfig(format="%(levelname)s:%(message)s", level=numeric_level)

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object

    # Set up logging
    set_loglevel(options.loglevel)
        
    # Make Avro-format input file readers. This one is for allele groups.
    allele_group_reader = avro.datafile.DataFileReader(
        options.allele_group_file, avro.io.DatumReader())
    # This one is for adjacencies
    adjacency_reader = avro.datafile.DataFileReader(options.adjacency_file,
        avro.io.DatumReader())
    
    # Make a non-strict graph (i.e. a multigraph) to populate
    graph = pygraphviz.AGraph(strict=False)
    
    for allele_group in allele_group_reader:
        if (options.samples is not None and 
            allele_group["sample"] not in options.samples):
            
            # We want to skip this, it's not in the right sample
            continue
    
        # Make a node for each end
        graph.add_node(allele_group["fivePrime"], label="5'", shape="point")
        graph.add_node(allele_group["threePrime"], label="3'", shape="point")
        
        # Work out what to label the sequence node
        label = "{}:{:,}-{:,}".format(allele_group["contig"], 
            allele_group["start"], allele_group["end"])
        
        if allele_group["sequence"] is not None:
            # This is a non-reference site (or an insert). Use the DNA as the
            # label.
            label = "{}\n{}".format(label, allele_group["sequence"])
        elif allele_group["end"] == allele_group["start"]:
            # It's 0-length, so it's the lack of an insert. Don't say anything.
            label = ""
        else:
            # It's a piece of reference DNA. Just say the coordinates.
            pass
            
        
        # Make a node to represent the sequence itself.
        graph.add_node(allele_group["id"], label=label, shape="box")
        
        # Make edges connecting the sequence to its ends, with the AlleleGroup
        # ID as the key.
        graph.add_edge(allele_group["fivePrime"], allele_group["id"], 
            key=allele_group["id"], color="blue")
        graph.add_edge(allele_group["id"], allele_group["threePrime"], 
            key=allele_group["id"], color="blue")
        
    for adjacency in adjacency_reader:
        if (options.samples is not None and 
            adjacency["sample"] not in options.samples):
            
            # We want to skip this, it's not in the right sample
            continue
    
        # Make an edge for this adjacency, with a unique key.
        graph.add_edge(adjacency["first"], adjacency["second"], 
            key=adjacency["id"], color="red")
        
            
    # Close up the files
    allele_group_reader.close()
    adjacency_reader.close()
    
    # Lay out the graph
    graph.layout()
    
    # Draw the graph
    graph.draw(options.graph_file)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
