#!/usr/bin/env python
"""
sequence_graph.py: An implementation of David Haussler's Sequence Graph data
model for representing potentially incompletely specified genomes.

A Genome is composed of AlleleGroups (which contain an Allele and a ploidy)
connected by Adjacencies (which contain a ploidy). Each AllleleGroup has two
distinct ends (5' and 3'). An AlleleGroup's Allele may be reference (matching
the reference genome sequence) or variant. Each AlleleGroup belongs to a single
Site, and each Site contains zero or more allele groups. Sites can be reference
(and have coordinates in the reference genome) or non-reference (and represent
novel insertions). Each adjacency belongs to exactly one Breakpoint, and a
Breakpoint contains zero or more adjacencies. All Adjacencies in a single
Breakpoint must connect the same ends of AlleleGroups in the same pair of Sites.
Breakpoints that connect two Reference sites in the reference order and
orientation are reference Breakpoints, whilke other breakpoints are non-
reference.

Contig Ends and Telomeres are special single-ended AlleleGroup-type things that
occupy their own Sites.

A Genome is incompletely specified if it contains any allele groups or
breakpoints with a ploidy that is unknown or greater than one. A Genome is fully
specified if all of the ploidies of its AlleleGroups and Breakpoints are either
0 or 1.

Zero or more Genomes can belong to a Model, which allows a set of Integer Linear
Programming constriants to be solved to generate sets of fully-specified
versions of each Genome that are mutually consistent.

"""

import sys, logging, math
import pulp, networkx

class PenaltyTree(object):
    """
    Maintains a tree of penalty terms, so that we can have arbitrarily many
    penalty terms without arbitrarily large constraint expressions.
    
    """
    
    def __init__(self, degree=100):
        """
        Make a new PenaltyTree. degree specifies the maximum number of terms to
        sum together at once. Only one PenaltyTree may be used on a given LP
        problem.
        
        """
        
        # This holds the number of children per node/number of terms to sum at
        # once.
        self.degree = degree
        
        # This holds all our leaf-level terms.
        self.terms = []
        
    def get_variable(self):
        """
        Return a fresh LpVariable with a unique name.
        
        """
        
        # Make the variable
        var = pulp.LpVariable("PenaltyTree_{}".format(get_id()))
        
        # Give the fresh variable to the caller.
        return var
        
    def add_term(self, term):
        """
        Add the given LP expression as a term in the tree.
        
        """
        
        self.terms.append(term)
        
    def set_objective(self, problem):
        """
        Add the sum of all terms as the given LP problem's objective. The
        PenaltyTree must have at least one term.
        
        """
        
        # Algorithm: Go through our leaves, making a variable for the sum of
        # each group of self.degree terms. Then repeat on the sums, making sums
        # of sums, and so on. When we only have one sum, make that the
        # objective.
        
        # This holds the list we're collecting
        collecting = self.terms
        
        # This holds the list of sum variables
        sums = []
        
        while len(collecting) > 1:
            logging.info("Collecting {} terms in groups of {}".format(len(
                collecting), self.degree))
            
            for i in xrange(0, len(collecting), self.degree):
                # This holds the terms we have collected for this sum
                collected = []
                for j in xrange(0, min(self.degree, len(collecting) - i)):
                    # Grab each term up to our degree that actually exists
                    collected.append(collecting[i + j])
                    
                # This holds the variable we use to represent the sum of the
                # things we have collected
                sum_var = self.get_variable()

                # Constrain this variable to equal the sum of all the collected
                # terms
                problem += sum(collected) == sum_var

                # Add this variable for the next level of collection
                sums.append(sum_var)
               
            # Move up a level in the tree
            collecting = sums
            sums = []
           
        # We have now collected everything down to one term, which is in the
        # collecting list. Use it as our objective function.
        problem += collecting[0]
        
class SequenceGraphLpProblem(object):
    """
    Represents an LP copy number problem. You can attach several models to them,
    constrain them together, and solve the problem.
    
    Internally, contains a pulp LpProblem, and a PenaltyTree.
    
    """
    
    def __init__(self):
        """
        Make a new SequenceGraphLpProblem that we can solve.
        
        """
        
        # We need an actual LpProblem
        self.problem = pulp.LpProblem("copynumber", pulp.LpMinimize)
        
        # We also need a PenaltyTree for organizing penalty terms
        self.penalties = PenaltyTree()
        
    def constrain_approximately_equal(self, var_a, var_b, penalty=1):
        """
        Constrain the two LP variables (or constants) var_a and var_b to be
        approximately equal, subject to the given penalty.
        
        Adds the appropriate constraints to the CopyNumberLpProblem's internal
        LpProblem, and penalties to the model's PenaltyTree.
        
        """
        
        # Make an LP variable for the amount that var_b is above var_a. Note
        # that str() on variables produces their names. Also, we have to make
        # sure that this starts with a letter.
        amount_over = pulp.LpVariable("over_{}".format(get_id()), 0)
            
        # Add the constraint for not being more than that much over
        self.add_constraint(var_b <= var_a + amount_over)
        
        # Make an LP variable for the amount that var_b is below var_a
        amount_under = pulp.LpVariable("under_{}".format(get_id()), 0)
            
        # Add the constraint for not being more than that much under
        self.add_constraint(var_b >= var_a - amount_under)
        
        # Apply an equal penalty in each direction
        self.add_penalty((penalty * amount_over) + (penalty * amount_under))
            
    def add_penalty(self, term):
        """
        Add the given penalty term to the problem's objective.
        
        """
        
        # Just put the term in the PenaltyTree
        self.penalties.add_term(term)
        
    def add_constraint(self, constraint):
        """
        Add the given (exact) constraint to the problem. For approximate
        constraints, use constrain_approximately_equal() instead.
        
        """
        
        # Just add the constraint to the internal problem.
        self.problem += constraint
        
    def solve(self, save=None):
        """
        Solve the LP problem with the best solver we can find. After solving,
        you can get_ploidy() on all the AlleleGroups in Models attached to this
        problem.
        
        If save is specified, it is a filename to which to save the LP problem
        in LP format.
        
        You may only solve a SequenceGraphLpProblem once.
        
        """
        
        # Set up the penalties described by the penalty tree
        logging.info("Setting up penalties")
        self.penalties.set_objective(self.problem)
        
        if save is not None:
            logging.info("Saving problem to {}".format(save))
            self.problem.writeLP(save)
        
        logging.info("Looking for solver")
        # This is a list of solvers to use in reverse order of priority (last is
        # best)
        candidates = [
            pulp.GLPK(),
            # We're not actually using 132 threads here, just 32; the 100 means
            # use deterministic parallelism, working around the segfault issue
            # described here: <http://list.coin-
            # or.org/pipermail/cbc/2013-March/001044.html>
            pulp.solvers.COIN_CMD(threads=132, msg=1)
        ]
        
        # This is the solver we actually use
        solver = None
        for candidate in candidates:
            if candidate.available():
                logging.info("{} is available".format(
                    candidate.__class__.__name__))
                solver = candidate
            else:
                logging.info("{} is unavailable".format(
                    candidate.__class__.__name__))
        
        if solver is None:
            logging.critical("No solver found!")
            raise Exception("No LP solver available")
        
        logging.info("Solving with {}...".format(solver.__class__.__name__))
        
        # Solve the problem
        status = self.problem.solve(solver)
        logging.info("Solution status: {}".format(pulp.LpStatus[status]))
        
        if len(self.problem.variables()) < 20:
            # It's short enough to look at.
            for var in self.problem.variables():
                logging.debug("\t{} = {}".format(var.name, pulp.value(var)))
                
        # Report the total penalty we got when we solved.
        logging.info("Penalty: {}".format(pulp.value(self.problem.objective)))
        
        if status != pulp.constants.LpStatusOptimal:
            raise Exception("Unable to solve problem optimally.")
            
def get_id():
    """
    Return a unique integer ID (for this execution). Use this to ensure your LP
    variable names are unique.
    
    """
    
    if not hasattr(get_id, "last_id"):
        # Start our IDs at 0, which means the previous ID was -1. Static
        # variables in Python are hard.
        setattr(get_id, "last_id", -1)
        
    # Advance the ID and return the fresh one.
    get_id.last_id += 1
    return get_id.last_id

class Genome(object):
    """
    Represents an (incompletley specified) genome as a sequence graph. Wraps all
    access to the internal graph data structures, which the user never gets to
    hear about.
    """
    
    def __init__(self):
        """
        Make a new Genome sequence graph.
        
        """
        
        # To actually implement this API, including the ability to look up, say,
        # Adjacencies by both Breakpoint they belong to and their own ID, we're
        # going to need a proper database.
        
    def add_site(self, reference_contig=None, reference_start=None, 
        reference_end=None):
        """
        Add a new Site to the sequence graph, with an optional reference
        position.
        
        Returns the ID of the Site added.
        """
        
    def add_breakpoint(self, site_a, side_a, site_b, side_b):
        """
        Add a new Breakpoint between the specified site IDs, on the specified
        side (0 or 1) of each Site. Note that flipping the A and B sites (and
        their sides) still indicates the same Breakpoint. Note that you can have
        a Breakpoitn between the two sides of a single Site, or attached to only
        one side of one Site (in which case it is a self-loop in the underlying
        graph).
        
        Returns the ID of the Breakpoint added. 
        """
        
    def add_allele_group(self, site, sequence=None):
        """
        Add a new AlleleGroup to the Site with the given ID. If specified,
        sequence is a DNA sequence to associate with the allele group in the
        side 0 to side 1 direction.
        
        Returns the ID of the AlleleGroup added.
        
        """
        
    def add_adjacency(self, breakpoint, allele_group_a, allege_group_b):
        """
        Add a new Adjacency to the Breakpoint with the given ID, connecting the
        given pair of Allele Groups (also specified by ID). The order of the
        AlleleGroups does not matter. The ends that the Adjacency connects are
        specified by the Breakpoint to which it is added.
        
        Returns the ID of the Adjacency added.
        
        """
        
    def get_sites(self):
        """
        Iterate over Site IDs in an arbitrary order.
        
        """
    
    def get_allele_groups(self, site):
        """
        Iterate over AlleleGroup IDs in a Site specitfied by ID in arbitrary
        order.
        
        """
        
    def get_breakpoints(self):
        """
        Iterate over Breakpoint IDs in an arbitrary order.
        
        """
        
    def get_adjacencies(self, breakpoint):
        """
        Iterate over Adjacency IDs in a Breakpoint specified by ID in arbitrary
        order.
        
        """
        
    def set_constraint(self, to_constrain, min_ploidy, max_ploidy):
        """
        Set the minimum and maximum integer ploidy on an AlleleGroup or
        Adjacency.
        
        """
        
    def solutions():
        """
        Iterate over solutions to the ploidy constraints that have been
        specified. Solutions are Solution objects with a get_ploidy method that
        returns the ploidies of AlleleGroups and Adjacencies by ID.
        
        """
        
class Solution(object):
    """
    Represents a complete ploidy specification for an incompletely specified
    genome. Just a collection of ploidies by AlleleGroup or Adjacency ID.
    
    """
    
    def __init__(self):
        """
        Make a new Solution.
        
        """
        
    def get_ploidy(self, to_get):
        """
        Given the ID of an AlleleGroup or Adjacency, get the ploidy of that
        AlleleGroup or Adjacency under this Solution.
        
        """
        
# Test code
if __name__ == "__main__":
    # Test out the API
    
    # Make a Genome
    genome = Genome()
    
    # Add some Sites
    site_a = genome.add_site("chr1", 0, 9)
    site_b = genome.add_site()
    site_c = genome.add_site("chr1", 10, 19)
    
    # Add some Breakpoints
    break_ab = genome.add_breakpoint(site_a, 1, site_b, 0)
    break_bc = genome.add_breakpoint(site_b, 1, site_c, 0)
    break_ac = genome.add_breakpoint(site_a, 1, site_c, 0)
    
    # Add some AlleleGroups
    group_a = genome.add_allele_group(site_a, "AACTGGCACC")
    group_b = genome.add_allele_group(site_b, "CCAG")
    group_c = genome.add_allele_group(site_c, "CTTACGGATG")
    
    # Add some Adjacencies
    adj_ab = genome.add_adjacency(break_ab, group_a, group_b)
    adj_bc = genome.add_adjacency(break_bc, group_b, group_c)
    adj_ac = genome.add_adjacency(break_ac, group_a, group_c)
    
    # Add some ploidy constraints
    genome.set_constraint(adj_ab, 1, 2)
    genome.set_constraint(adj_bc, 1, 2)
    genome.set_constriant(adj_ac, 1, 2)
    
    for solution in genome.solutions():
        # For each solution
        for site in genome.get_sites():
            # For each site
            print "Site {}".format(site)
            for allele_group in genome.get_allele_groups(site):
                # Print the solution's ploidy for each AlleleGroup
                print "AlleleGroup {}: {}".format(allele_group, 
                    solution.get_ploidy(site))
