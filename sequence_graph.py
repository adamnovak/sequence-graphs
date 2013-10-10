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
import pulp

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

class Entity(object):
    """
    Represents an Entity in a sequence graph. It has a series of Components that
    determine what it is actually like. It also has a unique ID.
    
    An AlleleGroup has a SequenceComponent, two EndComponents, and a
    PloidyComponent.
    
    An Adjacency has a PloidyComponent and two EndComponents.
    
    Telomeres and Contig Ends both have a PloidyComponent and one EndComponent.
    """
    
    def __init__(self):
        """
        Make a new Entity.
        
        """
        
        # Keep our ID integer
        self.id = get_id()
        
        # Make a place to store our components. It needs to be in an order so we
        # can have a 5' and 3' end if necessary.
        self.components = []
        
    def add_component(self, component):
        """
        Add the given Component to this Entity.
        
        """
        
        # Put the component in our list
        self.components.append(component)
        
        # Make the component know it belongs to us
        component.entity = self
        
    def __str__(self):
        """
        Represent this Entity as a human-readable string.
        
        """
        
        return "<Entity #{}>".format(self.id)
        
class Component(object):
    """
    Represents a component that can be part of an Entity. Knows which Entity it
    has been attached to, after it is attached.
    
    Must be attached to an Entity with Entity.add_component()
    
    """
    
    def __init__(self):
        """
        Make a new Component not attached to any Entity.
        
        """
        
        # When we are attached to an entity, this will hold that entity. For now
        # it holds None.
        self.entity = None
        
class PloidyComponent(Component):
    """
    Allows an Entity to have a ploidy value (represented by an LP variable).
    
    """
    
    def __init__(self):
        """
        Make a new PloidyComponent attached to no entity.
        
        """
        
        # Initialize the base class
        super(PloidyComponent, self).__init__()
        
        # Make the LP variable
        self.variable = pulp.LpVariable("Ploidy_".format(get_id()), 
            cat="Integer")
            
    def get_ploidy(self):
        """
        Return the ploidy as an integer. The relavent
        Linear Programing problem must be solved first.
        
        """
        
        # Go get the value from pulp
        return pulp.value(self.variable)

class EndComponent(Component):
    """
    Allows an Entity to connect to another Entity. Contains a reference to one
    or more other EndComponents that it is connected to, and to the "other"
    EndComponent of the entity it belongs to, which may or may not actually be
    this EndComponent.
    
    EndComponents on AlleleGroups, Contig Ends, and Telomeres should connect to
    EndComponents on one or more Adjacencies. EndComponents on Adjacencies
    should connect to exactly one EndComponent on an AlleleGroup, Contig End, or
    Telomere.
    
    """
    
    def __init__(self):
        """
        Make a new EndComponent attached to no entity.
        
        """
        # Initialize the base class
        super(EndComponent, self).__init__()
        
        # Make a set of connected EndComponents
        self.connections = set()
        
        # Make a spot to put the other end of the thing we are attached to. For
        # now assume it's us.
        self.other_end = self
        
        
        
    def connect(self, other):
        """
        Connect this EndComponent to another EndComponent. Only needs to be
        called one way.
        
        """
        
        # Save our reference to them
        self.connections.add(other)
        # Give them a reference to us
        other.connections.add(self)
                
        
        
class AlleleGroup(Entity):
    """
    Represents an AlleleGroup. An Entity with a Ploidy and two Ends.
    
    """
    
    def __init__(self):
        """
        Make a new AlleleGroup. Add all the components.
        
        """
        
        # Initialize the base class
        super(AlleleGroup, self).__init__()
        
        # Add the ploidy component
        self.add_component(PloidyComponent())
        
        # TODO: We're not supposed to really have any logic here. We should have
        # just one Graph Member component or something to represent a two-sided
        # graph node.
        
        # Make the two EndComponents
        end_1 = EndComponent()
        end_2 = end_component()
        
        # Point them at each other
        end_1.other_end = end_2
        end_2.other_end = end_1
        
        # Add the end components
        self.add_component(end_1)
        self.add_component(end_2)
        
class Adjacency(Entity):
    """
    Represents an Adjacency. An Entity with a Ploidy and two Ends.
    
    """
    
    def __init__(self):
        """
        Make a new AlleleGroup. Add all the components.
        
        """
        
        # Initialize the base class
        super(Adjacency, self).__init__()
        
        # Add the ploidy component
        self.add_component(PloidyComponent())
        
        # Make the two EndComponents
        end_1 = EndComponent()
        end_2 = EndComponent()
        
        # Point them at each other
        end_1.other_end = end_2
        end_2.other_end = end_1
        
        # Add the end components
        self.add_component(end_1)
        self.add_component(end_2)
    
            
    
        
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
