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

A Genome is incompletely specified if it contains any allele groups or
breakpoints with a ploidy that is unknown or greater than one. A Genome is fully
specified if all of the ploidies of its AlleleGroups and Breakpoints are either
0 or 1.

Zero or more Genomes can belong to a Model, which allows a set of Integer Linear
Programming constriants to be solved to generate sets of fully-specified
versions of each Genome that are mutually consistent.

"""

import sys, logging
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
        

class MeasurementLocation(object):
    """
    Represents a location where copy number has been (ambiguously) measured.
    Keeps track of a string measurement identifier, the mapping location in the
    genome, and the relative contribution of this location to the measurement.
    
    """
    
    def __init__(self, measurement, contig, start, end, affinity):
        """
        Make a new MeasurementLocation, representing the mapping of the given
        probe to the region between the given start and end positions, with the
        given relative affinity.
        
        If affinity is a callable, it will be called once if the affinity of
        this MeasurementLocation is requested. Otherwise, it will
        be used as the actual affinity value.
        
        """
        
        # Save the measurement name
        self.measurement = measurement
        # And the contig
        self.contig = contig
        # And the start position
        self.start = start
        # And the end position
        self.end = end
        
        if callable(affinity):
            # We've been given a function to lazily call when someone asks for
            # our affinity. Remember it.
            self.affinity_function = affinity
            
            # Make a place to cache the affinity value we get when we call the
            # above function.
            self.affinity = None
            
        else:
            # Just use the given affinity as our actual affinity value.
            self.affinity = affinity
            
    def get_affinity(self):
        """
        Return the affinity of this MeasurementLocation: the weight with which
        the AlleleGroup it is in contributes to the total measurement.
        
        If there is no affinity stored, call the affinity function and cache its
        result. Otherwise, return the stored affinity.
        
        """
        
        if self.affinity is None:
            # Go calculate the affinity. This is expensive and probably goes and
            # calls out to R
            logging.debug("Calculating affinity for {} at {}:{}-{}".format(
                self.measurement, self.contig, self.start, self.end))
            self.affinity = self.affinity_function()
            
        return self.affinity
        
class AlleleGroup(object):
    """
    Represents an AlleleGroup. Has a ploidy, and maintains a list of
    MeasurementLocations inside it.
    
    """
    
    def __init__(self, contig, start, end, min_ploidy=0):
        """
        Make a new AlleleGroup, with an integer LP variable to represent its
        ploidy. 
        
        contig, start, and end specify this AlleleGroup's position in the linear
        reference; they must all be None for contigs that are not present in the
        reference.
        
        min_ploidy, if specified, gives a minimum ploidy for the allele
        group.
        
        """
        
        # Save our reference chromosome
        self.contig = contig
        # And reference start position (0-based)
        self.start = start
        # And reference end position (1-past-the-end)
        self.end = end
        
        # This holds the ploidy LP variable
        self.ploidy = pulp.LpVariable("AG_{}".format(get_id()), 
            min_ploidy, cat="Integer")
            
        # Keep a list of measurements that measure us
        self.measurements = []
            
            
            
    def get_ploidy(self):
        """
        Return the ploidy of this allele group as an integer. The relavent
        Linear Programing problem must be solved first.
        
        """
        
        # Go get the value from pulp
        return pulp.value(self.ploidy)
        
    def add_measurement(self, measurement_location):
        """
        Add a MeasurementLocation to this AlleleGroup's measurements list.
        
        The MeasurementLocation should be contained within the AlleleGroup in
        the genome.
        
        """
        
        # Add it to the list
        self.measurements.append(measurement_location)
        
class Model(object):
    """
    Represents a Sequence Graph model of the copy number of a genome. It
    contains lists of AlleleGroups with reference genome locations, separated
    out by contig. New AlleleGroups can be inserted, and, when all AlleleGroups
    have been made, approximate-equality constraints can be generated tying
    AlleleGroup copy numbers to those of their neighbors.
    
    All constraints and penalty terms are automatically added to the
    SequenceGraphLpProblem to which the Model is attached (specified on
    construction).
    
    """
    
    def __init__(self, problem):
        """
        Make a new Model attached to the given SequenceGraphLpProblem. Once the
        problem is solved, the Model will have copy number values available from
        its AlleleGroups.
        
        """
        
        # We need a place to keep our AlleleGroup lists
        self.allele_groups = collections.defaultdict(list)
        
        # We need to remember the problem we're attached to
        self.problem = problem
        
    def constrain_approximately_equal(self, var_a, var_b, penalty=1):
        """
        Constrain the two LP variables (or constants) var_a and var_b to be
        approximately equal, subject to the given penalty.
        
        Adds the appropriate constraints and penalties to the Model's
        SequenceGraphLpProblem.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        
        # Just forward the call on to the problem.
        self.problem.constrain_approximately_equal(var_a, var_b, 
            penalty=penalty)
            
    def add_constraint(self, constraint):
        """
        Add the given constraint to the problem this Model is attached to.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        
        self.problem.add_constraint(constraint)
    
    def add_penalty(self, penalty): 
        """
        Add the given penalty term to the problem this Model is attached to.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        
        self.problem.add_penalty(penalty)
        
    def add_allele_group(self, allele_group):
        """
        Given an AlleleGroup with its reference position defined, put it in the
        Model.
        
        Sticks it at the end of the appropriate allele group list, so the list
        may no longer be sorted.
        
        """
        
        # Put the allele group in the right list
        self.allele_groups[allele_group.contig].append(allele_group)
        
    def get_allele_groups(self):
        """
        Iterate through all AlleleGroups in the model.
        
        """
        
        for allele_group_list in self.allele_groups.itervalues():
            # For each list of AlleleGroups on a contig
            for allele_group in allele_group_list:
                # For each AlleleGroup on that contig, yield it.
                yield allele_group
        
        
        
    def sort_allele_groups(self):
        """
        Sort all the Model's allele group lists into genome coordinate order.
        
        """
        
        for allele_group_list in self.allele_groups.itervalues():
            # Sort each list in place, by start coordinate and then end
            # coordinate
            allele_group_list.sort(key=lambda group: (group.start, group.end))
            
    def build_aliasing_model(self):
        """
        Takes a Model and adds measurement aliasing LP variables to it.
        
        Returns a dict of LP variables representing total hybridization (where
        1.0 = 1 perfect match copy) by probe name.
        
        Assumes we're using microarray data, where MeasurementLocations are
        probe mappings.
        
        """
        
        # This dict holds probe hybridization expressions by probe name
        hybridization_expressions = collections.defaultdict(list) 
        
        logging.info("Scanning for mapping locations")
        
        for allele_group in self.get_allele_groups():
            for measurement in allele_group.measurements:
                # Work out this measurement location's hybridization
                # contribution expression
                contribution = (measurement.get_affinity() *
                    allele_group.ploidy)
                
                # Add it in to the appropriate probe hybridization expression
                hybridization_expressions[measurement.measurement] = \
                    hybridization_expressions[measurement.measurement] + \
                    contribution
            
        # This holds hybridization LP variables by probe name
        hybridizations = {}
        
        logging.info("Creating hybridization variables")
                
        for probe, hybridization_expression in \
            hybridization_expressions.iteritems():
            
            # Make an LP variable for each of these hybridization expressions
            probe_hybridization = pulp.LpVariable(
                "{}_hybridization_{}".format(probe, get_id()), 0)
                
            # Set it equal to its expression
            self.add_constraint(hybridization_expression == 
                probe_hybridization)
            
            # Add the LP variable to the dict we're building
            hybridizations[probe] = probe_hybridization
        
        # Give back the completed dict of probe hybridization LP variables
        return hybridizations
            
    def add_constraints(self, genome_length, dna_penalty=0, end_penalty=1E9,
        default_ploidy=2, breakpoint_weight=1):
        """
        Must be called after all add_allele_group() calls, and after the allele
        group lists have been sorted with sort_allele_groups(). Constrains
        adjacent allele groups together, and puts those constraints in the
        attached SequenceGraphLpProblem.
        
        This function may only be called once.
        
        genome_length is the number of potential breakpoint positions (i.e.
        bases) in the genome.
        
        dna_penalty is how much to charge per base of gained/lost DNA.
        
        end_penalty is the penalty to charge per duplicated/deleted AlleleGroup
        at the end of a contig.
        
        default_ploidy specifies what copy number to charge for variations
        from.
        
        breakpoint_weight is a factor by which to multiply the distance-
        determined breakpoint penalty.
        
        """
        
        for allele_group_list in self.allele_groups.itervalues():
            # This holds the previous allele group we visited
            last_group = None
            
            for group in allele_group_list:
                if last_group is None:
                    # This is the first group. Constrain it to the default copy
                    # number (i.e. the assumed number of telomeres).
                    
                    self.constrain_approximately_equal(group.ploidy, 
                        default_ploidy, penalty=end_penalty)
                    
                else:
                
                    # What's the end-end distance from the last thing? Since the
                    # end location is really 1-past-the-end, we need to subtract
                    # 1 from it to get the last base actually included in the
                    # previous region.
                    distance = group.start - (last_group.end - 1)
                
                    if distance <= 0:
                        # There is at least 1 base of overlap. They must be
                        # equal.
                        self.add_constraint(group.ploidy == last_group.ploidy)
                    else:
                        # There is no overlap. Allow them to vary.
                        
                        # We want to charge for the gained/lost DNA.
                        
                        # How much do we need to charge for a copy number change
                        # on either side of this stretch of DNA? Assume we have
                        # to charge for half the bases between the centers.
                        penalty = float(group.start + group.end - 
                            last_group.start - last_group.end) / 2 * dna_penalty
                        
                        # If the last probe changes copy number from 2, charge
                        # it for half the DNA between it and here.
                        self.constrain_approximately_equal(last_group.ploidy, 
                            default_ploidy, penalty=penalty)
                        # Charge us for the other half
                        self.constrain_approximately_equal(group.ploidy, 
                            default_ploidy, penalty=penalty)
                        
                        # The copy number change penalty is the negative log
                        # likelihood of a break between them, which is the
                        # negative log of the probability that a breakpoint
                        # would hit here rather than anywhere else in the
                        # genome. We're really just adding/subtracting a
                        # constant to the breakpoint_penalty, so we take it from
                        # the user so that our breakpoint penalty is consistent
                        # for different probe sets. TODO: We don't properly
                        # account for longer genomes picking up more
                        # breakpoints.
                        breakpoint_penalty = -math.log(float(distance) / 
                            genome_length) * breakpoint_weight
                        
                        # Add the breakpoint constriant: charge the breakpoint
                        # penalty for any copy number difference between the
                        # previous allele group and this one.
                        self.constrain_approximately_equal(last_group.ploidy, 
                            group.ploidy, penalty=breakpoint_penalty)
                
                # Now we're done with this group, so it's the last
                last_group = group
                
            if len(allele_group_list) > 1:
                # Constrain the last allele group in a contig (if it isn't also
                # the first) to the default ploidy, so contigs have a cost to
                # duplicate or delete in their entirety.
                self.constrain_approximately_equal(last_group.ploidy, 
                    default_ploidy, penalty=end_penalty)
                    
    def save_map(self, stream, measured_ratios=None, 
        reference_hybridizations=None):
        """
        Save this Model to the given stream as a human-readable copy number map.
        
        The SequenceGraphLpProblem that the Model is attached to must be solved.
        
        measured_ratios is an optional dict of the sample/normal ratios for each
        MeasurementLocation measurement name (i.e. probe name).
        
        reference_hybridizations is an optional dict of the total effective
        perfect match copy number of the reference genome for each probe.
        
        """
        
        for contig, allele_groups in self.allele_groups.iteritems():
            # Print the normal copy number for each allele group
            stream.write("{}:\n".format(contig))
            for allele_group in allele_groups:
                
                # Print the copy number for the region. These may disagree where
                # probes overlap.
                stream.write("\t{} ({}:{}-{})".format(allele_group.get_ploidy(),
                    allele_group.contig, allele_group.start, allele_group.end))
                
                for mapping in allele_group.measurements:
                    # Put the probe that maps here    
                    stream.write(" ")
                    stream.write(mapping.measurement)
                    
                    
                    if measured_ratios is not None: 
                        #  Put the ratio
                        stream.write(" R: {}".format(measured_ratios[
                            mapping.measurement]))
                    
                    if reference_hybridizations is not None:
                        # Put the hybridization here vs. the total
                        stream.write(" A:{}/{}".format(mapping.get_affinity(),
                        reference_hybridizations[mapping.measurement]))
                    
                stream.write("\n")
                
    def save_bedgraph(self, stream, track_name="copyNumber",
        track_description=None):
        """
        Save the Model's copy numbers as a BEDgraph-format file (BED file where
        the feature names are really values to plot, with an appropriate browser
        track header) written to the given stream.
        
        The SequenceGraphLpProblem that the Model is attached to must be solved.
        
        track_name is an optional name for the track in the browser.
        
        track_description is an optional description for the track, also
        displayed in the browser.
        
        Automatically quotes name and description.
        
        """
        
        # If there's no description, use the track's name
        if track_description is None:
            track_description = track_name
        
        # Find the maximum copy number used, or use 2 if it's smaller than that.
        max_copy_number = 2
        
        for allele_groups in self.allele_groups.itervalues():
            for allele_group in allele_groups:
                # Scan through all the allele groups looking for the highest
                # ploidy we called anywhere.
                if allele_group.get_ploidy() > max_copy_number:
                    # This is the new highest copy number observed.
                    max_copy_number = allele_group.get_ploidy()
        
        # Save the copy number map as a BEDgraph.
        # Start with the appropriate header.
        stream.write("track type=bedGraph name=\"{}\" description=\"{}\" "
            "viewLimits=0:{} autoScale=off yLineMark=2 yLineOnOff=on\n".format(
            track_name, track_description, max_copy_number))
        
        for contig, allele_groups in self.allele_groups.iteritems():
            for group in allele_groups:
                
                # Say that the copy number under the group is what we've
                # guessed.
                stream.write("{}\t{}\t{}\t{}\n".format(contig, 
                    group.start, group.end, group.get_ploidy()))
        
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
