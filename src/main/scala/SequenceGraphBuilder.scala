package edu.ucsc.genome
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

// We want to write GraphViz graphs
import org.kohsuke.graphviz._

// We want to write files
import java.io._

/**
 * This object handles providing sequential global IDs.
 */
object IDMaker {
    // This holds the next available ID
    var next : Long = 0
    
    /**
     * Get a unique ID.
     */
    def get() : Long = {
        // Make an ID to return
        val id = next
        // Advance so we generate a different ID next time.
        next += 1
        // Return the ID we generated
        id
    }
}

/**
 * Dummy class for Anchors, which are not yet in the Avro spec. Anchors really
 * represent long runs of the reference haplotype, and are most like
 * AlleleGroups, but we base them off Adjacencies here since that class has all
 * the same members.
 *
 * TODO: Add Anchor back to the Avro spec.
 */
class Anchor(a: Edge, b: PloidyBounds, c: String) extends Adjacency(a, b, c)

/**
 *
 * SequenceGraphBuilder: a class to build up a sequence graph for a diploid
 * genome.
 * 
 * This class contains collections of all the parts needed to build a proper
 * sequence graph.
 * 
 * It keeps track of the last Side for each copy of each chromosome (which may
 * be the same Side for both copies of a chromosome in an area with no phasing
 * information). New AlleleGroups can be added to the end of each copy of each
 * chromosome.
 * 
 * Operates on a given genome/sample name, and a given reference name.
 *
 */
class SequenceGraphBuilder(sample: String, reference: String) {
    // This holds all the Sides we have created, by ID
    val sides = HashMap.empty[Long, Side]
    
    // This holds all the AlleleGroups we have created.
    val alleleGroups = HashSet.empty[AlleleGroup]
    
    // This holds all the Adjacencies we have created.
    val adjacencies = HashSet.empty[Adjacency]
    
    // This holds all the Anchors we have created.
    val anchors = HashSet.empty[Anchor]
    
    // This holds all the IDs of the Sides at the ends of chromosomes, by contig
    // name and phase number (usually 0 or 1).
    val ends = HashMap.empty[(String, Int), Long]
    
    /**
     * Attach the given AlleleGroup to the end of the given contig's given
     * phase.
     */
    def addAlleleGroup(contig: String, phase: Int, alleleGroup: AlleleGroup) = {
        
        ends.get((contig, phase)).foreach { (end) => 
            // If we do have something at the end of this phase of this contig
            // already, make an Adjacency to this AlleleGroup's first Side.
            val newAdjacency = Adjacency.newBuilder()
                // Attach the Edge
                .setEdge(new Edge(IDMaker.get(), end, alleleGroup.edge.left))
                // Set ploidy to exactly 1
                .setPloidy(new PloidyBounds(1, null, null))
                // Attach to our genome
                .setGenome(sample)
                .build()
                
            // Add the Adjacency to our collection
            adjacencies += newAdjacency
            
        }
        
        // Put this AlleleGroup's second Side as the new trailing end of
        // the chromosome.
        ends((contig, phase)) = alleleGroup.edge.right
        
        // Remember the AlleleGroup in our set of AlleleGroups
        alleleGroups += alleleGroup
        
        println("Added %s to %s:%d".format(alleleGroup.allele, contig, phase))
    }
    
    /**
     * Attach the given Anchor to the end of the given contig's given
     * phase.
     * 
     * TODO: Unify somehow with addAlleleGroup.
     */
    def addAnchor(contig: String, phase: Int, 
        anchor: Anchor) : Unit = {
        
        ends.get((contig, phase)).foreach { (end) => 
            // If we do have something at the end of this phase of this contig
            // already, make an Adjacency to this AlleleGroup's first Side.
            val newAdjacency = Adjacency.newBuilder()
                // Attach the Edge
                .setEdge(new Edge(IDMaker.get(), end, anchor.edge.left))
                // Set ploidy to exactly 1
                .setPloidy(new PloidyBounds(1, null, null))
                // Attach to our genome
                .setGenome(sample)
                .build()
                
            // Add the Adjacency to our collection
            adjacencies += newAdjacency
            
        }
        
        // Put this Anchor's second Side as the new trailing end of
        // the chromosome.
        ends((contig, phase)) = anchor.edge.right
        
        // Remember the Anchor in our set of Anchors
        anchors += anchor
        
        println("Added anchor to %s:%d".format(contig, phase))
    }
    
    /**
     * Append a new AlleleGroup holding the given allele to all of the given
     * phases of the given contig. The total ploidy is exactly equal to the
     * number of phases to which the allele is being appended. The
     * referenceLength parameter specifies the length of the reference Site
     * which the AlleleGroup is to occupy (i.e. how far its ending Side should
     * be from its starting Side); by default this is the number of bases in the
     * given Allele.
     * 
     * Phases must not be empty. The first phase specified determines the
     * starting position of the Site this AlleleGroup occupies. All phases
     * specified must currently end with `Face.RIGHT` Sides.
     */
    def addAllele(contig: String, phases: List[Int], allele: Allele, 
        referenceLength: Int = -1) = {
        
        // Ploidy is number of phases to append to.
        val ploidy = new PloidyBounds(phases.size, null, null)
        // Reference length defaults to number of bases in allele.
        val actualReferenceLength  = if(referenceLength == -1) {
            // If it's -1 (the default), just use the length of the allele.
            allele.bases.size
        } else {
            // Otherwise use the specified value.
            referenceLength
        }
        
        // Get the left Side for the new AlleleGroup. It can be generated from
        // any phase. We know it will be a Face.LEFT Side because of the
        // preconditions on this method.
        val leadingSide = getNextSide(contig, phases.head)
        
        // Make a right Side for the AlleleGroup (non-reference)
        val trailingSide = new Side(IDMaker.get(), new Position(contig, 
            leadingSide.position.base + actualReferenceLength, Face.RIGHT),
            false)
            
            
        // Make an AlleleGroup with the correct ploidy to be added to that many
        // phases.
        val alleleGroup = new AlleleGroup(new Edge(IDMaker.get(), 
            leadingSide.id, trailingSide.id), allele, ploidy, sample)
        
        // Add the new Sides
        addSide(leadingSide)
        addSide(trailingSide)
        
        phases map { (phase) =>
            // Add the AlleleGroup into each Phase with a ploidy-1 Adjacency.
            addAlleleGroup(contig, phase, alleleGroup)
        }
    }
    
    /**
     * Add an Anchor, with ploidy equal to the number of phases specified here,
     * to the given phases of the given contigs. The anchor will run for the
     * specified number of bases.
     *
     * If referenceLength is 0, does nothing.
     *
     * Phases must not be empty. The first phase specified determines the
     * starting position of the Site this AlleleGroup occupies. All phases
     * specified must currently end with `Face.RIGHT` Sides.
     *
     * TODO: Unify somewhow with addAllele.
     */
    def addAnchor(contig: String, phases: List[Int], 
        referenceLength: Int) : Unit = {
        
        if(referenceLength > 0) {
        
            // Ploidy is number of phases to append to.
            val ploidy = new PloidyBounds(phases.size, null, null)
            
            // Get the left Side for the new Anchor. It can be generated from
            // any phase. We know it will be a Face.LEFT Side because of the
            // preconditions on this method.
            val leadingSide = getNextSide(contig, phases.head)
            
            // Make a right Side for the Anchor
            val trailingSide = new Side(IDMaker.get(), new Position(contig, 
                leadingSide.position.base + referenceLength, Face.RIGHT), false)
                
                
            // Make an Anchor with the correct ploidy to be added to that many
            // phases.
            val anchor = new Anchor(new Edge(IDMaker.get(), 
                leadingSide.id, trailingSide.id), ploidy, sample)
            
            // Add the new Sides
            addSide(leadingSide)
            addSide(trailingSide)
            
            phases map { (phase) =>
                // Add the Anchor into each Phase with a ploidy-1 Adjacency.
                addAnchor(contig, phase, anchor)
            }
        }
    }
    
    /**
     * Remember the given Side by its ID. Can be called repeatedly on the same
     * Side.
     */
    def addSide(side: Side) = {
        // Put the Side in the HashMap
        sides(side.id) = side
    }
    
    /**
     * Get the last Side on the given contig in the given phase, or add and
     * remember a new leading telomere if none is found.
     */
    def getLastSide(contig: String, phase: Int) : Side = {
        ends.get((contig, phase)).flatMap(sides.get) getOrElse {
            // We couldn't find anything there already. Make a new (vacuously
            // phased) telomer Side.
            val telomere = new Side(IDMaker.get(), new Position(contig, 0, 
                Face.RIGHT), false)
                
            // Remember the telomere
            addSide(telomere)
            
            // Stick it at the end where it goes
            ends((contig, phase)) = telomere.id
            
            // Return it
            telomere
        }
    }
    
    /**
     * Get the last AlleleGroup on the given phase of the given contig, or None
     * if there is no AlleleGroup there.
     */
    def getLastAlleleGroup(contig: String, phase: Int) : Option[AlleleGroup] = {
        // Go get the last Side ID on that phase of that contig.
        ends.get((contig, phase)).flatMap({ (end) =>
            // Find an AlleleGroup that has this ID as one of its edge's ends.
            // There ought to be at most one.
            alleleGroups.find({ (alleleGroup) =>
                (alleleGroup.edge.left == end) || 
                (alleleGroup.edge.right == end)
            })
        })
    }
    
    /**
     * Create and return (but do not remember) a Side corresponding to the 5'
     * face of the next unaccounted-for base of the given phase of the given
     * contig.
     *
     * If there is nothing at the end of that phase of that contig, adds a
     * telomere first.
     */
    def getNextSide(contig: String, phase: Int) : Side = {
        // Go get the last Side, adding a telomere if necessary.
        val end = getLastSide(contig, phase)
        
        // Flip to the opposite face
        val newFace = end.position.face match {
            case Face.LEFT => Face.RIGHT
            case Face.RIGHT => Face.LEFT
        }
        
        // Advance or un-advance the base
        val newBase = end.position.face match {
            case Face.LEFT => end.position.base - 1
            case Face.RIGHT => end.position.base + 1
        }
         
        // Make and return a new Side that comes directly after the last one
        // (which may have been a leading telomere)
        new Side(IDMaker.get(), new Position(contig, newBase, newFace), false)
    }
    
    /**
     * Add a trailing telomere to the given phase of the given contig. Attach it
     * to the last thing we have there with an Adjacency. Will create a new
     * leading telomere if nothing is in that phase of that contig already.
     */
    def close(contig: String, phase: Int) {
        // Get the side we have to come after, or a new leading telomere
        val end = getLastSide(contig, phase)
        
        // Create a new side that ought to be at the end of the given contig.
        val telomere = getNextSide(contig, phase)

        // Make an Adjacency to link them up
        val adjacency = Adjacency.newBuilder()
            // Attach the Edge
            .setEdge(new Edge(IDMaker.get(), end.id, telomere.id))
            // Set ploidy to exactly 1
            .setPloidy(new PloidyBounds(1, null, null))
            // Attach to our genome
            .setGenome(sample)
            .build()
            
        // Remember everything
        addSide(telomere)
        adjacencies += adjacency
    }
    
    /**
     * Write out the graph we have built to the given GraphViz dot file.
     */
    def writeDotFile(file: String) {
        // Make a new graph
        val graph = new org.kohsuke.graphviz.Graph()
        
        // Set up edge styles
        graph.edgeWith(new Style().attr("arrowsize", "0"))
        
        // And node styles
        graph.nodeWith(new Style().attr("shape", "point").attr("label", ""))
        
        // Make a Node for every Side, organized by ID, and add to the graph
        val nodes = sides.mapValues({ (side) =>
            val node = new Node().id(side.id.toString)
            graph.node(node)
            node
        })
        
        alleleGroups map { (alleleGroup : AlleleGroup) =>
            
            // What should we label it?
            val label = alleleGroup.allele match {
                // Put bases if we have actual bases
                case allele: Allele => allele.bases
                // Put something else if we're referencing an allele by index
                case index: java.lang.Integer => "#%d".format(index)
                // It's something else
                case _ => "<unknown>"
            }
            
            // Make an edge for every AlleleGroup.
            graph.edge(nodes(alleleGroup.edge.left),
                nodes(alleleGroup.edge.right), new Style()
                    // Label it with the bases
                    .attr("label", label)
                    // Give it an arrow head
                    .attr("arrowsize", "1")
                    .attr("color", "#ff0000"))
        }
        
        adjacencies map { (adjacency) =>
            // Make an edge for every Adjacency
            graph.edge(nodes(adjacency.edge.left),
                nodes(adjacency.edge.right))
        }
        
        anchors map { (anchor) =>
            // What are the two Sides of the anchor?
            val leftPos = sides(anchor.edge.left).position
            val rightPos = sides(anchor.edge.right).position
            
            // What should we write on the anchor?
            val label = if(leftPos.contig == rightPos.contig) {
                // The anchor is properly on a single contig, so we can
                // determine its length.
                val length = rightPos.base - leftPos.base
                "%sbp anchor".format(length)
            } else {
                // We have no idea how long the anchor is
                "anchor across contigs"
            }
            
            // Make an edge for every Anchor
            graph.edge(nodes(anchor.edge.left),
                nodes(anchor.edge.right), new Style()
                    .attr("label", label)
                    .attr("color", "#0000ff"))
        }
        
        // Write the graph to the file
        graph.writeTo(new FileOutputStream(file))
        
        println("Wrote %s".format(file))
        
    }
    
    
}








