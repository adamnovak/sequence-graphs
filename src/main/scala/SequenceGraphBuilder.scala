package edu.ucsc.genome
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

// We want to write files
import java.io._

// We need to work with Avro things
import org.apache.avro.generic.IndexedRecord

// import spark
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

// import parquet
import parquet.hadoop.{ParquetOutputFormat, ParquetInputFormat}
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, 
                     AvroReadSupport, AvroParquetWriter}

// We need to make Paths for Parquet output.
import org.apache.hadoop.fs.Path

// And we use a hack to get at the (static) schemas of generic things
import org.apache.avro.Schema

// import hadoop job
import org.apache.hadoop.mapreduce.Job

// And Spark stuff
import org.apache.spark.rdd.RDD

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
 * SequenceGraphBuilder: an interface that allows the streaming production of
 * sequence graphs. Keeps track of the current end of every chromosome in any
 * number of "phases", and allows new Alleles and Anchors to be appended to
 * chromosomes and phases. Also allows the most recently added things to be
 * peeked at, and handles the creation of telomere Sides.
 */
trait SequenceGraphBuilder {

    /**
     * Append a new AlleleGroup holding the given allele to all of the given
     * phases of the given contig. The total ploidy is exactly equal to the
     * number of phases to which the allele is being appended. The AlleleGroup
     * starts at the given reference start position, but the caller is
     * responsible for adding any Anchors that are necessary beforehand. The
     * referenceLength parameter specifies the length of the reference Site
     * which the AlleleGroup is to occupy (i.e. how far its ending Side should
     * be from its starting Side); by default this is the number of bases in the
     * given Allele.
     * 
     * Phases must not be empty. The first phase specified determines the
     * starting position of the Site this AlleleGroup occupies. All phases
     * specified must currently end with `Face.RIGHT` Sides.
     */
    def addAllele(contig: String, phases: Seq[Int], start: Long, allele: Allele, 
        referenceLength: Long = -1)
        
    /**
     * Add an Anchor, with ploidy equal to the number of phases specified here,
     * to the given phases of the given contigs. The anchor will start at the
     * given start position, and run for the specified number of bases.
     *
     * If referenceLength is 0, does nothing.
     *
     * Phases must not be empty. The first phase specified determines the
     * starting position of the Site this AlleleGroup occupies. All phases
     * specified must currently end with `Face.RIGHT` Sides.
     *
     * TODO: Unify somewhow with addAllele.
     */
    def addAnchor(contig: String, phases: Seq[Int], start: Long,
        referenceLength: Long)
        
    /**
     * Attach Anchors to the ends of some of the specified phases until all the
     * phases end at the same reference position. Returns that reference
     * position.
     */
    def squareOff(contig: String, phases: Seq[Int]) : Long
        
    /**
     * Add a trailing telomere to the given phase of the given contig, closing
     * it off to the given length. Attach to the last thing we have there
     * with an Adjacency. Will create a new leading telomere if nothing is in
     * that phase of that contig already.
     */
    def close(contig: String, phase: Int, length: Long)
        
    /**
     * Get the last Side on the given contig in the given phase, or add and
     * remember a new leading telomere if none is found.
     */
    def getLastSide(contig: String, phase: Int) : Side
    
    /**
     * Called when the sequence graph is completed, and should be made ready for
     * reading from.
     */
    def finish()
}

/**
 * Represents a SequenceGraphBuilder that, in addition to taking Alleles and
 * Anchors on phases, can support wiring up nonreference adjacencies.
 */
trait StructuralSequenceGraphBuilder extends SequenceGraphBuilder {

    /**
     * Add a deletion on the specified phases of the specified contig, running
     * from the given start position and consuming the specified number of
     * bases.
     *
     * The next thing added to any of the given phases will be placed after the
     * end of the deletion. It is the caller's responsibility not to add the
     * things that should have been deleted, and deal with e.g. deletions which
     * end in the middle of where the caller might want to put, say, an unphased
     * AlleleGroup.
     */
    def addDeletion(contig: String, phases: Seq[Int], start: Long, length: Long)
    
    // Insertions are handled through addAllele, using the referenceLength
    // parameter. If the sequence is unknown, the caller has to make an ALlele
    // with a run of Ns.
    
    /**
     * Add an inversion of the specified length to the specified phases of the
     * specified contig. Should be called when the current end of the contig is
     * at the start of the inversion. When enough AlleleGroups and Anchors have
     * been added to reach the end of the inversion, the Adjacencies defining
     * the inversion will be added.
     *
     * Calls to addAnchor and addAllele will be intercepted and modified in
     * order to break AlleleGroups or Anchors that would cover the endpoint of
     * the inversion into two pieces, to provide an attachment point for the new
     * adjacencies.
     */
    def addInversion(contig: String, phases: Seq[Int], start: Long,
        length: Long)
    
    /**
     * Add a tandem duplication of the specified length to the specified phases
     * of the specified contig. Should be called when the cirrent end of the
     * contig is at the start of the region to be duplicated. When enough
     * AlleleGroups and Anchors have been added to reach the end of the
     * duplication, the Adjacencies defining the duplication will be added.
     *
     * Note that a tandem duplication constructed in this way cannot be
     * distinguished from an additional circular chromosome.
     *
     * Calls to addAnchor and addAllele will be intercepted and modified in
     * order to break AlleleGroups or Anchors that would cover the endpoint of
     * the duplication into two pieces, to provide an attachment point for the
     * new adjacencies.
     */
    def addDuplication(contig: String, phases: Seq[Int], start: Long,
        length: Long)
    
}

/**
 * EasySequenceGraphBuilder: Implement a StructuralSequenceGraphBuilder, writing
 * elements to the given SequenceGraphWriter.
 */
class EasySequenceGraphBuilder(sample: String, reference: String, 
    writer: SequenceGraphWriter) extends StructuralSequenceGraphBuilder {
    
    // This mutable HashMap holds all the Sides at the ends of chromosomes, by
    // contig name and phase number (usually 0 or 1).
    private val ends = HashMap.empty[(String, Int), Side]
    
    // These are the useful methods involved in implementing the
    // SequenceGraphBuilder interface.
    
    /**
     * Set the last Side of the given phase of the given contig to the given
     * Side. If Side is null, means a trailing telomere has been added and
     * nothing more can be added to that phase.
     */
    protected def setLastSide(contig: String, phase: Int, side: Side) {
        ends((contig, phase)) = side
    }
    
    /**
    * Make a Side that could be the next Side on the given contig, advancing
    * ahead to the given reference position.
    */
    def getNextSide(contig: String, phase: Int, position: Long) : Side = {
        // Make and return a new Side that comes at the given position.
        new Side(IDMaker.get(), new Position(contig, position, Face.LEFT),
            false)
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
        
        if(end.position.face != Face.RIGHT) {
            // Sanity-check the trailing base, so we don't accidentally go
            // backwards in the reference.
            throw new Exception("A 5' Face is trailing contig %s phase %d"
                .format(contig, phase))
        }
        
        // Advance the base
        val newBase = end.position.base + 1
         
        getNextSide(contig, phase, newBase)
    }
    
    /**
     * Attach the given AlleleGroup to the end of the given contig's given
     * phase. The caller is responsible for also setting the last Side of that
     * contig and phase with `setLastSide`, and writing the AlleleGroup to the
     * SequenceGraphWriter.
     */
    protected def connectAlleleGroup(contig: String, phase: Int, 
        alleleGroup: AlleleGroup) : Unit = {
        
        ends.get((contig, phase)).foreach { (end) => 
            // If we do have something at the end of this phase of this contig
            // already, make an Adjacency to this AlleleGroup's first Side.
            val newAdjacency = Adjacency.newBuilder()
                // Attach the Edge
                .setEdge(new Edge(IDMaker.get(), end.id, alleleGroup.edge.left))
                // Set ploidy to exactly 1
                .setPloidy(new PloidyBounds(1, null, null))
                // Attach to our genome
                .setGenome(sample)
                .build()
                
            // Add the Adjacency to our collection
            writer.writeAdjacency(newAdjacency)
            
        }
        
    }
    
    /**
     * Attach the given Anchor to the end of the given contig's given phase. The
     * caller is responsible for also setting the last Side of that contig and
     * phase with `setLastSide`, and writing the Anchor to the
     * SequenceGraphWriter.
     * 
     * TODO: Unify somehow with connectAlleleGroup.
     */
    protected def connectAnchor(contig: String, phase: Int, 
        anchor: Anchor) : Unit = {
        
        ends.get((contig, phase)).foreach { (end) => 
            // If we do have something at the end of this phase of this contig
            // already, make an Adjacency to this AlleleGroup's first Side.
            val newAdjacency = Adjacency.newBuilder()
                // Attach the Edge
                .setEdge(new Edge(IDMaker.get(), end.id, anchor.edge.left))
                // Set ploidy to exactly 1
                .setPloidy(new PloidyBounds(1, null, null))
                // Attach to our genome
                .setGenome(sample)
                .build()
                
            // Add the Adjacency to our collection
            writer.writeAdjacency(newAdjacency)
        }
    }
    
    // These are the SequenceGraphBuilder method implementations
    
    def squareOff(contig: String, phases: Seq[Int]) : Long = {
        // Where do all the phases end?
        val ends = phases.map(getLastSide(contig, _).position.base)
        
        // Work out what position we need to square off up to.
        val targetPosition = ends.max
            
        // Put phased Anchors up to there.
        phases.zip(ends).map { (pair) =>
            // Unpack the pair
            val (phase, end) = pair
            
            if(end < targetPosition) {
                // Pad up to TargetPosition with an anchor
                addAnchor(contig, List(phase), end + 1, targetPosition - end)
            }
        }
        
        // Return the position where all those phases now end.
        targetPosition
    }
    
    def getLastSide(contig: String, phase: Int) : Side = {
        ends.get((contig, phase)) getOrElse {
            // We couldn't find anything there already. Make a new (vacuously
            // phased) telomer Side.
            val telomere = new Side(IDMaker.get(), new Position(contig, 0, 
                Face.RIGHT), false)
                
            // Remember the telomere
            writer.writeSide(telomere)
            
            // Stick it at the end where it goes
            setLastSide(contig, phase, telomere)
            
            // Return it
            telomere
        }
    }
    
    def addAllele(contig: String, phases: Seq[Int], start: Long, allele: Allele, 
        referenceLength: Long = -1) = {
        
        if(phases.size > 1) {
            // Make sure we have an even starting point to add on to (so we
            // don't create deletions that aren't real).
            squareOff(contig, phases)
        }
        
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
        val leadingSide = getNextSide(contig, phases.head, start)
        
        // Make a right Side for the AlleleGroup (non-reference). We subtract 1
        // from the position because a 1-base AlleleGroup should have the same
        // numbers on its two Sides (one is left, and the other is right)
        val trailingSide = new Side(IDMaker.get(), new Position(contig, 
            leadingSide.position.base + actualReferenceLength - 1, Face.RIGHT),
            false)
            
            
        // Make an AlleleGroup with the correct ploidy to be added to that many
        // phases.
        val alleleGroup = new AlleleGroup(new Edge(IDMaker.get(), 
            leadingSide.id, trailingSide.id), allele, ploidy, sample)
        
        // Add the new Sides
        writer.writeSide(leadingSide)
        writer.writeSide(trailingSide)
        
        // Add the AlleleGroup itself
        writer.writeAlleleGroup(alleleGroup)
        
        phases map { (phase) =>
            // Connect the AlleleGroup into each Phase with a ploidy-1 Adjacency.
            connectAlleleGroup(contig, phase, alleleGroup)
            // Set the last Side for that phase
            setLastSide(contig, phase, trailingSide)
        }
        
        
    }
    
    def addAnchor(contig: String, phases: Seq[Int], start: Long,
        referenceLength: Long) : Unit = {
        
        if(phases.size > 1) {
            // Make sure we have an even starting point to add on to (so we
            // don't create deletions that aren't real).
            squareOff(contig, phases)
        }
        
        if(referenceLength >= 0) {
            // 0-length Anchors are allowed, but negative-length ones are not.
        
            // Ploidy is number of phases to append to.
            val ploidy = new PloidyBounds(phases.size, null, null)
            
            // Ensure a telomere exists for all phases
            phases.map(getLastSide(contig, _))
            
            // Get the left Side for the new Anchor. It can be generated from
            // any phase. We know it will be a Face.LEFT Side because of the
            // preconditions on this method.
            val leadingSide = getNextSide(contig, phases.head, start)
            
            // Make a right Side for the Anchor. We advance by 1 less than the
            // specified reference length because of the left/right distinction
            // for Sides: a 1-base anchor has both its Sides at the same
            // numerical position, on different faces.
            val trailingSide = new Side(IDMaker.get(), new Position(contig, 
                leadingSide.position.base + referenceLength - 1, Face.RIGHT),
                false)
                
                
            // Make an Anchor with the correct ploidy to be added to that many
            // phases.
            val anchor = new Anchor(new Edge(IDMaker.get(), 
                leadingSide.id, trailingSide.id), ploidy, sample)
            
            // Add the new Sides
            writer.writeSide(leadingSide)
            writer.writeSide(trailingSide)
            
            // Add the Anchor itself
            writer.writeAnchor(anchor)
            
            phases map { (phase) =>
                // Connect the Anchor into each Phase with a ploidy-1 Adjacency.
                connectAnchor(contig, phase, anchor)
                // Set the last Side for that phase
                setLastSide(contig, phase, trailingSide)
            }
        }
    }
    
    def close(contig: String, phase: Int, length: Long) {
        // Get the side we have to come after, or a new leading telomere
        val end = getLastSide(contig, phase)
        
        // Create a new side that ought to be at the end of the given contig.
        val telomere = getNextSide(contig, phase, length)

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
        writer.writeSide(telomere)
        writer.writeAdjacency(adjacency)
        setLastSide(contig, phase, null)
    }
    
    // When we finish, close the writer.
    def finish() = writer.close
    
    // Structural variants (StructuralSequenceGraphBuilder methods)
    
    def addDeletion(contig: String, phases: Seq[Int], start: Long, 
        length: Long) = {
        // The next thing we add to any of these Phases needs to be at least
        // that far away. We accomplish this by putting a dummy goes-backwards
        // empty Anchor at the end of the deletion.
        
        // TODO: Do it without the anchor.
        
        phases.map { (phase) =>
        
            // Make sure we have a leading telomere.
            getLastSide(contig, phase)
            
            // Get the left Side for the new Anchor, skipping the length of the
            // deletion.
            println("Skipping %d to %d".format(start, start + length))
            val leadingSide = getNextSide(contig, phase, start + length)
            
            // Make a right Side for the Anchor. It will have the numerical
            // position of 1 *before* that of the left side of the Anchor.
            val trailingSide = new Side(IDMaker.get(), new Position(contig, 
                leadingSide.position.base - 1, Face.RIGHT), false)
                
            // Make an Anchor with ploidy 1
            val ploidy = new PloidyBounds(1, null, null)
            val anchor = new Anchor(new Edge(IDMaker.get(), 
                leadingSide.id, trailingSide.id), ploidy, sample)
            
            // Add the new Sides
            writer.writeSide(leadingSide)
            writer.writeSide(trailingSide)
            
            // Add the Anchor itself
            writer.writeAnchor(anchor)
            
            // Connect the Anchor into this phase with a ploidy-1 Adjacency.
            connectAnchor(contig, phase, anchor)
            // Set the last Side for this phase
            setLastSide(contig, phase, trailingSide)
        }
    }
    
    def addInversion(contig: String, phases: Seq[Int], start: Long,
        length: Long) = {}
    
    def addDuplication(contig: String, phases: Seq[Int], start: Long,
        length: Long) = {}
}








