package edu.ucsc.genome
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

// We want to write GraphViz graphs
import org.kohsuke.graphviz._

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
                     AvroReadSupport}

// And we use a hack to get at the (static) schemas of generic things
import org.apache.avro.Schema

// import hadoop job
import org.apache.hadoop.mapreduce.Job

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
    def addAllele(contig: String, phases: Seq[Int], allele: Allele, 
        referenceLength: Int = -1)
        
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
    def addAnchor(contig: String, phases: Seq[Int], 
        referenceLength: Int)
        
    /**
     * Add a trailing telomere to the given phase of the given contig. Attach it
     * to the last thing we have there with an Adjacency. Will create a new
     * leading telomere if nothing is in that phase of that contig already.
     */
    def close(contig: String, phase: Int)
        
    /**
     * Get the last Side on the given contig in the given phase, or add and
     * remember a new leading telomere if none is found.
     */
    def getLastSide(contig: String, phase: Int) : Side
    
    /**
     * Called when the sequence graph is completed, and should be made ready for
     * reading from.
     */
    def finish() = {}
}

/**
 * EasySequenceGraphBuilder: Implement most of the functionality of a
 * SequenceGraphBuilder, leaving just the part about actually providing a way to
 * write graph elements.
 */
abstract class EasySequenceGraphBuilder(sample: String, reference: String)
    extends SequenceGraphBuilder {
    
    // This mutable HashMap holds all the Sides at the ends of chromosomes, by
    // contig name and phase number (usually 0 or 1).
    private val ends = HashMap.empty[(String, Int), Side]
    
    // These are the methods implementations need to define.
    
    /**
     * Add the given AlleleGroup to the graph. Called exactly once per
     * AlleleGroup.
     */
    protected def addAlleleGroup(alleleGroup: AlleleGroup) : Unit
    
    /**
     * Add the given Adjacency to the graph. Called exactly once per Adjacency.
     */
    protected def addAdjacency(adjacency: Adjacency) : Unit
    
    /**
     * Add the given Anchor to the graph. Called exactly once per Anchor.
     */
    protected def addAnchor(anchor: Anchor) : Unit
    
    /**
     * Add the given Side to the graph. May be called repeatedly on the same
     * Side.
     */
    protected def addSide(side: Side) : Unit
    
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
     * Attach the given AlleleGroup to the end of the given contig's given
     * phase. The caller is responsible for also setting the last Side of that
     * contig and phase with `setLastSide`, and adding the AlleleGroup to the
     * graph with `addAlleleGroup(AlleleGroup)`.
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
            addAdjacency(newAdjacency)
            
        }
        
    }
    
    /**
     * Attach the given Anchor to the end of the given contig's given phase. The
     * caller is responsible for also setting the last Side of that contig and
     * phase with `setLastSide`, and adding the Anchor to the graph with
     * `addAnchor(Anchor)`.
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
            addAdjacency(newAdjacency)
        }
    }
    
    // These are the SequenceGraphBuilder method implementations
    
    def getLastSide(contig: String, phase: Int) : Side = {
        ends.get((contig, phase)) getOrElse {
            // We couldn't find anything there already. Make a new (vacuously
            // phased) telomer Side.
            val telomere = new Side(IDMaker.get(), new Position(contig, 0, 
                Face.RIGHT), false)
                
            // Remember the telomere
            addSide(telomere)
            
            // Stick it at the end where it goes
            setLastSide(contig, phase, telomere)
            
            // Return it
            telomere
        }
    }
    
    def addAllele(contig: String, phases: Seq[Int], allele: Allele, 
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
        
        // Add the AlleleGroup itself
        addAlleleGroup(alleleGroup)
        
        phases map { (phase) =>
            // Connect the AlleleGroup into each Phase with a ploidy-1 Adjacency.
            connectAlleleGroup(contig, phase, alleleGroup)
            // Set the last Side for that phase
            setLastSide(contig, phase, trailingSide)
        }
        
        
    }
    
    def addAnchor(contig: String, phases: Seq[Int], 
        referenceLength: Int) : Unit = {
        
        if(referenceLength > 0) {
        
            // Ploidy is number of phases to append to.
            val ploidy = new PloidyBounds(phases.size, null, null)
            
            // Ensure a telomere exists for all phases
            phases.map(getLastSide(contig, _))
            
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
            
            // Add the Anchor itself
            addAnchor(anchor)
            
            phases map { (phase) =>
                // Connect the Anchor into each Phase with a ploidy-1 Adjacency.
                connectAnchor(contig, phase, anchor)
                // Set the last Side for that phase
                setLastSide(contig, phase, trailingSide)
            }
        }
    }
    
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
        addAdjacency(adjacency)
        setLastSide(contig, phase, null)
    }
}

/**
 * GraphvizSequenceGraphBuilder: build a sequence graph on disk, writing to a
 * GraphViz dot file in a streaming fashion.
 *
 * TODO: Make streaming. Actually write to disk instead of the GarphViz library.
 *
 * Build a graph for the given genome and sample names, and write it to the
 * given file.
 *
 */
class GraphvizSequenceGraphBuilder(sample: String, reference: String, 
    file: String) extends EasySequenceGraphBuilder(sample, reference) {
    
    // We need to write escaped strings. See
    // <http://stackoverflow.com/a/9914380/402891> and
    // <http://commons.apache.org/proper/commons-
    // lang/apidocs/org/apache/commons/lang3/StringEscapeUtils.html>
    import org.apache.commons.lang.StringEscapeUtils.escapeJava
    
    // Open the file for writing. See
    // <http://www.tutorialspoint.com/scala/scala_file_io.htm>
    val graphWriter = new java.io.PrintWriter(new File(file))
    
    // Write a header
    graphWriter.write("digraph %s {\n".format(sample))
    
    // Set up edge styles
    graphWriter.write("edge [arrowsize=\"0\"];\n")
    
    // And node styles
    graphWriter.write("node [shape=\"point\", label=\"\"];\n")
    
    // This holds all the Nodes we have created, by ID
    private val nodes = HashMap.empty[Long, Node]
    
    // Implementations of what we need to be an EasySequenceGraphBuilder
    
    protected def addAlleleGroup(alleleGroup: AlleleGroup) : Unit = {
        // What should we label it?
        val label = (alleleGroup.allele match {
            // Put bases if we have actual bases
            case allele: Allele => allele.bases
            // Put something else if we're referencing an allele by index
            case index: java.lang.Integer => "#%d".format(index)
            // It's something else
            case _ => "<unknown>"
        }) + "\nx%d".format(alleleGroup.ploidy.lower)
        
        // Make an edge for every AlleleGroup.
        // TODO: escape labels
        graphWriter.write(
            "%d -> %d [label=\"%s\",arrowsize=\"1\",color=\"#ff0000\"];\n"
            .format(alleleGroup.edge.left, alleleGroup.edge.right, 
            escapeJava(label)))
    }
    
    protected def addAdjacency(adjacency: Adjacency) : Unit = {
        // Make an edge for every Adjacency
        graphWriter.write("%d -> %d;\n".format(
            adjacency.edge.left, adjacency.edge.right))
    }
    
    protected def addAnchor(anchor: Anchor) : Unit ={
        // What should we write on the anchor?
        val label = "Anchor x%d".format(anchor.ploidy.lower)
        
        // Make an edge for every Anchor
        graphWriter.write(
            "%d -> %d [label=\"%s\",color=\"#0000ff\"];\n"
            .format(anchor.edge.left, anchor.edge.right, escapeJava(label)))
    }
    
    protected def addSide(side: Side) : Unit = {
        // We don't really need to define sides, but this should let any sides
        // we somehow haven't connected show up.
        graphWriter.write(
            "%d;\n"
            .format(side.id))
    }
    
    
    // SequenceGraphBuilder hooks we use
    
    /**
     * Close the GraphViz graph we have been writing.
     */
    override def finish() {
        // Close the graph block
        graphWriter.write("}\n")
        
        // Close the file
        graphWriter.close()
        
        println("Wrote %s".format(file))
    }
}


/**
 *
 * InMemorySequenceGraphBuilder: a class to build up a sequence graph for a
 * genome in memory.
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
class InMemorySequenceGraphBuilder(sample: String, reference: String) 
    extends EasySequenceGraphBuilder(sample, reference) {

    // This holds all the Sides we have created, by ID
    private val sides = HashMap.empty[Long, Side]
    
    // This holds all the AlleleGroups we have created.
    private val alleleGroups = HashSet.empty[AlleleGroup]
    
    // This holds all the Adjacencies we have created.
    private val adjacencies = HashSet.empty[Adjacency]
    
    // This holds all the Anchors we have created.
    private val anchors = HashSet.empty[Anchor]
    
    // Implementations of what we need to be an EasySequenceGraphBuilder
    
    protected def addAlleleGroup(alleleGroup: AlleleGroup) : Unit = 
        alleleGroups += alleleGroup
    
    protected def addAdjacency(adjacency: Adjacency) : Unit = 
        adjacencies += adjacency
    
    protected def addAnchor(anchor: Anchor) : Unit =
        anchors += anchor
    
    protected def addSide(side: Side) : Unit = {
        // Put the Side in the HashMap
        sides(side.id) = side
    }
    
    // Extra methods we support
    
    /**
     * Get the last AlleleGroup on the given phase of the given contig, or None
     * if there is no AlleleGroup there.
     */
    def getLastAlleleGroup(contig: String, phase: Int) : Option[AlleleGroup] = {
        // Go get the last Side ID on that phase of that contig.
        val end = getLastSide(contig, phase).id
        
        // Find an AlleleGroup that has this ID as one of its edge's ends.
        // There ought to be at most one.
        alleleGroups.find({ (alleleGroup) =>
            (alleleGroup.edge.left == end) || 
            (alleleGroup.edge.right == end)
        })
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
            val label = (alleleGroup.allele match {
                // Put bases if we have actual bases
                case allele: Allele => allele.bases
                // Put something else if we're referencing an allele by index
                case index: java.lang.Integer => "#%d".format(index)
                // It's something else
                case _ => "<unknown>"
            }) + "\nx%d".format(alleleGroup.ploidy.lower)
            
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
                "%sbp anchor x%d".format(length, anchor.ploidy.lower)
            } else {
                // We have no idea how long the anchor is
                "anchor across contigs x%d".format(anchor.ploidy.lower)
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
    
    /**
     * Given an iteravle of Avro records of type RecordType, writes them to a
     * Parquet file in the specified directory. Any other Parquet data in that
     * directory will be overwritten.
     */
    def writeCollectionToParquet[RecordType <: IndexedRecord](
        things: Iterable[RecordType], directory: String, sc: SparkContext)(implicit m: Manifest[RecordType]) = {
        
        // Every timne we switch Avro schemas, we need a new Job. So we create
        // our own.
        val job = new Job()
        
        try {
            // That SparkContext needs to have the slf4j StaticLoggerBinder we
            // are using on its classpath, or it will complain to standard error
            // about not having that. Since it doesn't inherit our full
            // classpath by default, we need to add the JAR we would load
            // org.slf4j.impl.StaticLoggerBinder from. This is complicated by
            // the fact that that class my not necessarily have actually been
            // provided by the application developer, so we can't just load it
            // and look where it came from.
            
            // Get the class object for the binder we want to send out, via
            // reflection. If it does not exist, an exception will be produced.
            val binderClass = Class.forName("org.slf4j.impl.StaticLoggerBinder")
            
            // Get the .jar it lives in. If it somehow doesn't live in a .jar,
            // this will error out. See
            // <http://stackoverflow.com/q/1983839/402891>
            val binderJar = binderClass.getProtectionDomain.getCodeSource
                .getLocation
            
            // Add the .jar to the Spark context
            sc.addJar(binderJar.toString)
            
        } catch {
            case e: Exception => {
                println("""Failed to find an slf4j StaticLoggerBinder .jar to 
                ship out with Spark jobs. Are you sure you have set up slf4j 
                correctly in your application?""".stripMargin)
                println(e.getMessage)
            }
        }
        
        // TODO: The above shipping out of the log binder .jar *should* work,
        // but does not appear to: the angry message still gets printed by the
        // worker. Figure out why this is and fix it.
                
        // Set up Parquet to write using Avro for this job
        ParquetOutputFormat.setWriteSupportClass(job, classOf[AvroWriteSupport])
        
        // Go get the static SCHEMA$ for the Avro record type. We can't do it
        // the normal way (RecordType.SCHEMA$) because Scala keeps static things
        // in companion objects with the same names as the types they really
        // belong to, and type parameters don't automatically bring them along.
        
        // First we need the class of the record, which we need anyway for doing
        // the writing.
        val recordClass = m.erasure
        
        println("Writing out a collection of %s".format(
            recordClass.getCanonicalName))
        
        // Get the schema Field object    
        val recordSchemaField = recordClass.getDeclaredField("SCHEMA$")
        // Get the static value of this field. Argument is ignored for static
        // fields, so pass null. See <http://docs.oracle.com/javase/7/docs/api/j
        // ava/lang/reflect/Field.html#get%28java.lang.Object%29>
        val recordSchema = recordSchemaField.get(null).asInstanceOf[Schema]
        
        // Set the Avro schema to use for this job
        AvroParquetOutputFormat.setSchema(job, recordSchema)
        
        // Make an RDD of serializeable-wrapped Avro records, with null keys.
        // Put it all into one partition. We have to make sure we make a seq of
        // the values and map on that; if we are not careful with the types, the
        // keys from the HashMap can end up being our RDD keys, and we for some
        // reason need to definitely have a PairRDD (i.e. an RDD of key, value
        // pairs) that definitely has null keys.
        val rdd = sc.makeRDD(things.toSeq map { r => 
            (null, r)
        }, 1)
        
        // Save the RDD to a Parquet file in our output directory. The keys are
        // void, the values are RecordTypes, the output format is a
        // ParquetOutputFormat for RecordTypes, and the configuration of the job
        // we set up is used.
        rdd.saveAsNewAPIHadoopFile(directory, classOf[Void], 
            recordClass, classOf[ParquetOutputFormat[RecordType]],
            job.getConfiguration)
                    
    }
    
    /**
     * Write Parquet Avro files with all of the parts of the graph to the
     * specified directory.
     */
    def writeParquetFiles(directory: String, sc: SparkContext) = {
        // Put the sides under "/Sides"
        writeCollectionToParquet(sides.values, directory + "/Sides", sc)
        // AlleleGroups go under /AlleleGroups
        writeCollectionToParquet(alleleGroups, directory + "/AlleleGroups", sc)
        // Adjacencies go under /Adjacencies
        writeCollectionToParquet(adjacencies, directory + "/Adjacencies", sc)
        // Anchors go under /Anchors
        writeCollectionToParquet(anchors, directory + "/Anchors", sc)
    }
    
    
}

/**
 * ParquetSequenceGraphBuilder: a class that builds a sequence graph into
 * Parquet files in the specified directory, using the given SparkContext to run
 * its jobs.
 *
 * TODO: Make streaming. For now just uses an InMemorySequenceGraphBuilder and
 * writes out from that.
 */
class ParquetSequenceGraphBuilder(sample: String, reference: String, 
    directory: String, sc: SparkContext) 
    extends InMemorySequenceGraphBuilder(sample, reference) {
    
    // Write out the Parquet files at the end
    override def finish() = writeParquetFiles(directory, sc)
    
}








