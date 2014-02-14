package edu.ucsc.genome

// We want to write files
import java.io._



/**
 * SequenceGraphWriter: something that provides a place to put sequence graph
 * components. Write methods must not be called multiple times on the same
 * object.
 */
trait SequenceGraphWriter {
    /**
     * Write the given AlleleGroup.
     */
    def writeAlleleGroup(alleleGroup: AlleleGroup) : Unit
    
    /**
     * Write the given Adjacency.
     */
    def writeAdjacency(adjacency: Adjacency) : Unit
    
    /**
     * Write the given Anchor.
     */
    def writeAnchor(anchor: Anchor) : Unit
    
    /**
     * Add the given Side to the graph. Called exactly once per Side.
     */
    def writeSide(side: Side) : Unit
    
    /**
     * Finish writing the sequence graph.
     */
    def close() : Unit
}

/**
 * GraphvizSequenceGraphWriter: save a sequence graph on disk, writing to a
 * GraphViz dot file in a streaming fashion.
 *
 * Write it to the given file.
 *
 */
class GraphvizSequenceGraphWriter(file: String) extends SequenceGraphWriter {
    
    // We need to write escaped strings. See
    // <http://stackoverflow.com/a/9914380/402891> and
    // <http://commons.apache.org/proper/commons-
    // lang/apidocs/org/apache/commons/lang3/StringEscapeUtils.html>
    import org.apache.commons.lang.StringEscapeUtils.escapeJava
    
    // Open the file for writing. See
    // <http://www.tutorialspoint.com/scala/scala_file_io.htm>
    val graphWriter = new java.io.PrintWriter(new File(file))
    
    // Write a header
    graphWriter.write("digraph sample {\n")
    
    // Set up edge styles
    graphWriter.write("edge [arrowsize=\"0\"];\n")
    
    // And node styles
    graphWriter.write("node [shape=\"point\"];\n")
    
    // Implementations of what we need to be an EasySequenceGraphBuilder
    
    def writeAlleleGroup(alleleGroup: AlleleGroup) : Unit = {
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
    
    def writeAdjacency(adjacency: Adjacency) : Unit = {
        // Make an edge for every Adjacency
        graphWriter.write("%d -> %d;\n".format(
            adjacency.edge.left, adjacency.edge.right))
    }
    
    def writeAnchor(anchor: Anchor) : Unit ={
        // What should we write on the anchor?
        val label = "Anchor x%d".format(anchor.ploidy.lower)
        
        // Make an edge for every Anchor
        graphWriter.write("%d -> %d [label=\"%s\",color=\"#0000ff\"];\n"
            .format(anchor.edge.left, anchor.edge.right, escapeJava(label)))
    }
    
    def writeSide(side: Side) : Unit = {
        // What should we write on the Side?
        val label = "%s:%d-%s".format(side.position.contig, side.position.base, 
            side.position.face match {
                case Face.LEFT => "L"
                case Face.RIGHT => "R"
            })
            
        // Getting the label on the point is tricky. See 
        // <http://marc.info/?l=graphviz-interest&m=126029250410816>
            
        // Write out a node to label the side (because point nodes can't have
        // real labels)
        graphWriter.write("L%d [shape=\"plaintext\",label=\"%s\"];\n"
            .format(side.id, escapeJava(label)))

        // Write out the Side without a label
        graphWriter.write("%d;\n"
            .format(side.id))
            
        // Connect them with an invisible edge
        graphWriter.write(
            "{rank=same %d -> L%d [dir=\"none\",style=\"dotted\"];}\n"
            .format(side.id, side.id))
    }
    
    def close() {
        // Close the graph block
        graphWriter.write("}\n")
        
        // Close the file
        graphWriter.close()
        
        println("Wrote %s".format(file))
    }
}
