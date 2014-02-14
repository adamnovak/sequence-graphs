package edu.ucsc.genome

import org.scalatest._

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

class SequenceGraphTests extends FunSuite with SparkSuite {

    // Make a pen and draw a graph
    val pen = new SequenceGraphPen
    
    // Make a SequenceGraphChunk to colect drawn things
    var chunk = new SequenceGraphChunk

    // Draw a couple sides    
    val leftSide = pen.drawSide("chr1", 10, Face.LEFT)
    chunk += leftSide
    val rightSide = pen.drawSide("chr1", 19, Face.RIGHT)
    chunk += rightSide
    
    // Draw an edge
    chunk += pen.drawAlleleGroup(leftSide, rightSide,
        new Allele("AAAAAAAAA", "ACTGTCGGG"))
    
    var graph: SequenceGraph = null

    test("can be created") {
        // Set up the graph. sc is in scope since we are a SparkSuite
        graph = new SequenceGraph(sc.parallelize(List(chunk)))
    }
    
    test("contains the sides and edges that went into it") {
        val collectedSides = graph.sides.collect
        assert(collectedSides.size == 2)
        assert(collectedSides contains leftSide)
        assert(collectedSides contains rightSide)
        
        val collectedAlleleGroups = graph.alleleGroups.collect
        assert(collectedAlleleGroups.size == 1)
        assert(collectedAlleleGroups(0).allele.asInstanceOf[Allele]
            .bases == "AAAAAAAAA")
    }
    
    
    test("inRange can be executed") {
        graph.inRange(new BaseRange("chr10", 0, 100)).sides.collect
    }
    
    test("inRange excludes things on wrong contig") {
        assert(graph.inRange(new BaseRange("chr10", 0, 100)).sides.collect
            .size == 0)
    }
    
    test("inRange excludes things upstream") {
        assert(graph.inRange(new BaseRange("chr1", 0, 9)).sides.collect
            .size == 0)
    }
    
    test("inRange excludes things downstream") {
        assert(graph.inRange(new BaseRange("chr1", 20, 100)).sides.collect
            .size == 0)
    }
    
    test("inRange gets things overlapping on the left") {
    
        val subgraph = graph.inRange(new BaseRange("chr1", 0, 15))
        
        val collectedSides = subgraph.sides.collect
        assert(collectedSides.size == 2)
        assert(collectedSides contains leftSide)
        assert(collectedSides contains rightSide)
        
        val collectedAlleleGroups = subgraph.alleleGroups.collect
        assert(collectedAlleleGroups.size == 1)
        assert(collectedAlleleGroups(0).allele.asInstanceOf[Allele]
            .bases == "AAAAAAAAA")
    }
    
    test("inRange gets things overlapping on the right") {
    
        val subgraph = graph.inRange(new BaseRange("chr1", 16, 50))
        
        val collectedSides = subgraph.sides.collect
        assert(collectedSides.size == 2)
        assert(collectedSides contains leftSide)
        assert(collectedSides contains rightSide)
        
        val collectedAlleleGroups = subgraph.alleleGroups.collect
        assert(collectedAlleleGroups.size == 1)
        assert(collectedAlleleGroups(0).allele.asInstanceOf[Allele]
            .bases == "AAAAAAAAA")
    }
    
    test("inRange gets things overlapping in the middle") {
    
        val subgraph = graph.inRange(new BaseRange("chr1", 17, 17))
        
        val collectedSides = subgraph.sides.collect
        assert(collectedSides.size == 2)
        assert(collectedSides contains leftSide)
        assert(collectedSides contains rightSide)
        
        val collectedAlleleGroups = subgraph.alleleGroups.collect
        assert(collectedAlleleGroups.size == 1)
        assert(collectedAlleleGroups(0).allele.asInstanceOf[Allele]
            .bases == "AAAAAAAAA")
    }
    
    test("inRange gets things completly contained") {
    
        val subgraph = graph.inRange(new BaseRange("chr1", 0, 100))
        
        val collectedSides = subgraph.sides.collect
        assert(collectedSides.size == 2)
        assert(collectedSides contains leftSide)
        assert(collectedSides contains rightSide)
        
        val collectedAlleleGroups = subgraph.alleleGroups.collect
        assert(collectedAlleleGroups.size == 1)
        assert(collectedAlleleGroups(0).allele.asInstanceOf[Allele]
            .bases == "AAAAAAAAA")
    }
    
    
}
