package edu.ucsc.genome

import org.scalatest._
import scala.collection.JavaConversions._

import org.apache.avro.ipc.{NettyServer, NettyTransceiver}
import org.apache.avro.ipc.specific.{SpecificRequestor, SpecificResponder}
import java.net.InetSocketAddress

/**
 * Tests for the SequenceGraphServer. Basically a translation and testification
 * of <http://gbif.blogspot.com/2011/06/getting-started-with-avro-rpc.html>.
 */
class SequenceGraphServerTests extends FunSuite with SparkSuite {

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
    
    // Make a SequenceGraph. sc is in scope since we are a SparkSuite
    val graph = new SequenceGraph(sc.parallelize(List(chunk)))

    // Have a place to store the server
    var server: NettyServer = null
    
    // And a place to put the client transciever
    var transciever: NettyTransceiver = null
    
    // And a place to put the client proxy
    var client: SequenceGraphAPI = null

    test("can be created") {
        // Set up the server to serve our graph
        server = new NettyServer(
            new SpecificResponder(classOf[SequenceGraphAPI],
            new SequenceGraphServer(graph)), new InetSocketAddress(7001))
    }
    
    test("can have a client") {
        // Set up the client
        
        // First the connection
        transciever = new NettyTransceiver(new InetSocketAddress(
            server.getPort))
        
        // Then the client proxy itself
        client = SpecificRequestor.getClient(classOf[SequenceGraphAPI],
            transciever)
    }
    
    test("can look up an empty range") {
        val page = client.search("chr1", 20, 25)
        
        assert(page.total == 0)
    }
    
    test("can look up a full range and get all pages") {
        var page = client.search("chr1", 11, 15)
        
        assert(page.total == 3)
        
        var allParts: List[GraphElement] = Nil ++ page.elements
        
        while(page.total > 0) {
            // Download all the pages until we fall off the end. TODO: have a
            // better way to signal this that isn't an empty page meaning the
            // end.
            page = client.getNextPage(page.nextPageToken)
            
            allParts ++= page.elements
        }
        
        // Now we downloaded the whole thing. Make sure it has what it should.
        assert(allParts.size == 3)
        assert(allParts.flatMap(_.element match {
            case side: Side => Some(side)
            case _ => None
        }).size == 2)
        assert(allParts.flatMap(_.element match {
            case alleleGroup: AlleleGroup => Some(alleleGroup)
            case _ => None
        }).size == 1)
    
    }
    
    test("can shut down") {
        transciever.close();
        server.close();
    }

}
 
      
