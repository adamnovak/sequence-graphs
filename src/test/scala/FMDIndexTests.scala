package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._

/**
 * Tests for making and using FMDIndex objects.
 */
class FMDIndexTests extends RLCSASuite {

    var index: FMDIndex = null;

    test("FMDIndex can be created") {
        index = new FMDIndex(basename)
    }
    
    test("pairToPosition works at start of forward strand") {
        val position = index.pairToPosition(0, 0)
        
        assert(position.contig === "seq1")
        assert(position.base === 1)
        assert(position.face === Face.LEFT)
    }
    
    test("pairToPosition works at start of a reverse strand") {
        val position = index.pairToPosition(1, 0)
        assert(position.contig === "seq1")
        assert(position.base === 10)
        assert(position.face === Face.RIGHT)
    }
    
    test("pairToPosition works at end of a forward strand") {
        val position = index.pairToPosition(0, 9)
        
        assert(position.contig === "seq1")
        assert(position.base === 10)
        assert(position.face === Face.LEFT)
    }
    
    test("pairToPosition works at end of a reverse strand") {
        val position = index.pairToPosition(1, 9)
        
        assert(position.contig === "seq1")
        assert(position.base === 1)
        assert(position.face === Face.RIGHT)
    }
    
    test("positionToBWT inverts bwtToPosition") {
        // BWT space is ((# of contigs) + (totoal contig length)) * 2 in size.
        // We know that RLCSASuite has 2 sequences, one of 10 bases and one of
        // 11.
        val bwtSize = (2 + (10 + 11)) * 2
        // The actual letters in the BWT start after the 2 * # of texts text end
        // characters
        val bwtStart = 2 * 2
        
        for(bwt <- bwtStart until bwtSize) {
            // For a bunch of positions in the BWT
            assert(index.positionToBWT(index.bwtToPosition(bwt)) === bwt)
        }
    }
    
    test("bwtToPosition works by index") {
        // Get the location of the first non-end character in the BWT. 0-3 are
        // $.
        val position = index.bwtToPosition(4)
        // The first suffix is going to be the one that starts AAG: all of seq2.
        assert(position.contig === "seq2")
        // Sequence position is 1-based.
        assert(position.base === 1)
        assert(position.face === Face.LEFT)
        
        // And the second one ought to be AAT: all of seq1
        val position2 = index.bwtToPosition(5)
        assert(position2.contig === "seq1")
        assert(position2.base === 1)
        assert(position2.face === Face.LEFT)
    }
    
    test("positionToBWT works correctly") {
        // Make sure we match the test above
        assert(index.positionToBWT(new Position("seq2", 1, Face.LEFT)) === 4)
        assert(index.positionToBWT(new Position("seq1", 1, Face.LEFT)) === 5)
    }
    
    test("bwtToPosition inverts positionToBWT") {
        // Let's try every base on seq1 (length 10), for both strands
        for(i <- 1 until 11) {
            // Check the left sides
            val pos = new Position("seq1", i, Face.LEFT)
            assert(index.bwtToPosition(index.positionToBWT(pos)) === pos)
            
            // And the right sides
            val pos2 = new Position("seq1", i, Face.RIGHT)
            assert(index.bwtToPosition(index.positionToBWT(pos2)) === pos2)
        }
    }
    
    test("does not find missing substring") {
        assert(index.count("AAAAA") == 0)
    }
    
    test("finds entire contig") {
        assert(index.count("AATCTACTGC") == 1)
    }
    
    test("finds reverse complement") {
        assert(index.count("GCTAGTAGCTT") == 1)
    }
    
    test("finds common subsequence") {
        assert(index.count("AA") == 2)
    }
    
    test("fails to left-map left-ambiguous thing") {
        val mapping = index.map("AATCTACTGC", 2, Face.LEFT)
        assert(mapping == None)
    }
    
    test("right-maps left-ambiguous thing") {
        val position = index.map("AATCTACTGC", 2, Face.RIGHT).get
        
        assert(position.contig === "seq1")
        assert(position.base === 2)
        assert(position.face === Face.RIGHT)
    }
    
    test("left-maps last base of contig") {
        val position = index.map("AATCTACTGC", 10, Face.LEFT).get
        
        assert(position.contig === "seq1")
        assert(position.base === 10)
        assert(position.face === Face.LEFT)
    }
    
    test("left-maps last base of reverse complement") {
        val position = index.map("GCTAGTAGCTT", 11, Face.LEFT).get
        
        assert(position.contig === "seq2")
        assert(position.base === 1)
        assert(position.face === Face.RIGHT)
    }
    
    test("left-maps unambiguous bases in contig") {
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Position]] = index.map(pattern, Face.LEFT)
        
        // The last 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 8)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left
                assert(mapping.face == Face.LEFT)
                // And on this contig
                assert(mapping.contig == "seq1")
            }
            case None => {}
        }
        
        
    }
    
    test("left-maps unambiguous bases in reverse complement") {
        val pattern = "GCAGTAGATT"
        val mappings: Seq[Option[Position]] = index.map(pattern, Face.LEFT)
        
        // The last 8 characters ought to map this way too.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 8)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left
                assert(mapping.face == Face.RIGHT)
                // And on this contig
                assert(mapping.contig == "seq1")
            }
            case None => {}
        }
        
        
    }
    
    test("right-maps unambiguous bases in contig") {
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Position]] = index.map(pattern, Face.RIGHT)
        
        // The first 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 8)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the right
                assert(mapping.face == Face.RIGHT)
                // And on this contig
                assert(mapping.contig == "seq1")
            }
            case None => {}
        }
        
        
    }
    
    test("left-maps correctly to 1-item ranges") {
        import fi.helsinki.cs.rlcsa.{RangeEncoder, RangeVector}
        
        // BWT space is ((# of contigs) + (totoal contig length)) * 2 in size.
        // We know that RLCSASuite has 2 sequences, one of 10 bases and one of
        // 11.
        val bwtSize = (2 + (10 + 11)) * 2
        
        // Make a RangeVector that says every base in BWT space is in a
        // different range. Use a block size of 32.
        val rangeEncoder = new RangeEncoder(32)
        // Set it to 1 everywhere
        rangeEncoder.addRun(0, bwtSize)
        rangeEncoder.flush()
        
        // Make a RangeVector from it
        val rangeVector = new RangeVector(rangeEncoder, bwtSize)
        
        // Now try mapping with the range vector
        val pattern = "AATCTACTGC"
        val mappings: Seq[Long] = index.map(rangeVector, pattern, Face.LEFT)
        
        println((pattern, mappings).zipped.mkString("\n"))
        
        // The last 8 characters ought to map
        assert(mappings(0) === -1)
        assert(mappings(1) === -1)
        
        assert(mappings.map {
            case -1 => 0
            case _ => 1
        }.sum === 8)
        
        // There should be 9 distint values: 8 mapping locations and two -1s
        assert(mappings.distinct.size === 9)
        
    }
    
    test("right-maps correctly to 1-item ranges") {
        import fi.helsinki.cs.rlcsa.{RangeEncoder, RangeVector}
        
        // BWT space is ((# of contigs) + (totoal contig length)) * 2 in size.
        // We know that RLCSASuite has 2 sequences, one of 10 bases and one of
        // 11.
        val bwtSize = (2 + (10 + 11)) * 2
        
        // Make a RangeVector that says every base in BWT space is in a
        // different range. Use a block size of 32.
        val rangeEncoder = new RangeEncoder(32)
        // Set it to 1 everywhere
        rangeEncoder.addRun(0, bwtSize)
        rangeEncoder.flush()
        
        // Make a RangeVector from it
        val rangeVector = new RangeVector(rangeEncoder, bwtSize)
        
        // Now try mapping with the range vector
        val pattern = "AATCTACTGC"
        val mappings: Seq[Long] = index.map(rangeVector, pattern, Face.RIGHT)
        
        // The first 8 characters ought to map
        assert(mappings.map {
            case -1 => 0
            case _ => 1
        }.sum === 8)
        
        assert(mappings(8) === -1)
        assert(mappings(9) === -1)
        
        // There should be 9 distint values: 8 mapping locations and two -1s
        assert(mappings.distinct.size === 9)
    }
    
    test("left-maps correctly to 1-item ranges on reverse complement") {
        import fi.helsinki.cs.rlcsa.{RangeEncoder, RangeVector}
        
        // BWT space is ((# of contigs) + (total contig length)) * 2 in size.
        // We know that RLCSASuite has 2 sequences, one of 10 bases and one of
        // 11.
        val bwtSize = (2 + (10 + 11)) * 2
        
        // Make a RangeVector that says every base in BWT space is in a
        // different range. Use a block size of 32.
        val rangeEncoder = new RangeEncoder(32)
        // Set it to 1 everywhere
        rangeEncoder.addRun(0, bwtSize)
        rangeEncoder.flush()
        
        // Make a RangeVector from it
        val rangeVector = new RangeVector(rangeEncoder, bwtSize)
        
        // TODO: Fixturize this range thing.
        
        // Now try mapping with the range vector
        val pattern = "GCAGTAGATT"
        val mappings: Seq[Long] = index.map(rangeVector, pattern, Face.LEFT)
        
        // The last 8 characters ought to map this way around also.
        assert(mappings(0) === -1)
        assert(mappings(1) === -1)
        
        assert(mappings.map {
            case -1 => 0
            case _ => 1
        }.sum === 8)
        
        // There should be 9 distint values: 8 mapping locations and two -1s
        assert(mappings.distinct.size === 9)
    }
    
}
