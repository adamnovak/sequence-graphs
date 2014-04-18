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
    
    test("positionToContigNumber works") {
        // Ends of contig 0
        assert(index.positionToContigNumber(0) === 0)
        assert(index.positionToContigNumber(9) === 0)
        // Ends of contig 1
        assert(index.positionToContigNumber(10) === 1)
        assert(index.positionToContigNumber(20) === 1)
    }
    
    test("contigNumberToPosition works") {
        // Should return first position in each contig.
        assert(index.contigNumberToPosition(0) === 0)
        assert(index.contigNumberToPosition(1) === 10)
    }
    
    test("pairToSide works at start of forward strand") {
        val side = index.pairToSide(0, 0)
        
        assert(side.coordinate === 0)
        assert(side.face === Face.LEFT)
    }
    
    test("pairToSide works at start of a reverse strand") {
        val side = index.pairToSide(1, 0)
        assert(side.coordinate === 9)
        assert(side.face === Face.RIGHT)
    }
    
    test("pairToSide works at end of a forward strand") {
        val side = index.pairToSide(0, 9)
        
        assert(side.coordinate === 9)
        assert(side.face === Face.LEFT)
    }
    
    test("pairToSide works at end of a reverse strand") {
        val side = index.pairToSide(1, 9)
        
        assert(side.coordinate === 0)
        assert(side.face === Face.RIGHT)
    }
    
    test("sideToPair inverts pairToSide") {
        val pairs = Seq((0, 0), (1, 0), (0, 9), (1, 9), (2, 0), (2, 10), (3, 0), 
            (3, 10))
        
        for(pair <- pairs) {
            // Convert to a pair and back.
            val side = index.pairToSide(pair._1, pair._2)
            val newPair = index.sideToPair(side)
            
            assert(newPair === pair)
        }
    }
    
    test("sideToBWT inverts bwtToSide") {
        // BWT space is ((# of contigs) + (totoal contig length)) * 2 in size.
        // We know that RLCSASuite has 2 sequences, one of 10 bases and one of
        // 11.
        val bwtSize = (2 + (10 + 11)) * 2
        // The actual letters in the BWT start after the 2 * # of texts text end
        // characters
        val bwtStart = 2 * 2
        
        for(bwt <- bwtStart until bwtSize) {
            // For a bunch of sides in the BWT
            assert(index.sideToBWT(index.bwtToSide(bwt)) === bwt)
        }
    }
    
    test("bwtToSide works by index") {
        // Get the location of the first non-end character in the BWT. 0-3 are
        // $.
        val side = index.bwtToSide(4)
        // The first suffix is going to be the one that starts AAG: all of seq2.
        assert(side.coordinate === 10)
        assert(side.face === Face.LEFT)
        
        // And the second one ought to be AAT: all of seq1
        val side2 = index.bwtToSide(5)
        assert(side2.coordinate === 0)
        assert(side2.face === Face.LEFT)
    }
    
    test("sideToBWT works correctly") {
        // Make sure we match the test above
        
        // Start of seq 2
        assert(index.sideToBWT(new Side(10, Face.LEFT)) === 4)
        // Start of seq 1
        assert(index.sideToBWT(new Side(0, Face.LEFT)) === 5)
    }
    
    test("bwtToSide inverts sideToBWT") {
        // Let's try every base on seq1 and seq2, for both strands
        for(i <- 0 until 21) {
            // Check the left sides
            // Starting at 0: first ID in seq1
            val pos = new Side(i, Face.LEFT)
            
            assert(index.bwtToSide(index.sideToBWT(pos)) === pos)
            
            // And the right sides
            val pos2 = new Side(i, Face.RIGHT)
            assert(index.bwtToSide(index.sideToBWT(pos2)) === pos2)
        }
    }
    
    test("displays characters correctly") {
        val seq1 = "AATCTACTGC"
        
        val got1 = seq1.zipWithIndex.map { case (char, offset) =>
            // Display the appropriate base and see what we get.
            val displayed = index.display(new Side(
                index.contigNumberToPosition(0) + offset, Face.LEFT))
            
            displayed
        }.mkString
        
        val seq2 = "AAGCTACTAGC"
        
        val got2 = seq2.zipWithIndex.map { case (char, offset) =>
            // Display the appropriate base and see what we get.
            val displayed = index.display(new Side(
                index.contigNumberToPosition(1) + offset, Face.LEFT))
            
            displayed
        }.mkString
        
        assert(seq1 === got1)
        assert(seq2 === got2)
        
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
        val side = index.map("AATCTACTGC", 2, Face.RIGHT).get
        
        assert(side.coordinate === index.contigNumberToPosition(0) + 1)
        assert(side.face === Face.RIGHT)
    }
    
    test("left-maps last base of contig") {
        val side = index.map("AATCTACTGC", 10, Face.LEFT).get
        
        assert(side.coordinate === index.contigNumberToPosition(0) + 9)
        assert(side.face === Face.LEFT)
    }
    
    test("left-maps last base of reverse complement") {
        val side = index.map("GCTAGTAGCTT", 11, Face.LEFT).get
        
        assert(side.coordinate === index.contigNumberToPosition(1))
        assert(side.face === Face.RIGHT)
    }
    
    test("left-maps unambiguous bases in contig") {
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Side]] = index.map(pattern, Face.LEFT)
        
        // The last 8 characters ought to map.
        assert(mappings.map {
            case Some(_) => 1
            case None => 0
        }.sum == 8)
        
        mappings.foreach {
            case Some(mapping) => {
                // Everything mapped should be mapped on the left
                assert(mapping.face === Face.LEFT)
                // And on this contig
                assert(index.positionToContigNumber(mapping.coordinate) === 0)
            }
            case None => {}
        }
        
        
    }
    
    test("left-maps unambiguous bases in reverse complement") {
        val pattern = "GCAGTAGATT"
        val mappings: Seq[Option[Side]] = index.map(pattern, Face.LEFT)
        
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
                assert(index.positionToContigNumber(mapping.coordinate) === 0)
            }
            case None => {}
        }
        
        
    }
    
    test("right-maps unambiguous bases in contig") {
        val pattern = "AATCTACTGC"
        val mappings: Seq[Option[Side]] = index.map(pattern, Face.RIGHT)
        
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
                assert(index.positionToContigNumber(mapping.coordinate) === 0)
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
