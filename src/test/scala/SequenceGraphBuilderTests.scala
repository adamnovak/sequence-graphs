package edu.ucsc.genome

import org.scalatest._

class SequenceGraphBuilderTests extends FunSuite {

    val builder = new InMemorySequenceGraphBuilder("bob", "hg19")
    
    test("starts with no AlleleGroups") {
        assert(builder.getLastAlleleGroup("chr1", 0) === None)
        assert(builder.getLastAlleleGroup("chr1", 1) === None)
    }

    test("can create a leading telomere") {
        val telomere = builder.getLastSide("chr1", 0)
        assert(telomere.position.contig === "chr1")
        assert(telomere.position.base === 0)
        assert(telomere.position.face === Face.RIGHT)
    }
    
    test("gives each phase its own leading telomere") {
        val telomere0 = builder.getLastSide("chr1", 0)
        val telomere1 = builder.getLastSide("chr1", 1)
        
        assert(telomere0.id != telomere1.id)
        assert(telomere0.position === telomere1.position)
    }
    
    test("can add a phased pair of Alleles") {
        val allele0 = new Allele("A")
        val allele1 = new Allele("T")
    
        val startSide0 = builder.getLastSide("chr1", 0)
        val startSide1 = builder.getLastSide("chr1", 1)
        
        builder.addAllele("chr1", List(0), allele0)
        val midSide0 = builder.getLastSide("chr1", 0)
        val midSide1 = builder.getLastSide("chr1", 1)
        
        assert(midSide0 != startSide0)
        assert(midSide0 != midSide1)
        assert(midSide1 === startSide1)
        assert(builder.getLastAlleleGroup("chr1", 0).get.allele === allele0)
        
        builder.addAllele("chr1", List(1), allele1)
        val endSide0 = builder.getLastSide("chr1", 0)
        val endSide1 = builder.getLastSide("chr1", 1)
        
        assert(endSide0 != endSide1)
        assert(endSide1 != midSide1)
        assert(endSide0 === midSide0)
        assert(builder.getLastAlleleGroup("chr1", 0).get.allele === allele0)
        assert(builder.getLastAlleleGroup("chr1", 1).get.allele === allele1)
    }
    
    test("can add a phased pair of Anchors") {
        val startSide0 = builder.getLastSide("chr1", 0)
        val startSide1 = builder.getLastSide("chr1", 1)
        
        builder.addAnchor("chr1", List(0), 100)
        val midSide0 = builder.getLastSide("chr1", 0)
        val midSide1 = builder.getLastSide("chr1", 1)
        
        assert(midSide0 != startSide0)
        assert(midSide0 != midSide1)
        assert(midSide1 === startSide1)
        assert(builder.getLastAlleleGroup("chr1", 0) === None)
        
        builder.addAnchor("chr1", List(1), 100)
        val endSide0 = builder.getLastSide("chr1", 0)
        val endSide1 = builder.getLastSide("chr1", 1)
        
        assert(endSide0 != endSide1)
        assert(endSide1 != midSide1)
        assert(endSide0 === midSide0)
        assert(builder.getLastAlleleGroup("chr1", 0) === None)
        assert(builder.getLastAlleleGroup("chr1", 1) === None)
    }
    
    test("can add an unphased Allele") {
        val allele = new Allele("G")
        
        val startSide0 = builder.getLastSide("chr1", 0)
        val startSide1 = builder.getLastSide("chr1", 1)
        
        builder.addAllele("chr1", List(0, 1), allele)
        val endSide0 = builder.getLastSide("chr1", 0)
        val endSide1 = builder.getLastSide("chr1", 1)
        
        assert(endSide0 === endSide1)
        assert(endSide0 != startSide0)
        assert(endSide1 != startSide1)
        assert(builder.getLastAlleleGroup("chr1", 0).get.allele === allele)
        assert(builder.getLastAlleleGroup("chr1", 1).get.allele === allele)
    }
    
    test("can add an unphased Anchor") {
        val startSide0 = builder.getLastSide("chr1", 0)
        val startSide1 = builder.getLastSide("chr1", 1)
        
        builder.addAnchor("chr1", List(0, 1), 10)
        val endSide0 = builder.getLastSide("chr1", 0)
        val endSide1 = builder.getLastSide("chr1", 1)
        
        assert(endSide0 === endSide1)
        assert(endSide0 != startSide0)
        assert(endSide1 != startSide1)
        assert(builder.getLastAlleleGroup("chr1", 0) === None)
        assert(builder.getLastAlleleGroup("chr1", 1) === None)
    }
}
