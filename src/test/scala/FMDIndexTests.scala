package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._
import org.ga4gh.FMDIndex

/**
 * Tests for making and using FMDIndex objects.
 */
class FMDIndexTests extends FMDSuite {

    var index: FMDIndex = null;

    test("FMDIndex can be created") {
        index = new FMDIndex(basename)
    }
    
    
    
}
