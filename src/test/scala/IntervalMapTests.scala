package edu.ucsc.genome

import org.scalatest._

/**
 * Tests for the IntervalMap class, to make sure it doesn't screw up containment
 * or something.
 */
class IntervalMapTests extends FunSuite {

    var map: IntervalMap[String] = null

    test("can be created") {
        map = new IntervalMap[String]
    }
    
    test("can insert") {
        map((5, 10)) = "Hello"
        map((1, 4)) = "World"
    }
    
    test("can retrieve exactly") {
        assert(map.get(5, 10) === Some("Hello"))
        assert(map.get(1, 4) === Some("World"))
        
        assert(map(5, 10) === "Hello")
        assert(map(1, 4) === "World")
    }
    
    test("can retrieve by subinterval") {
        assert(map.get(6, 6) === Some("Hello"))
        assert(map.get(2, 4) === Some("World"))
        
        assert(map(7, 8) === "Hello")
        assert(map(1, 2) === "World")
    }
    
    test("ignores non-subintervals") {
        assert(map.get(0, 10) === None)
        assert(map.get(4, 6) === None)
        assert(map.get(-100, 100) === None)
    }


}
