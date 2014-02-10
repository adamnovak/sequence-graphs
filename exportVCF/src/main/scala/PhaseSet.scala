package edu.ucsc.genome.ExportVCF

/**
 * A global phase set. I.e., a set of called variants that exist on a single phased
 * haplotype.
 */
case class PhaseSet (contig: String, start: Long, end: Long) {
}

/**
 * Builder for a phase set.
 */
class PhaseSetBuilder (contig: String) {

  // "state machine" signals
  var countLast = 1
  var endLast = 0L
  var startLast = 0L

  var sets = List[(Long, Long)]()

  /**
   * Adder for a tuple pair.
   *
   * @param kv Edge set tuple pair to add. Tuple is (count, (start, end)).
   */
  def add (kv: (Int, Long, Long)) {
    add (kv._1, kv._2, kv._3)
  }

  /**
   * Add method for a new position pair.
   *
   * @param count Number of anchors at this position.
   * @param start Start position of anchor.
   * @param end End position of anchor.
   */
  def add (count: Int, start: Long, end: Long) {
    // pseudo-"fsm" for adding phase sets:
    // - look at current number of anchors vs. last number of anchors
    // - if we see a set with only one anchor, that set is not phased
    // - if we see a set with more than one anchor, that set is phased
    // whenever we see a change from:
    // - unphased to phased --> we add a new phase set
    // - phased to unphased --> we close the current phase set
    if (count == 1 && countLast != 1) {
      // end last seen phase set
      sets = (startLast, start) :: sets
      countLast = count
    } else if (count != 1 && countLast == 1) {
      // log start of new phase set
      startLast = start
      countLast = count
    }

    endLast = end
  }

  /**
   * From this builder, returns a list of phase sets.
   *
   * @param return A list of phase sets.
   */
  def toSets (): List[PhaseSet] = {
    sets.map(p => new PhaseSet(contig, p._1, p._2))
  }
}
