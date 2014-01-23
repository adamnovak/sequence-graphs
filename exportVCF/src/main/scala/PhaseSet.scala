package edu.ucsc.genome.ExportVCF

class PhaseSet (val contig: String,
                val start: Long,
                val end: Long) {
}

class PhaseSetBuilder (contig: String) {

  var countLast = 1
  var endLast = 0L
  var startLast = 0L

  var sets = List[(Long, Long)]()

  def add (kv: (Int, Long, Long)) {
    add (kv._1, kv._2, kv._3)
  }

  def add (count: Int, start: Long, end: Long) {
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

  def toSets (): List[PhaseSet] = {
    sets.map(p => new PhaseSet(contig, p._1, p._2))
  }
}
