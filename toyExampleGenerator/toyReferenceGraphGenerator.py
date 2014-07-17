import sys
from sonLib.bioio import addNodeToGraph, addEdgeToGraph, setupGraphFile, finishGraphFile
from optparse import OptionParser
import random

"""Script to calculate and display reference genome hierarchies and mappings
"""

class BasePosition:
    def __init__(self, id, base):
        self.base = base
        self.id = id
        
    def getDotNodeName(self):
        return "n%sn" % self.id

class Side:
    def __init__(self, basePosition, orientation):
        self.adjacencies = set()
        self.otherSide = None
        self.basePosition = basePosition
        self.orientation = orientation #A True orientation means on the left, otherwise on the right, by convention
        self.mappedSides = [ self ]
       
    @staticmethod 
    def getReverseComplement(string):
        def fn(i):
            m = { 'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            if i in m:
                return m[i]
            return i
        return "".join([ fn(i) for i in string[::-1] ])
    
    def nonNAdjacencies(self):
        return [ adjSide for adjSide, ns in list(self.adjacencies) if ns == '' ]
    
    def base(self):
        if self.orientation:
            return self.basePosition.base
        return self.getReverseComplement(self.basePosition.base)
    
    def enumerateThreads(self, fn):
        stack = [ (side, side.otherSide.base()) for side in self.mappedSides ] #Add preceding base
        while len(stack) > 0:
            side, prefix = stack.pop()
            i = 0
            while len(side.nonNAdjacencies()) == 1 and i < 10:
                adjSide = side.nonNAdjacencies()[0]
                prefix += adjSide.base()
                side = adjSide.otherSide
                i += 1
            if fn(prefix) or len(side.nonNAdjacencies()) == 0:
                continue
            #assert len(side.nonNAdjacencies()) > 1
            for adjSide in side.nonNAdjacencies():
                stack.append((adjSide.otherSide, prefix + adjSide.base()))  

def hamDist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

def sharedPrefix(str1, str2):
    """Return shared prefix of the two strings"""
    for i in xrange(min(len(str1), len(str2))):
        if str1[i] != str2[i]:
            return str1[:i]
    return str1[:min(len(str1), len(str2))]

class ContextSet:
    def __init__(self):
        self.minimalUniqueStrings = set()
    
    def addString(self, string):
        """Add a string to the context set"""
        for i in xrange(1, len(string)):
            prefix = string[:i]
            if prefix in self.minimalUniqueStrings:
                return
        self.minimalUniqueStrings.add(string)

    def prefixInContextSet(self, string):
        for uniqueString in self.minimalUniqueStrings:
            if len(uniqueString) <= len(string) and hamDist(string, uniqueString) <= 0:
                return True
        return False  
    
    def maxLength(self):
        return max([ 0 ] + [ len(i) for i in self.minimalUniqueStrings ])
    
    def maxSharedPrefixLength(self, otherContextSet):
        i = ""
        for uniqueString in self.minimalUniqueStrings:
            for otherUniqueString in otherContextSet.minimalUniqueStrings:
                j = sharedPrefix(uniqueString, otherUniqueString)
                if len(j) > len(i):
                    i = j
        return len(i) 
    
    def tooShort(self, string):
        for uniqueString in self.minimalUniqueStrings:
            if len(uniqueString) > len(string):
                return True
        return False  
        
        
class SequenceGraph:
    def __init__(self, usePhasedContexts, label=""):
        self.id = 0
        self.sides = []
        self.contextSets = {}
        self.usePhasedContexts = usePhasedContexts
        if self.usePhasedContexts:
            self.mappedSequenceGraph = SequenceGraph(usePhasedContexts=False)
        else:
            self.mappedSequenceGraph = self
        self.label = label
            
    def getUniquePrefix(self, string, mismatches, dissimilarity=1):
        assert dissimilarity >= 1
        startPoints = sum([ [ (side, mappedSide, 0)  for mappedSide in side.mappedSides ] for side in self.sides ], [])
        index = 0
        
        while len(startPoints) > 0 and index < len(string):
            l = []
            for side, mappedSide, diff in startPoints:
                diff += mappedSide.base() != string[index]
                if diff <= mismatches + dissimilarity - 1:
                    l.append((side, mappedSide, diff))
            index += 1
            
            s2 = set([ side for side, mappedSide, diff in l if diff <= mismatches + dissimilarity - 1 ]) 
            s = set([ side for side, mappedSide, diff in l if diff <= mismatches ]) 
            if len(s2) == 1 and len(s) == 1: #We have a unique match to one node
                return string[:index], s.pop().otherSide
               
            startPoints = []
            for side, mappedSide, diff in l:
                for adjMappedSide in mappedSide.otherSide.nonNAdjacencies():
                    startPoints.append((side, adjMappedSide, diff)) 
                    
        return (None, None)
        
    def computeContextSets(self, minContextLength, maxContextLength, dissimilarity=1):
        """This function should be called on a graph before matching or printing is done"""
        self.contextSets = {}
        
        for side in self.sides:
            #Enumerate the threads
            contextSet = ContextSet()
            def fn(string):
                uniquePrefix, otherSide = self.getUniquePrefix(string, 0, dissimilarity=dissimilarity)
                if uniquePrefix == None:
                    return len(string) > maxContextLength #
                assert otherSide == side
                if len(uniquePrefix) < minContextLength:
                    if len(string) >= minContextLength:
                        contextSet.addString(string[:minContextLength])
                    else:
                        return False
                else:
                    contextSet.addString(uniquePrefix)
                return True #Return true stops the traversal on that thread
            side.enumerateThreads(fn)
            self.contextSets[side] = contextSet
    
    def computeContextSetsForMmers(self, m):
        """Computes contexts of upstream m-mers. """
        self.contextSets = {}
        for side in self.sides:
            #Enumerate the threads
            contextSet = ContextSet()
            def fn(string):
                if len(string) >= m:
                    contextSet.addString(string[:m])
                    return True
                return False 
            side.enumerateThreads(fn)
            self.contextSets[side] = contextSet
        
    def getMatches(self, side):
        #assert side not in self.sides
        matches = set()
        def fn(string):
            longEnough = True
            for otherSide in self.sides:
                if self.contextSets[otherSide].prefixInContextSet(string):
                    matches.add(otherSide)
                if self.contextSets[otherSide].tooShort(string):
                    longEnough = False
            return longEnough
        side.enumerateThreads(fn)
        return list(matches)
    
    def merge(self, side1, side2):
        if side1 == side2:
            return
        assert side1 != side2.otherSide
        assert side1.otherSide != side2
        self.sides.remove(side2)
        self.sides.remove(side2.otherSide)
        def fn(side1, side2):
            for adjSide, ns in list(side2.adjacencies):
                assert (side2, ns) in adjSide.adjacencies
                adjSide.adjacencies.remove((side2, ns))
                assert (side2, ns) not in adjSide.adjacencies
                if adjSide == side2: #Self loop on one side
                    adjSide = side1
                elif adjSide == side2.otherSide: #Self loop connecting sides
                    adjSide = side1.otherSide
                adjSide.adjacencies.add((side1, ns))
                side1.adjacencies.add((adjSide, ns))
        fn(side1, side2)
        fn(side1.otherSide, side2.otherSide)
        if side1.basePosition.id > side2.basePosition.id:
            side1.basePosition.id = side2.basePosition.id
        
        if self.usePhasedContexts:
            side1.mappedSides = side1.mappedSides + side2.mappedSides
            side1.otherSide.mappedSides = side1.otherSide.mappedSides + side2.otherSide.mappedSides
        
        for side in self.sides:
            for adjSide, ns in side.adjacencies:
                assert side != side2
                assert side != side2.otherSide
                assert adjSide != side2
                assert adjSide != side2.otherSide
                
    def positiveSides(self):
        return [ side for side in self.sides if side.orientation ]
    
    def mergeSequenceGraphs(self, sG2):
        sG2.renumber(startID = len(self.positiveSides()))
        self.sides += sG2.sides
        if self.usePhasedContexts:
            assert sG2.usePhasedContexts
            self.mappedSequenceGraph.mergeSequenceGraphs(sG2.mappedSequenceGraph)
        
    def addString(self, string):
        pSide = None
        phasedPSide = None
        firstSide = None
        circular = string[-1:] == '!'
        if circular:
            string = string[:-1]
        def tokenise(string):
            string = string.upper()
            ns = ""
            for i in string:
                if i != 'N':
                    yield i, ns
                    ns = ""
                else:
                    ns += 'N'
        for base, ns in tokenise(string):            
            def makeLinkedSides(base, ns, pSide, sides):
                bP = BasePosition(len(sides)/2, base)
                leftSide = Side(bP, 1)
                rightSide = Side(bP, 0)
                leftSide.otherSide = rightSide
                rightSide.otherSide = leftSide
                if pSide != None:
                    leftSide.adjacencies.add((pSide, ns))
                    pSide.adjacencies.add((leftSide, ns))
                sides.append(leftSide)
                sides.append(rightSide)
                return leftSide, rightSide
                    
            leftSide, rightSide = makeLinkedSides(base, ns, pSide, self.sides)
            if firstSide == None:
                firstSide = leftSide
            pSide = rightSide
    
            if self.usePhasedContexts:
                mappedLeftSide, mappedRightSide = makeLinkedSides(base, ns, phasedPSide, self.mappedSequenceGraph.sides)
                leftSide.mappedSides = [ mappedLeftSide ]
                rightSide.mappedSides = [ mappedRightSide ]
                phasedPSide = mappedRightSide
        
        if circular: #Need to add an extra adjacency
            firstSide.adjacencies.add((pSide, ''))
            pSide.adjacencies.add((firstSide, ''))
            if self.usePhasedContexts:
                firstSide.mappedSides[0].adjacencies.add((phasedPSide, ''))
                phasedPSide.adjacencies.add((firstSide.mappedSides[0], ''))
    
    def renumber(self, startID=0, prefix=""):
        ids = [ int(side.basePosition.id) for side in self.positiveSides() ]
        assert len(ids) == len(set(ids))
        ids.sort()
        for side in self.positiveSides():
            i = int(side.basePosition.id)
            assert i in ids
            i = ids.index(i) + startID
            if prefix != "":
                side.basePosition.id = str(prefix) + "_" + str(i)
            else:
                side.basePosition.id = i
        if self.usePhasedContexts:
            self.mappedSequenceGraph.renumber(startID=startID, prefix=prefix)
    
    def printDotFile(self, graphVizFileHandle, showContextSets, showIDs, number, displayAsSubgraph):
        if displayAsSubgraph:
            graphVizFileHandle.write('subgraph cluster_%s {\nstyle=filled;\ncolor=lightgrey;label = "%s";\nsplines=false;\n' % (number,self.label))
        else:
            graphVizFileHandle.write('subgraph cluster_%s {\nlabel = "%s";\n' % (number,self.label))
        #Add nodes
        for side in self.positiveSides():
            #Graph vis
            if showContextSets:
                leftContextString = ",".join([ Side.getReverseComplement(i[1:]) for i in self.contextSets[side].minimalUniqueStrings ])
                if len(self.contextSets[side].minimalUniqueStrings) == 0:
                    leftContextString = "None"
                rightContextString = ",".join([ i[1:] for i in self.contextSets[side.otherSide].minimalUniqueStrings ])
                if len(self.contextSets[side.otherSide].minimalUniqueStrings) == 0:
                    rightContextString = "None"
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, 
                               shape="record", label="ID=%s | L=%s | %s | R=%s" % (side.basePosition.id, 
                                                                                   leftContextString, side.basePosition.base, rightContextString))
            elif showIDs:
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, shape="record", label="ID=%s | %s" % (side.basePosition.id, side.basePosition.base))
            else:
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, shape="record", label="%s" % (side.basePosition.base))
        #Add edges
        seen = set()
        for side in self.sides:
            for adjSide, ns in side.adjacencies:
                if not (adjSide, side, ns) in seen:
                    assert (side, adjSide, ns) not in seen
                    def arrowShape(side):
                        if side.orientation:
                            return "normal"
                        return "crow"
                    addEdgeToGraph(parentNodeName=side.basePosition.getDotNodeName(), 
                                   childNodeName=adjSide.basePosition.getDotNodeName(), 
                                   graphFileHandle=graphVizFileHandle, colour="black", #weight="1", 
                                   dir="both, arrowtail=%s, arrowhead=%s" % 
                                   (arrowShape(side), arrowShape(adjSide)), style="solid", length="10", label=ns)
                    seen.add((side, adjSide, ns))
        #if displayAsSubgraph:
        graphVizFileHandle.write("}\n")
        
    

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    usage = "usage: %prog [options] <query> <target> <output>\n\n" + \
            "    <genome>:  strings representing genome\n" + \
            "    <string>: string to match\n" 
    description = "Script to demonstrate context sets and unique matching"
    parser = OptionParser(usage=usage, description=description)

    #output stuff
    parser.add_option("--graphVizFile", dest="graphVizFile", type="string",
                     help="File to do dump graphViz representation in (dot format)",
                     default="/dev/null")
    
    parser.add_option("--showContextSets", dest="showContextSets", type="string",
                     help="Show the context sets for the selected graphs (enumerated starting from 0)",
                     default="")
    
    parser.add_option("--showIDs", dest="showIDs", type="string",
                     help="Show the IDs for the selected graphs (enumerated starting from 0)",
                     default="")
    
    parser.add_option("--mergeContigs", dest="mergeContigs", type="string",
                     help="Merge the contig in the selected graphs (enumerated starting from 0)",
                     default="")
    
    #This feature is not working
    parser.add_option("--mismatches", dest="mismatches", type="int", 
                     help="Number of mismatches to allow between context strings",
                     default=0)
    
    parser.add_option("--dissimilarity", dest="dissimilarity", type="int", 
                     help="Number of differences to require to be defined as unique (e.g. dissimilarity=2 means a context string must be a hamming distance of 2 from all other context strings). Default=1",
                     default=1)
    
    parser.add_option("--minContextLength", dest="minContextLength", type="int", 
                     help="Minimum length of a string in a context set",
                     default=0)
    
    parser.add_option("--maxContextLength", dest="maxContextLength", type="int", 
                     help="Maximum length of a string in a context set",
                     default=10)
    
    parser.add_option("--showMultiMaps", dest="showMultiMaps", action="store_true",
                     help="Show doubly mapped match edges",
                     default=False)
    
    parser.add_option("--showOnlyLowestMaps", dest="showOnlyLowestMaps", action="store_true",
                     help="Show maps to one level in the hierarchy",
                     default=False)
    
    parser.add_option("--usePhasedContexts", dest="usePhasedContexts", action="store_true",
                     help="For each graph only allow context sets to be defined by the underlying sequences that serve as input to construct the sequence graph",
                     default=False)
    
    parser.add_option("--mergeForM", dest="mergeForM", type="string",
                     help="For each sequence graph indexed collapse to given m, requiring matching only on one side",
                     default="")
    
    parser.add_option("--mergeSymmetric", dest="mergeSymmetric", action="store_true",
                     help="When merging for m require mapping on both sides",
                     default=False)
    
    parser.add_option("--mapSymmetric", dest="mapSymmetric", action="store_true",
                     help="When mapping require matching on both sides",
                     default=False)
    
    parser.add_option("--targetSequenceGraphs", dest="targetSequenceGraphs", type="string",
                     help="List of sequence graph indices to map to",
                     default="")
    
    parser.add_option("--noMapping", dest="mapping", action="store_false",
                     help="Do not show mapping",
                     default=True)
    
    options, args = parser.parse_args()
    
    if len(args) == 0:
        parser.print_help()
        return 1
    
    mergeContigs = [ int(i) for i in options.mergeContigs.split() ]  
    print "We got merge-contigs", mergeContigs
    mergeForM = {} 
    for i in options.mergeForM.split():
        mergeForM[int(i.split("=")[0])] = int(i.split("=")[1])
    
    #First create the sequence graphs for each input graph
    sequenceGraphs = []
    
    for index in xrange(0, len(args), 2):
        assembly = args[index]
        label = args[index + 1]
        seqGraphIndex = index/2
        print "Processing sequence graph", index, options.mismatches, label
        sG = SequenceGraph(options.usePhasedContexts, label=label)
        sequenceGraphs.append(sG)
        
        matches = []
        m = {}
        if seqGraphIndex in mergeContigs:
            #First identify merges.
            sGs = []
            for string in assembly.split():
                print "Adding string:", string, " of assembly:", index
                sG2 = SequenceGraph(options.usePhasedContexts)
                sG2.addString(string)
                sG2.computeContextSets(minContextLength=options.minContextLength, maxContextLength=options.maxContextLength, 
                                      dissimilarity=options.dissimilarity)
                for sG3 in sGs:
                    for side in sG2.positiveSides():
                        leftMatches = sG3.getMatches(side)
                        rightMatches = sG3.getMatches(side.otherSide)
                        if len(leftMatches) == 1:
                            if rightMatches == [] or (leftMatches[0].otherSide == rightMatches[0] and len(rightMatches) == 1):
                                matches.append((leftMatches[0], side))
                        elif len(rightMatches) == 1 and len(leftMatches) == 0:
                            matches.append((rightMatches[0].otherSide, side))
                sGs.append(sG2)
            for sG2 in sGs:
                sG.mergeSequenceGraphs(sG2)
                print "Graph now has %i nodes" % len(sG.sides)
        else:
            for string in assembly.split():
                sG.addString(string)
            print "Graph now has %i nodes" % len(sG.sides)
        
        for side in sG.positiveSides():
            m[side] = side
            m[side.otherSide] = side.otherSide
               
        if seqGraphIndex in mergeForM:
            sG.computeContextSetsForMmers(m=mergeForM[seqGraphIndex])
            for side in sG.positiveSides():
                leftMatches = sG.getMatches(side)
                rightMatches = sG.getMatches(side.otherSide)
                if options.mergeSymmetric:
                    for matchingSide in leftMatches:
                        if matchingSide.otherSide in rightMatches:
                            matches.append((matchingSide, side))
                else:
                    for matchingSide in leftMatches + [ rightMatch.otherSide for rightMatch in rightMatches ]:
                        matches.append((matchingSide, side))
            
        for targetSide, inputSide in matches:
            while targetSide != m[targetSide]:
                targetSide = m[targetSide]
            print inputSide, m[inputSide]
            while inputSide != m[inputSide]:
                inputSide = m[inputSide]
            assert targetSide in sG.sides
            assert inputSide in sG.sides
            m[inputSide] = targetSide
            m[inputSide.otherSide] = targetSide.otherSide
            sG.merge(targetSide, inputSide)
        sG.renumber()   
        if seqGraphIndex in mergeForM:
            sG.computeContextSetsForMmers(m=mergeForM[seqGraphIndex])
        else:
            sG.computeContextSets(minContextLength=options.minContextLength, maxContextLength=options.maxContextLength, 
                                  dissimilarity=options.dissimilarity)
     
    #Now reindex them and print them
    showContextSets = [ int(i) for i in options.showContextSets.split() ]  
    showIDs = [ int(i) for i in options.showIDs.split() ]  
    graphVizFileHandle = open(options.graphVizFile, 'w')      
    setupGraphFile(graphVizFileHandle)
    graphVizFileHandle.write("splines=false;\n")   
    graphVizFileHandle.write("rankdir=LR;\n") 
       
    i = 1
    for index in xrange(len(sequenceGraphs)):
        sG = sequenceGraphs[index]
        sG.renumber(startID=i)
        i += len(sG.positiveSides())
        print "Renumbered graph %i with %i sides" % (index, len(sG.positiveSides()))
        #i += len(sG.positiveSides())
        sG.printDotFile(graphVizFileHandle, index in showContextSets and options.mismatches == 0, index in showIDs, index, len(sequenceGraphs) > 1)
        
    #Now print the matching edges between the graphs
    if options.mapping:
        targetSequenceGraphs = [ int(i) for i in options.targetSequenceGraphs.split() ]  
        if len(targetSequenceGraphs) == 0:
            targetSequenceGraphs = range(len(sequenceGraphs))
        for index in xrange(1, len(sequenceGraphs)):
            sGInput = sequenceGraphs[index]
            for side in sGInput.positiveSides():
                haveMatched = False
                
                for pIndex in xrange(index-1, -1, -1):
                    if pIndex in targetSequenceGraphs:
                        sGTarget = sequenceGraphs[pIndex]
                        
                        leftMatches = sGTarget.getMatches(side)
                        rightMatches = sGTarget.getMatches(side.otherSide)
                        
                        def addMatchEdge(colour, label, matchingSide):
                            if not options.showOnlyLowestMaps or not haveMatched:
                                addEdgeToGraph(parentNodeName=matchingSide.basePosition.getDotNodeName(), 
                                               childNodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, colour=colour, 
                                               weight="100", label=label, dir="both, arrowtail=normal, arrowhead=none", style="solid", length="1")
                        if options.mapSymmetric:
                            matches = [ leftMatch for leftMatch in leftMatches if leftMatch.otherSide in rightMatches ]
                            if len(matches) == 1:
                                addMatchEdge("red", "", matches[0])
                                haveMatched = True
                            elif not haveMatched and options.showMultiMaps:
                                for matchingSide in matches:
                                    addMatchEdge("orange", "", matchingSide)
                        else:
                            if len(leftMatches) == 1:
                                if len(rightMatches) == 1 and leftMatches[0].otherSide == rightMatches[0]: 
                                        addMatchEdge("red", "B", leftMatches[0])
                                        haveMatched = True
                                elif len(rightMatches) == 0:
                                    addMatchEdge("blue", "L", leftMatches[0])
                                    haveMatched = True
                            elif len(leftMatches) == 0 and len(rightMatches) == 1:
                                addMatchEdge("blue", "L", rightMatches[0])
                                haveMatched = True
                        
                            if not haveMatched and options.showMultiMaps:
                                for matchingSide in set(leftMatches):
                                    addMatchEdge("orange", "L", matchingSide)
                                for matchingSide in set(rightMatches):
                                    addMatchEdge("orange", "R", matchingSide)
        
    finishGraphFile(graphVizFileHandle)
    graphVizFileHandle.close()

if __name__ == '__main__':
    exit(main())  
    