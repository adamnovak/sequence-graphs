#!/bin/sh

###This script generates examples that are used in the theory paper (after postprocessing).

set -e

outDir=/Users/benedictpaten/Desktop

g0Label="Linear Sequence Graph"
g0='AAGCTACTGCC AGGCT!' 

g1Label="Sequence Graph"
g1='AAGCTACTGCC AATCTACTCC'

#Figure 1: Oriented sequence graph, with Ns and circular contig?
echo python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig1.dot "${g0}" "${g0Label}"  --showIDs 0 
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig1.dot "${g0}" "${g0Label}"  --showIDs 0 
dot -Tpdf ${outDir}/fig1.dot > ${outDir}/fig1.pdf

#Figure 2: Sequence graph (not linear)  - developed from first example.
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig2.dot "${g1}" "${g1Label}"  --showIDs 0 --mergeContigs 0 --usePhasedContexts
dot -Tpdf ${outDir}/fig2.dot > ${outDir}/fig2.pdf

#Figure 3: Linear genome context sets for simplicity shown for single scaffold.
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig3.dot "${g0}" "${g0Label} with Unique Context Strings"  --showContextSets 0 
dot -Tpdf ${outDir}/fig3.dot > ${outDir}/fig3.pdf

g2Label="Target Linear Sequence Graph"
g2='AAGCTACTGCC' 

g3Label="Input Linear Sequence Graph"
g3='AATCTACTCC'

g4Label="Input Linear Sequence Graph"
g4='AAGCTACTGCC' 

#Figure 4: Unique matching examples - isomorphic and variant mapping of linear input sequence grapgs to linear target sequence graphs
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig4.dot "${g2}" "${g2Label}" "${g3}" "${g3Label}" "${g4}" "${g4Label}" --showContextSets 0 --showMultiMaps --targetSequenceGraphs 0
dot -Tpdf ${outDir}/fig4.dot > ${outDir}/fig4.pdf

g5Label="Input Linear Sequence Graph"
g5='AACCTACTGCC'

g6Label="Target Sequence Graph"
g6='AAGCTACTGCC AATCTACTCC' 

#Figure 5: Mapping of a linear sequence graph to a non-linear graph. 
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig5.dot "${g6}" "${g6Label}" "${g3}" "${g3Label}" "${g4}" "${g4Label}" "${g5}" "${g5Label}" --showContextSets 0 --mergeContigs 0 --usePhasedContexts --targetSequenceGraphs 0
dot -Tpdf ${outDir}/fig5.dot > ${outDir}/fig5.pdf

g7Label="Input Sequence Graph"
g7='AAGCTACTGCC AATCTACTCC AACCTACTGCC'

#Figure 6:Mapping of a non-linear graph to a non-linear graph
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig6.dot "${g6}" "${g6Label}" "${g7}" "${g7Label}"  --showContextSets 0 --mergeContigs="0 1" --usePhasedContexts --targetSequenceGraphs 0
dot -Tpdf ${outDir}/fig6.dot > ${outDir}/fig6.pdf

g8Label="Sequence Graph"
g8='AAGCTACTGCC AATCTACTCC'

#Figure 7:Non symmetric merging.
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig7.dot "${g8}" "${g8Label} p=1" "${g8}" "${g8Label} p=3" "${g8}" "${g8Label} p=5" "${g8}" "${g8Label}" --mergeForM "0=1 1=3 2=5" --showContextSets "0 1 2 3" --noMapping --usePhasedContexts 
dot -Tpdf ${outDir}/fig7.dot > ${outDir}/fig7.pdf

g9Label="Input Sequence Graph"
g9='AAGCTACTCC'

#Figure 8: Hierarchy example
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig8.dot "${g8}" "${g8Label} p=3" "${g8}" "${g8Label}" "${g8}" "${g8Label}" "${g9}" "${g9Label}" --showMultiMaps --mergeContigs "0 1" --mergeForM "0=3" --showContextSets '0 1 2' --showIDs '0 1' --usePhasedContexts --showOnlyLowestMaps  
dot -Tpdf ${outDir}/fig8.dot > ${outDir}/fig8.pdf

#Figure 9:Nonsymmetric merging - hierarchy
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig9.dot "${g8}" "${g8Label} p=1" "${g8}" "${g8Label} p=3" "${g8}" "${g8Label} p=5" "${g8}" "${g8Label}" --mergeForM "0=1 1=3 2=5" --showContextSets "0 1 2 3" --usePhasedContexts --showOnlyLowestMaps  
dot -Tpdf ${outDir}/fig9.dot > ${outDir}/fig9.pdf

###Symmetric merging examples

g3Label="Input Linear Sequence Graph"
g3='AATCTACTCC!'

g4Label="Input Linear Sequence Graph"
g4='AAGCTACTGCC!' 

g5Label="Input Linear Sequence Graph"
g5='AACCTACTGCC!'

g8Label="Sequence Graph"
g8='AAGCTACTGCC! AATCTACTCC!'

#Figure 10:Symmetric merging - hierarchy
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig10.dot "${g8}" "${g8Label} p=1" "${g8}" "${g8Label} p=2" "${g8}" "${g8Label} p=3" "${g3}" "${g3Label}" "${g4}" "${g4Label}" "${g5}" "${g5Label}" --mergeForM "0=1 1=2 2=3" --mergeSymmetric --mapSymmetric --showContextSets "0 1 2" --usePhasedContexts --showOnlyLowestMaps --targetSequenceGraphs "0 1 2" 
dot -Tpdf ${outDir}/fig10.dot > ${outDir}/fig10.pdf

#Figure 11:Symmetric merging
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig11.dot "${g8}" "${g8Label} p=2" "${g3}" "${g3Label}" "${g4}" "${g4Label}" "${g5}" "${g5Label}" --mergeForM "0=2" --mergeSymmetric --mapSymmetric --showContextSets "0" --usePhasedContexts --showOnlyLowestMaps --targetSequenceGraphs 0
dot -Tpdf ${outDir}/fig11.dot > ${outDir}/fig11.pdf

g7Label="Input Sequence Graph"
g7='AAGCTACTGCC! AATCTACTCC!'

#Figure 12:Symmetric merging and mapping non-linear input graph
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig12.dot "${g8}" "${g8Label} p=2" "${g7}" "${g7Label}" --mergeForM "0=2" --mergeSymmetric --mapSymmetric --showContextSets "0" --usePhasedContexts --mergeContigs 1 --showOnlyLowestMaps --targetSequenceGraphs 0
dot -Tpdf ${outDir}/fig12.dot > ${outDir}/fig12.pdf

#Figure 13:Symmetric merging - hierarchy mapping a non-linear input
python src/scripts/toyReferenceGraphGenerator.py --graphVizFile ${outDir}/fig13.dot "${g8}" "${g8Label} p=1" "${g8}" "${g8Label} p=2" "${g8}" "${g8Label} p=3" "${g7}" "${g7Label}" --mergeForM "0=1 1=2 2=3" --mergeSymmetric --mapSymmetric --showContextSets "0 1 2" --mergeContigs 3 --usePhasedContexts --showOnlyLowestMaps --targetSequenceGraphs "0 1 2" 
dot -Tpdf ${outDir}/fig13.dot > ${outDir}/fig13.pdf


