#!/usr/bin/env bash
set -e
shopt -s extglob

# makePlots.sh: make plots in a directory of scheme comparison stats output.

OUTFILE="precision.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls truth.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | cut -f 1 | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${OUTFILE}
    done
    
done

boxplot.py ${OUTFILE} --x_label "Merging Scheme" --y_label "Precision" --title "Precision vs. Merging Scheme" --x_sideways --save precisionVsScheme.png

OUTFILE="recall.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls truth.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | cut -f 2 | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${OUTFILE}
    done
    
done

boxplot.py ${OUTFILE} --x_label "Merging Scheme" --y_label "Recall" --title "Recall vs. Merging Scheme" --x_sideways --save recallVsScheme.png

OUTFILE="coverages.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls coverage.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | cut -f 2 | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${OUTFILE}
    done
    
done
    
boxplot.py ${OUTFILE} --x_label "Merging Scheme" --y_label "Portion Aligned to Reference" --title "Coverage vs. Merging Scheme" --x_sideways --save coverageVsScheme.png

mkdir -p spectrums

for SPECTRUM in spectrum.*
do
    # Grab the extension. See <http://unix.stackexchange.com/a/1574/44348>
    SCHEME=${SPECTRUM##*.}
    histogram.py --min 0 --max 40 --bins 40 --logCounts ${SPECTRUM} --save spectrums/${SCHEME}.spectrum.png --x_label "Adjacency Component Size" --y_label "Occurrences" --title "${SCHEME} Component Spectrum" --label
done

OUTFILE="indelLengths.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls indels.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" 
    done >> ${OUTFILE}
    
done

barchart.py ${OUTFILE} --x_label "Merging Scheme" --y_label "Mean Indel Length" --title "Mean Length vs. Merging Scheme" --x_sideways --save meanIndels.png

OUTFILE="tandemDupes.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls tandem.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${OUTFILE}
    done
    
done
    
boxplot.py ${OUTFILE} --x_label "Merging Scheme" --y_label "Tandem Duplication Count" --title "Tandem Duplications vs. Merging Scheme" --x_sideways --save tandemsVsScheme.png

OUTFILE="precisionRecall.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls truth.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    # Break it into scheme proper and parameter setting
    # Parameter is whatever's after all the alpha stuff (2p0)
    SCHEME_PARAMETER=${SCHEME##*([A-Za-z])}
    # Type is whatever's before that (ICMult)
    SCHEME_TYPE=${SCHEME%${SCHEME_PARAMETER}}
    
    echo "Scheme ${SCHEME_TYPE} parameter ${SCHEME_PARAMETER}"
    
    cat ${FILE} | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME_TYPE}\t${LINE}\n" >> ${OUTFILE}
    done
    
done

scatter.py ${OUTFILE} --x_label "Precision" --y_label "Recall" --title "Recall vs. Precision by Scheme" --sparse_ticks --save recallVsPrecision.png
