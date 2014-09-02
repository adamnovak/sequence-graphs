#!/usr/bin/env bash
set -e

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
