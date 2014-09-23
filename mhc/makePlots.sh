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
    
    # Fix up fractional values
    SCHEME_PARAMETER=$(printf "${SCHEME_PARAMETER}" | sed s/p/./)
    
    echo "Scheme ${SCHEME_TYPE} parameter ${SCHEME_PARAMETER}"
    
    # We want totals
    TOTAL_PRECISION=0
    TOTAL_RECALL=0
    LINE_COUNT=0
    
    while read LINE || [[ -n $LINE ]]
    do
        # Split on spaces
        PARTS=(${LINE})
        
        # Add in everything. Float math needs the bc tool.
        TOTAL_PRECISION=$(echo ${TOTAL_PRECISION} + ${PARTS[0]} | bc -l)
        TOTAL_RECALL=$(echo ${TOTAL_RECALL} + ${PARTS[1]} | bc -l)
        LINE_COUNT=$((${LINE_COUNT} + 1))
        
        #echo "At line ${LINE_COUNT} totals ${TOTAL_PRECISION}, ${TOTAL_RECALL}"
        
    done < ${FILE}
    
    # Calculate averages. We seem to get integer division unless we pass -l?
    MEAN_PRECISION=$(echo ${TOTAL_PRECISION} / ${LINE_COUNT} | bc -l)
    MEAN_RECALL=$(echo ${TOTAL_RECALL} / ${LINE_COUNT} | bc -l)
    
    #echo "Means ${MEAN_PRECISION}, ${MEAN_RECALL}"
    
    # Write a mean precision/recall point for this type of scheme. We'll get one
    # per parameter value.
    printf "${SCHEME_TYPE}\t${MEAN_PRECISION}\t${MEAN_RECALL}\n" >> ${OUTFILE}
    
done

scatter.py ${OUTFILE} --x_label "Precision" --y_label "Recall" --title "Recall vs. Precision by Scheme" --sparse_ticks --save recallVsPrecision.png
