#!/usr/bin/env bash
set -e
shopt -s extglob

# makeOrderPlots.sh: make plots in a directory of order comparison stats output.

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

boxplot.py precision.tsv --x_label "Merging Scheme" --y_label "Precision" --title "Precision vs. Merging Scheme" --x_sideways --save precisionVsScheme.png

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

boxplot.py recall.tsv --x_label "Merging Scheme" --y_label "Recall" --title "Recall vs. Merging Scheme" --x_sideways --save recallVsScheme.png

OUTFILE="orderAgreement.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls agreement.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${OUTFILE}
    done
    
done

boxplot.py orderAgreement.tsv --x_label "Merging Scheme" --y_label "Order Agreement (F score)" --title "Order Agreement vs. Merging Scheme" --x_sideways --save agreementVsScheme.png

OUTFILE="runtimes.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls runtime.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${OUTFILE}
    done
    
done

boxplot.py runtimes.tsv --x_label "Merging Scheme" --y_label "Runtime (seconds)" --title "Runtime vs. Merging Scheme" --x_sideways --save runtimeVsScheme.png

OUTFILE="refCoverages.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls refCoverage.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | cut -f 2 | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${OUTFILE}
    done
    
done
    
boxplot.py refCoverages.tsv --x_label "Merging Scheme" --y_label "Portion Aligned to refmhc" --title "refmhc Coverage vs. Merging Scheme" --x_sideways --save coverageVsScheme.png

mkdir -p spectrums

for SPECTRUM in spectrum.*
do
    # Grab the extension. See <http://unix.stackexchange.com/a/1574/44348>
    SCHEME=${SPECTRUM##*.}
    
    # Plot the histogram indication the tandem dupes
    histogram.py --min 0 --max 40 --bins 40 --logCounts ${SPECTRUM} --save spectrums/${SCHEME}.spectrum.png --x_label "Adjacency Component Size" --y_label "Occurrences" --title "${SCHEME} Component Spectrum" --label
done

OUTFILE="precisionRecall.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls truth.*`
do
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    # Type is the scheme name without the numbers (and any ps in numbers)
    SCHEME_TYPE=$(echo ${SCHEME} | sed 's/[0-9]p\?//g')
    
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

scatter.py ${OUTFILE} --x_label "Precision" --y_label "Recall" --title "Recall vs. Precision by Scheme" --save recallVsPrecision.png --max_x 1 --max_y 1

mkdir -p orderCoverage

for ORDERCOVERAGE in orderCoverage.*
do
    # Grab the extension. See <http://unix.stackexchange.com/a/1574/44348>
    SCHEME=${ORDERCOVERAGE##*.}
    
    # Plot a scatterplot
    scatter.py --save orderCoverage/${SCHEME}.orderCoverage.png --x_label "Addition Rank" --y_label "Coverage" --title "${SCHEME} Coverage vs. Order" --min_y 0.75 --max_y 1.0 "${ORDERCOVERAGE}"
done




