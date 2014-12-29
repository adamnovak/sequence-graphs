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
    
boxplot.py coverages.tsv --x_label "Merging Scheme" --y_label "Portion Aligned to Reference" --title "Coverage vs. Merging Scheme" --x_sideways --save coverageVsScheme.png

mkdir -p spectrums

for SPECTRUM in spectrum.*
do
    # Grab the extension. See <http://unix.stackexchange.com/a/1574/44348>
    SCHEME=${SPECTRUM##*.}
    
    # Work out the total tandem count
    TANDEM_COUNT=$(sum.sh < tandem.${SCHEME}) 
    
    # Plot the histogram indication the tandem dupes
    histogram.py --min 0 --max 40 --bins 40 --logCounts ${SPECTRUM} --save spectrums/${SCHEME}.spectrum.png --x_label "Adjacency Component Size" --y_label "Occurrences" --title "${SCHEME} Component Spectrum" --label --redWeight 0 --redWeight 0 --redWeight 0 --redWeight 0 --redWeight ${TANDEM_COUNT}
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
    
    # Type is the scheme name without the numbers (with ps in them)
    SCHEME_TYPE=$(echo ${SCHEME} | sed 's/[0-9]\(p\)\?//g')
    
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

# Go get the size-4 adjacency components from the spectrums, and the tandem
# dupe counts, and see what portion of size-4 adjacency components are tandem
# dupes per scheme.
OUTFILE="tandemPortion.tsv"
truncate -s 0 ${OUTFILE}
for FILE in `ls tandem.*`
do

    # Grab the scheme
    SCHEME=${FILE##*.}

    # Make the spectrum filename
    SPECTRUM_FILE="spectrum.${SCHEME}"
    
    # Sum up the tandem duplications
    TOTAL_TANDEMS=0
    while read LINE || [[ -n $LINE ]]
    do
        TOTAL_TANDEMS=$((${TOTAL_TANDEMS} + ${LINE}))
    done < ${FILE}
    
    
    # Sum up the 4-part rearrangements
    TOTAL_4PART=0
    while read LINE || [[ -n $LINE ]]
    do
        
        # Split on spaces
        PARTS=(${LINE})
        
        if [[ ${PARTS[0]} -eq 4 ]]
        then
            TOTAL_4PART=$((${TOTAL_4PART} + ${PARTS[1]}))
        fi
    done < ${SPECTRUM_FILE}
    
    # Work out the fraction of 4-part rearrangements (i.e. 2-break operations)
    # that are tandem duplications by our count. Will miss tandems with any
    # differences, or any other things aligned.
    FRACTION=$(echo ${TOTAL_TANDEMS} / ${TOTAL_4PART} | bc -l)
    
    printf "${SCHEME}\t${FRACTION}\n" >> ${OUTFILE}
    
    echo "Scheme ${SCHEME} is ${TOTAL_TANDEMS} / ${TOTAL_4PART} = ${FRACTION} tandems"

done

barchart.py ${OUTFILE} --x_label "Merging Scheme" --y_label "2-break tandem portion" --max 1 --title "Portion Tandem Dupes vs. Merging Scheme" --x_sideways --save tandemPortion.png

# Make a boxplot of alignment mappings by scheme for each category, to display
# the result of classifying mappings according to what they say about genes.
# Still not sure exactly what the right way to display this is.
rm -Rf checkgenes/
mkdir checkgenes/
for FILE in `ls checkgenes.*`
do
    # Sort things out by category

    # Grab the scheme
    SCHEME=${FILE##*.}
    
    while read LINE || [[ -n $LINE ]]
    do
        
        # Split on spaces
        PARTS=(${LINE})
        
        # Parse out the category and count for each line
        CATEGORY="${PARTS[0]}"
        COUNT="${PARTS[1]}"
        
        # Save the count and scheme name to the file for the category
        printf "${SCHEME}\t${COUNT}\n" >> checkgenes/${CATEGORY}.tsv
        
    done < ${FILE}
    
done

for FILE in `ls checkgenes/*.tsv`
do

    # Now make the actual charts

    # Pull out the filename without directory
    FILENAME="${FILE##*/}"
    # And the category name from that
    CATEGORY="${FILENAME%.*}"

    boxplot.py ${FILE} --x_label "Scheme" --x_sideways --y_label "Count" --title "Occurrences of ${CATEGORY} mappings" --save checkgenes/${CATEGORY}.checkgenes.png
    
done


