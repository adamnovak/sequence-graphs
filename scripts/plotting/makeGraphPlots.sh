#!/usr/bin/env bash
set -e
shopt -s extglob

# makeGraphPlots.sh: make plots in a directory of multi-region scheme comparison
# stats output.

mkdir -p plots

truncate -s 0 plots/precision.tsv
truncate -s 0 plots/recall.tsv
truncate -s 0 plots/precisionRecall.tsv
truncate -s 0 plots/coverages.tsv
rm -Rf plots/checkgenes
mkdir -p plots/checkgenes

for REGION in `ls`
do
    # For each per-region directory
    if [ "${REGION}" == "plots" ]
    then
        # Skip the place where we put our plots
        continue
    fi
    
    for FILE in `ls ${REGION}/truth.*`
    do
        # We're going to get the precision and recall values.
        SCHEME=${FILE##*.}
        
        cat ${FILE} | cut -f 1 | while read LINE
        do
            # Grab the precisions
            printf "${SCHEME}\t${LINE}\n" >> plots/precision.tsv
        done
        
        cat ${FILE} | cut -f 2 | while read LINE
        do
            # Grab the recalls
            printf "${SCHEME}\t${LINE}\n" >> plots/recall.tsv
        done
        
        # Type is the scheme name without the numbers (with ps in them)
        SCHEME_TYPE=$(echo ${SCHEME} | sed 's/[0-9]\(p\)\?//g')
        
        if [ "${SCHEME_TYPE}" == "EzipMinEdRange" ]
        then
            # Label Camel scheme more nicely.
            SCHEME_TYPE="Camel"
        fi
        
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
        printf "${SCHEME_TYPE} in ${REGION}\t${MEAN_PRECISION}\t${MEAN_RECALL}\n" >> plots/precisionRecall.tsv
        
    done
    
    for FILE in `ls ${REGION}/coverage.*`
    do
        # Also the coverage values
        
        SCHEME=${FILE##*.}
        
        cat ${FILE} | cut -f 2 | while read LINE
        do
            # Tag each entry with its scheme
            printf "${SCHEME}\t${LINE}\n" >> plots/coverages.tsv
        
        done
    done
    
    for FILE in `ls ${REGION}/checkgenes.*`
    do
        # Also the gene mapping category counts

        # Grab the scheme
        SCHEME=${FILE##*.}
        
        # We want total membership in all categories (except background) for
        # this file, so we can make the other things be fractions. TODO: This
        # won't work after we get total_columns, or if you try to do multiple
        # alignments per scheme & region.
        TOTAL_MEMBERSHIPS=0
        
        while read LINE || [[ -n $LINE ]]
        do
            # Split on spaces
            PARTS=(${LINE})
            
            if [ "${PARTS[0]}" == "background" ]
            then
                # Skip columns with no genes at all, to try to normalize out
                # gene density.
                continue
            fi
            
            # Add in everything. Float math needs the bc tool.
            TOTAL_MEMBERSHIPS=$(echo ${TOTAL_MEMBERSHIPS} + ${PARTS[1]} | bc -l)
            
        done < ${FILE}
        
        while read LINE || [[ -n $LINE ]]
        do
            
            # Split on spaces
            PARTS=(${LINE})
            
            # Parse out the category and count for each line
            CATEGORY="${PARTS[0]}"
            COUNT="${PARTS[1]}"
            
            if [ "${CATEGORY}" == "background" ]
            then
                # Skip columns with no genes at all.
                continue
            fi
            
            # Divide by total
            PORTION=$(echo ${COUNT} / ${TOTAL_MEMBERSHIPS} | bc -l)
            
            # Save the portion and scheme name to the file for the category
            printf "${SCHEME}\t${PORTION}\n" >> plots/checkgenes/${CATEGORY}.tsv
            
        done < ${FILE}
        
    done

done

if [ -e plots/precision.tsv ]
then
    boxplot.py plots/precision.tsv --x_label "Merging Scheme" --y_label "Precision" --title "Precision vs. Merging Scheme" --x_sideways --save plots/precisionVsScheme.png
fi

if [ -e plots/recall.tsv ]
then
    boxplot.py plots/recall.tsv --x_label "Merging Scheme" --y_label "Recall" --title "Recall vs. Merging Scheme" --x_sideways --save plots/recallVsScheme.png
fi

if [ -e plots/coverages.tsv ]
then
    boxplot.py plots/coverages.tsv --x_label "Merging Scheme" --y_label "Portion Aligned to Reference" --title "Coverage vs. Merging Scheme" --x_sideways --save plots/coverageVsScheme.png
fi

if [ -e plots/precisionRecall.tsv ]
then
    scatter.py plots/precisionRecall.tsv --x_label "Precision" --y_label "Recall" --title "Recall vs. Precision by Scheme" --save plots/recallVsPrecision.png --max_x 1 --max_y 1 --tsv
fi

for FILE in `ls plots/checkgenes/*.tsv`
do

    # Now make the actual charts

    # Pull out the filename without directory
    FILENAME="${FILE##*/}"
    # And the category name from that
    CATEGORY="${FILENAME%.*}"

    boxplot.py ${FILE} --x_label "Scheme" --x_sideways --y_label "Count" --title "Occurrences of ${CATEGORY} mappings" --save plots/checkgenes/${CATEGORY}.checkgenes.png
    
done



