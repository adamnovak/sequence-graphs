#!/usr/bin/env bash
# mhcPrecisionRecall.sh: Make the recall vs. precision plots for the MHC for the
# paper. Execute from the cluster run output directory.

# Die on errors
set -e

# We want to put output in this directory
OUTDIR="paper"
# We will make this data file
TSV="${OUTDIR}/mhcPrecisionRecall.tsv"
# And this plot image
GRAPH="${OUTDIR}/mhcPrecisionRecall.png"

function series {
    # Append a series to the given file with the given name, consisting of the
    # schemes specified by the rest of the arguments, connected in order by a
    # line.
    
    # What file name should we use?
    FILE_NAME="$1"; shift
    
    # What series name should we use?
    SERIES_NAME="$1"; shift
    
    echo "Series ${SERIES_NAME}: ${*}"
    
    while (( $# ))
    do
        # Make a point for this instantiation of a scheme.
        POINT_FILE="$1"; shift
        
        # We want totals for this file to make an average point
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
            
        done < ${POINT_FILE}
        
        # Calculate averages. We seem to get integer division unless we pass -l?
        MEAN_PRECISION=$(echo ${TOTAL_PRECISION} / ${LINE_COUNT} | bc -l)
        MEAN_RECALL=$(echo ${TOTAL_RECALL} / ${LINE_COUNT} | bc -l)
        
        # Write a mean precision/recall point for this series.
        printf "${SERIES_NAME}\t${MEAN_PRECISION}\t${MEAN_RECALL}\n" >> ${FILE_NAME}
        
    done
}

# Make sure the out directory exists
mkdir -p "${OUTDIR}"

# And make sure the data file is empty
truncate -s 0 ${TSV}

# Define all the series we want

# Vary the minimum length
series "${TSV}" "Min Length" \
    truth.ICnaturalMin20 \
    truth.ICnaturalMin50 \
    truth.ICnaturalMin100 \
    truth.ICnaturalMin150 \
    truth.ICnaturalMin200 \
    truth.ICnaturalMin250

# Do all the other series
for HAMMING_CLEARANCE in {1..6}
do
    # For each minimum Hamming clearance we used...
    
    # Make an array of points
    POINTS=("truth.ICnaturalHam${HAMMING_CLEARANCE}")
    
    for HAMMING_DISTANCE in $(seq 2 ${HAMMING_CLEARANCE})
    do
        # For each maximum Hamming distance we used...
        
        # Make a point at this maximum Hamming distance
        POINTS+=("${POINTS[0]}Mis${HAMMING_DISTANCE}")
    done
    
    # Make a series for this Hamming clearance
    series "${TSV}" "Hamming Clearance ${HAMMING_CLEARANCE}" ${POINTS[*]}
done

# Make the actual plot
scatter.py ${TSV} --tsv --x_label "Precision" --y_label "Recall" --title "Recall vs. Precision by Scheme" --save ${GRAPH} --max_x 1 --max_y 1 --lines
