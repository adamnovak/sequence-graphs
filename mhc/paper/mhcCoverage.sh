#!/usr/bin/env bash

# mhcCoverage.sh: Make the coverage (vs. precision) plot for the MHC for the
# paper. Execute from the cluster run output directory.

# Die on errors
set -e

# We want to put output in this directory
OUTDIR="paper"
# We will make this data file
TSV="${OUTDIR}/mhcCoverage.tsv"
# And this plot image
GRAPH="${OUTDIR}/mhcCoverage.png"

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
        # Make a point for this instantiation of a scheme. I need the file with
        # coverage values for the scheme, and the file with precision and recall
        # values.
        SCHEME="$1"; shift
        POINT_COVERAGE_FILE="coverage.${SCHEME}"
        POINT_PRECISION_FILE="truth.${SCHEME}"
        
        # We want totals for this file to make an average point
        TOTAL_COVERAGE=0
        LINE_COUNT=0
        
        while read LINE || [[ -n $LINE ]]
        do
            # Split on spaces
            PARTS=(${LINE})
            
            # Add in recall. Float math needs the bc tool.
            TOTAL_COVERAGE=$(echo ${TOTAL_COVERAGE} + ${PARTS[1]} | bc -l)
            LINE_COUNT=$((${LINE_COUNT} + 1))
            
        done < ${POINT_COVERAGE_FILE}
        
        # Calculate average coverage. We seem to get integer division unless we
        # pass -l?
        MEAN_COVERAGE=$(echo ${TOTAL_COVERAGE} / ${LINE_COUNT} | bc -l)
        
        # Reset to count up precision
        TOTAL_PRECISION=0
        LINE_COUNT=0
        
        while read LINE || [[ -n $LINE ]]
        do
            # Split on spaces
            PARTS=(${LINE})
            
            # Add in precision. Float math needs the bc tool.
            TOTAL_PRECISION=$(echo ${TOTAL_PRECISION} + ${PARTS[0]} | bc -l)
            LINE_COUNT=$((${LINE_COUNT} + 1))
            
        done < ${POINT_PRECISION_FILE}
        
        # Calculate average recall. We seem to get integer division unless we
        # pass -l?
        MEAN_PRECISION=$(echo ${TOTAL_PRECISION} / ${LINE_COUNT} | bc -l)
        
        # Write a mean coverage vsprecision point for this series.
        printf "${SERIES_NAME}\t${MEAN_PRECISION}\t${MEAN_COVERAGE}\n" >> ${FILE_NAME}
        
    done
}

# Make sure the out directory exists
mkdir -p "${OUTDIR}"

# And make sure the data file is empty
truncate -s 0 ${TSV}

# Define all the series we want

# Vary the minimum length
series "${TSV}" "Flat Length Threshold" \
    ICnaturalMin20 \
    ICnaturalMin50 \
    ICnaturalMin100 \
    ICnaturalMin150 \
    ICnaturalMin200 \
    ICnaturalMin250

# Do all the other series
for HAMMING_CLEARANCE in {1..6}
do
    # For each minimum Hamming clearance we used...
    
    # Make an array of points
    POINTS=("ICnaturalHam${HAMMING_CLEARANCE}")
    
    for HAMMING_DISTANCE in $(seq 1 $((HAMMING_CLEARANCE - 1)))
    do
        # For each maximum Hamming distance we used that's strictly smaller than
        # the clearance...
        
        # Make a point at this maximum Hamming distance
        POINTS+=("${POINTS[0]}Mis${HAMMING_DISTANCE}")
    done
    
    # Make a series for this Hamming clearance
    series "${TSV}" "Hamming Clearance ${HAMMING_CLEARANCE}" ${POINTS[*]}
done

# Make the actual plot
scatter.py ${TSV} --tsv --no_sort \
    --colors 'b' 'g' 'r' 'c' 'm' 'y' 'k' \
    --markers 'o' 's' '>' 'v' '+' '_' 'D'  \
    --x_label "Precision" --y_label "Coverage" \
    --title "Coverage vs. Precision by Filter Criterion" \
    --max_x 1 --max_y 1 \
    --lines \
    --legend_overlay best \
    --save ${GRAPH}
