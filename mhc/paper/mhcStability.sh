#!/usr/bin/env bash

# mhcStability.sh: Make the recall vs. precision plots showing how stability
# does nothing for MHCs. Execute from the cluster run output directory. Takes an
# optional argument, which is the format to generate graphs in.

# Die on errors
set -e

# We want to put output in this directory
OUTDIR="paper"
# And save it in SVG unless someone tells us different.
GRAPH_FORMAT=${1-svg}
# We will make this data file
TSV="${OUTDIR}/mhcStability.tsv"
# And this plot image
GRAPH="${OUTDIR}/mhcStability.${GRAPH_FORMAT}"

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
        POINT_FILE="truth.${1}"; shift
        
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

series "${TSV}" "Weakly Stable" \
    ICnaturalHam5Mis1U \
    ICnaturalHam5Mis2U \
    ICnaturalHam5Mis3U \
    ICnaturalHam5Mis4U
    
series "${TSV}" "Stable" \
    ICnaturalHam5Mis1 \
    ICnaturalHam5Mis2 \
    ICnaturalHam5Mis3 \
    ICnaturalHam5Mis4

# Make the actual plot
scatter.py ${TSV} --tsv --no_sort \
    --colors 'r' 'k'  \
    --markers 'd' 'D'  \
    --x_label "Precision" --y_label "Recall" \
    --title "$(printf 'Recall vs. Precision With Strong and Weak Stability\n$\\alpha = 5$')" \
    --max_x 1 --max_y 1 \
    --lines \
    --legend_overlay best \
    --no_n \
    --save ${GRAPH}
