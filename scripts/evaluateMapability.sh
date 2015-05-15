#!/usr/bin/env bash

# evaluateMapability.sh: Find and plot context lengths for the references for
# all regions.

if [ $# -ne 1 ]
then
    echo "Please specify output directory."
    exit 1
fi

# Grab the output directory
OUT_DIR="${1}"

# Make it exist
mkdir -p "${OUT_DIR}"

for REGION in SMA MHC LRC_KIR BRCA1 BRCA2
do
    # For each region
    
    # Evaluate its reference's mapability to itself and save the result
    createIndex/evaluateMapability "${OUT_DIR}/index" \
        "altRegions/${REGION}/ref.fa" "${OUT_DIR}/${REGION}.tsv"
        
    # Make a line-chart histogram
    histogram.py "${OUT_DIR}/${REGION}.tsv" \
        --title "${REGION} Minimum Context Length Distribution" \
        --x_label "Minimum Context Length" --y_label "Count" \
        --log --log_counts --bins 100 --line --save "${OUT_DIR}/${REGION}.png"

done
