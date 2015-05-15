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
    
    if [ ! -e "${OUT_DIR}/${REGION}.tsv" ]
    then
        # Only re-run the actual evaluation if the TSV is missing.
        
        # Evaluate its reference's mapability to itself and save the result
        createIndex/evaluateMapability "${OUT_DIR}/index" \
            "altRegions/${REGION}/ref.fa" "${OUT_DIR}/${REGION}.tsv"
    fi
        
    # Make a histogram such that we can compare total area and get a good idea
    # of where things are.
    histogram.py "${OUT_DIR}/${REGION}.tsv" \
        --title "${REGION} Minimum Context Length Distribution" \
        --x_label "Minimum Context Length" --y_label "Count" \
        --log --min 1 --max 10000 --save "${OUT_DIR}/${REGION}.png"

done
