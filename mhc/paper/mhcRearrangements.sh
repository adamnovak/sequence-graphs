#!/usr/bin/env bash

# mhcRearrangements.sh: Make a plot about rearrangements. Execute from the
# cluster run output directory. Takes an optional argument, which is the format
# to generate graphs in.

# Die on errors
set -e

# We want to put output in this directory
OUTDIR="paper"
# And save it in SVG unless someone tells us different.
GRAPH_FORMAT=${1-svg}
# We will make this data file
TSV="${OUTDIR}/mhcRearrangements.tsv"
# And this plot image
GRAPH="${OUTDIR}/mhcRearrangements.${GRAPH_FORMAT}"

# Make sure the out directory exists
mkdir -p "${OUTDIR}"

# And make sure the data file is empty
truncate -s 0 ${TSV}

# What scheme are we interested in?
SCHEME="ICnaturalHam6Mis5"

# Work out the total tandem count
TANDEM_COUNT=$(sum.sh < tandem.${SCHEME}) 

# Plot the histogram indication the tandem dupes
histogram.py "spectrum.${SCHEME}" \
    --logCounts \
    --x_min 0 --x_max 25 --bins 25 \
    --y_max 1000000 \
    --x_label "Adjacency Component Size" --y_label "Occurrences" \
    --title "$(printf 'Adjacency Component Size Spectrum\n6 Clearance, 5 Mismatches')" \
    --label \
    --redWeight 0 --redWeight 0 --redWeight 0 --redWeight 0 --redWeight "${TANDEM_COUNT}" \
    --no_n \
    --save "${GRAPH}"
