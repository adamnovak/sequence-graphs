#!/usr/bin/env bash

# readGenes.sh: Make the plots of bases and genes that are gene2wrong for the
# paper. Execute from the cluster run output directory. Takes an optional
# argument, which is the format to generate graphs in.

# Die on errors
set -e

# We want to put output in this directory
OUTDIR="paper"
# And save it in SVG unless someone tells us different.
GRAPH_FORMAT=${1-svg}
# We will make this data file (per-gene wrongness)
TSV="${OUTDIR}/readGenes.tsv"
# And this plot image
GRAPH="${OUTDIR}/readGenes.${GRAPH_FORMAT}"

# Make sure out directory exists.
mkdir -p ${OUTDIR}

# Start fresh
truncate -s 0 ${TSV}

for FILE in `ls *.genes`
do
    # Sort things out by category

    # Grab the scheme
    SCHEME=${FILE%.*}
    
    # What genomes exist?
    GENOMES=$(cat "${FILE}" | grep gene2 | cut -f5 | sort | uniq)
    
    for GENOME in ${GENOMES}
    do
        # For each genome
    
        # How many total genes are there?
        TOTAL_GENES=$(cat "${FILE}" | grep "${GENOME}" | grep gene2 | cut -f2 | sort | uniq | wc -l)
        
        # How many genes with gene2wrong mappings from them are there?
        WRONG_GENES=$(cat "${FILE}" | grep "${GENOME}" | grep gene2wrong | cut -f2 | sort | uniq | wc -l)
        
        # How many of the genes have gene2wrong mapings then?
        FRACTION=$(echo "${WRONG_GENES} / ${TOTAL_GENES}" | bc -l)
        
        echo "${WRONG_GENES} / ${TOTAL_GENES} wrong in ${GENOME} under ${SCHEME}"
        
        # Make a datapoint for this genome and scheme.
        printf "${SCHEME}\t${FRACTION}\n" >> "${TSV}"
        
    done
    
done

# We can specify "groupings" of "categories", where each category becomes a
# boxplot column and each grouping a collection of same. This is accomplished by
# three options, each of which can be specified multiple times, where the nth
# specification of each is associated with the nth grouping. We want to group up
# all the schemes, so we need to generate these options.

# This array stores the "--grouping" options
GROUPING_OPTS=()
# This array stores the "--categories" options
CATEGORY_OPTS=()
# This array stores the "--category_labels" options
LABEL_OPTS=()

# Do BWA as a group
GROUPING_OPTS+=("--grouping" "BWA" "--grouping" "BWAStrict")
CATEGORY_OPTS+=("--categories" "BWA" "--categories" "BWAStrict")
LABEL_OPTS+=("--category_labels" "BWA" "--category_labels" 'BWA $Q \geq 60$')

# And the weakly stable ones
GROUPING_OPTS+=("--grouping" "Weakly Stable")
CATEGORY_OPTS+=("--categories" \
    "ICnaturalHam3Mis2U" \
    "ICnaturalHam5Mis4U")
LABEL_OPTS+=("--category_labels" \
    '$\alpha=3, \beta=2$'
    '$\alpha=5, \beta=4$')
    
GROUPING_OPTS+=("--grouping" "Stable")
CATEGORY_OPTS+=("--categories" \
    "ICnaturalHam3Mis2" \
    "ICnaturalHam5Mis4")
LABEL_OPTS+=("--category_labels" \
    '$\alpha=3, \beta=2$'
    '$\alpha=5, \beta=4$')
    
boxplot.py "${TSV}" \
    --x_label "Scheme Parameters" \
    --y_label "Portion of Genes with Mappings to Paralogs" \
    --title "Paralog Mapping vs. Mapping Scheme" \
    "${GROUPING_OPTS[@]}" "${CATEGORY_OPTS[@]}" "${LABEL_OPTS[@]}" \
    --legend_overlay 'best' \
    --grouping_colors 'k' 'y' 'b' 'r' \
    --x_sideways \
    --no_n \
    --save "${GRAPH}"
    
# We put the arrays in quotes and use @ above because that uses the array
# elements as tokens. If we omit *either* of those, Bash will re-split
# everything on spaces.
