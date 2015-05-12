#!/usr/bin/env bash

# readCoverage.sh: Make the coverage plots for the paper. Execute from the
# cluster run output directory. Takes an optional argument, which is the format
# to generate graphs in.

# Die on errors
set -e

# We want to put output in this directory
OUTDIR="paper"
# And save it in SVG unless someone tells us different.
GRAPH_FORMAT=${1-svg}
# We will make this data file (per-base coverage)
TSV="${OUTDIR}/readCoverage.tsv"
# And this plot image
GRAPH="${OUTDIR}/readCoverage.${GRAPH_FORMAT}"

# Make sure out directory exists.
mkdir -p ${OUTDIR}

# Start fresh
truncate -s 0 ${TSV}

for FILE in `ls *.counts`
do
    # Sort things out by category

    # Grab the scheme
    SCHEME=${FILE%.*}
    
    # We're going to count up how many mapped and unmapped bases it has.
    TOTAL_MAPPED=0
    TOTAL_UNMAPPED=0
    
    while read LINE || [[ -n $LINE ]]
    do
        
        # Split on spaces
        PARTS=(${LINE})
        
        # Parse out the category and count for each line
        CATEGORY="${PARTS[0]}"
        COUNT="${PARTS[1]}"
        
        if [[ ${CATEGORY} == *2unmapped ]]
        then
            # This one was unmapped. Add it in.
            TOTAL_UNMAPPED=$((${TOTAL_UNMAPPED} + ${COUNT}))
        elif [[ ${CATEGORY} == "!*" ]]
        then
            # Skip this one
            echo "Metadata category ${CATEGORY}"
        else
            # This one was mapped. Add it in.
            TOTAL_MAPPED=$((${TOTAL_MAPPED} + ${COUNT}))
        fi
        
        if [[ ${CATEGORY} == "gene2unmapped" ]]
        then
            # This is the last category from each genome. Save stats per genome.
            # TODO: this is a hack.
            
            TOTAL_BASES=$((${TOTAL_MAPPED} + ${TOTAL_UNMAPPED}))
            
            echo "${SCHEME} mapped ${TOTAL_MAPPED}/${TOTAL_BASES}"
            
            # Calculate coverage
            COVERAGE=$(echo "${TOTAL_MAPPED} / ( ${TOTAL_MAPPED} + ${TOTAL_UNMAPPED} )" | bc -l)
    
            # Save coverage point
            printf "${SCHEME}\t${COVERAGE}\n" >> "${TSV}"
            
            # Restart for next scheme
            TOTAL_MAPPED=0
            TOTAL_UNMAPPED=0
        fi
        
    done < ${FILE}
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
LABEL_OPTS+=("--category_labels" "BWA" "--category_labels" "$(printf 'BWA\nMax MAPQ')")

# And the weakly stable ones
GROUPING_OPTS+=("--grouping" "Weakly Stable")
CATEGORY_OPTS+=("--categories" \
    "ICnaturalHam3Mis2U" \
    "ICnaturalHam5Mis4U")
LABEL_OPTS+=("--category_labels" \
    "$(printf 'Weak\n$\\alpha^\\prime=3, \\beta^\\prime=2$')"
    "$(printf 'Weak\n$\\alpha^\\prime=5, \\beta^\\prime=4$')")
    
GROUPING_OPTS+=("--grouping" "Stable")
CATEGORY_OPTS+=("--categories" \
    "ICnaturalHam3Mis2" \
    "ICnaturalHam5Mis4")
LABEL_OPTS+=("--category_labels" \
    '$\alpha^\prime=3, \beta^\prime=2$'
    '$\alpha^\prime=5, \beta^\prime=4$')
    
# Do LASTZ as a group
GROUPING_OPTS+=("--grouping" "LASTZ")
CATEGORY_OPTS+=("--categories" "Lastz")
LABEL_OPTS+=("--category_labels" "$(printf 'LASTZ\nscore$\\geq$17,500')")

    
boxplot.py "${TSV}" \
    --x_label "Scheme Parameters" \
    --y_label "$(printf 'Portion of Bases\nMapped to Reference')" \
    --title "$(printf 'Read Coverage vs.\nMapping Scheme')" \
    "${GROUPING_OPTS[@]}" "${CATEGORY_OPTS[@]}" "${LABEL_OPTS[@]}" \
    --no_legend \
    --grouping_colors 'k' 'y' 'b' 'r' 'g' \
    --x_sideways \
    --no_n \
    --max 1.01 \
    --width 4 --height 4 \
    --save "${GRAPH}"
    
# We put the arrays in quotes and use @ above because that uses the array
# elements as tokens. If we omit *either* of those, Bash will re-split
# everything on spaces.
