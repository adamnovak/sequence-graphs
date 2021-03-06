#!/usr/bin/env bash

# mhcCoverage.sh: Make the coverage vs. precision and coverage barchart for the
# MHC for the paper. Execute from the cluster run output directory. Takes an
# optional argument, which is the format to generate graphs in.

# Die on errors
set -e

# We want to put output in this directory
OUTDIR="paper"
# And save it in SVG unless someone tells us different.
GRAPH_FORMAT=${1-svg}
# We will make this data file
TSV="${OUTDIR}/mhcCoverage.tsv"
# And this plot image
GRAPH="${OUTDIR}/mhcCoverage.${GRAPH_FORMAT}"

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
    series "${TSV}" "\$\\\\alpha^\\\\prime = ${HAMMING_CLEARANCE}\$" ${POINTS[*]}
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
    --no_n \
    --save ${GRAPH}
    
# OK now we do the bar chart

# We will make this data file
TSV="${OUTDIR}/mhcCoverageBar.tsv"
# And this plot image
GRAPH="${OUTDIR}/mhcCoverageBar.${GRAPH_FORMAT}"

# Start fresh
truncate -s 0 ${TSV}

for FILE in `ls coverage.*`
do
    # Collect all the coverage info into one file.
    
    # Grab the scheme
    SCHEME=${FILE##*.}
    
    cat ${FILE} | cut -f 2 | while read LINE
    do
        # Tag each entry with its scheme
        printf "${SCHEME}\t${LINE}\n" >> ${TSV}
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

# Do just the min length one
GROUPING_OPTS+=("--grouping" "Minimum Length")
CATEGORY_OPTS+=("--categories" \
    "ICnaturalMin250" \
    "ICnaturalMin200" \
    "ICnaturalMin150" \
    "ICnaturalMin100" \
    "ICnaturalMin50" \
    "ICnaturalMin20")
LABEL_OPTS+=("--category_labels" '$l \geq 250$' '$l \geq 200$' '$l \geq 150$' '$l \geq 100$' '$l \geq 50$' '$l \geq 20$')

# Do all the other series
for HAMMING_CLEARANCE in {1..6}
do
    # For each minimum Hamming clearance we used...
    
    # Title the grouping
    GROUPING_OPTS+=("--grouping" "\$\\alpha^\\prime=${HAMMING_CLEARANCE}\$")
    
    # Add category for 0 mismatches
    CATEGORY_OPTS+=("--categories" "ICnaturalHam${HAMMING_CLEARANCE}")
    
    # And label
    LABEL_OPTS+=("--category_labels" "\$\\alpha^\\prime=${HAMMING_CLEARANCE}\$ \$\\beta^\\prime=0\$")
    
    for HAMMING_DISTANCE in $(seq 1 $((HAMMING_CLEARANCE - 1)))
    do
        # For each maximum Hamming distance we used that's strictly smaller than
        # the clearance...
        
        # Add and label a category for this number of mismatches
        CATEGORY_OPTS+=("ICnaturalHam${HAMMING_CLEARANCE}Mis${HAMMING_DISTANCE}")
        LABEL_OPTS+=("\$\\alpha^\\prime=${HAMMING_CLEARANCE}\$ \$\\beta^\\prime=${HAMMING_DISTANCE}\$")
    done
done
    
# Work out the average coverage for the GRC alignment We can't just directly use
# the MAF coverage output, since we don't want to count Ns in the query against
# the coverage, since those can naver be aligned.

# Count up the total coverage and number of genomes, to divide.
TOTAL_COVERAGE=0
TOTAL_GENOMES=0

while read LINE || [[ -n $LINE ]]
do
    # Grab the clipped and cut table of MAF coverage info from the
    # mafPairCoverage tool, and parse each line.

    # Split on spaces
    PARTS=(${LINE})
    
    # Pull out contig, character count, and aligned bases
    CONTIG=${PARTS[0]}
    LENGTH=${PARTS[1]}
    ALIGNED=${PARTS[2]}
    
    # Find the number of Ns
    N_COUNT=$(cat ../${CONTIG}.fa  | grep -v ">" | grep -o "N" | tr -d '[:space:]' | wc -c)
    
    # Find the coverage
    COVERAGE=$(echo "${ALIGNED} / (${LENGTH} - ${N_COUNT})" | bc -l)
    
    echo "Genome ${CONTIG} has N-corrected coverage ${COVERAGE} = ${ALIGNED} / (${LENGTH} - ${N_COUNT})"
    
    # Add in for averaging
    TOTAL_COVERAGE=$(echo "${TOTAL_COVERAGE} + ${COVERAGE}" | bc -l)
    TOTAL_GENOMES=$((${TOTAL_GENOMES} + 1))
    
done <<< "`mafPairCoverage --maf ../GRCAltAlignmentNoMismatch.maf --seq1 refmhc --seq2 \"GI*\" | tail -n +9 | sed 's/\s\s\+/ /g' | sed 's/^\s//g' | cut -f 1,3,4 -d ' '`"

GRC_AVERAGE=$(echo ${TOTAL_COVERAGE} / ${TOTAL_GENOMES} | bc -l)
    
echo "GRC average coverage: ${GRC_AVERAGE}"
    
boxplot.py "${TSV}" \
    --x_label "Scheme Parameter (Min Length or Clearance and Tolerance)" \
    --y_label "Portion Aligned to Reference" \
    --title "Coverage vs. Mapping Scheme" \
    "${GROUPING_OPTS[@]}" "${CATEGORY_OPTS[@]}" "${LABEL_OPTS[@]}" \
    --grouping_colors 'b' 'g' 'r' 'c' 'm' 'y' 'k' \
    --no_legend \
    --x_sideways \
    --no_n \
    --hline ${GRC_AVERAGE} \
    --save "${GRAPH}"
    
# We put the arrays in quotes and use @ above because that uses the array
# elements as tokens. If we omit *either* of those, Bash will re-split
# everything on spaces.
