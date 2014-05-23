#!/usr/bin/env bash
# traversalScaling.sh: Produce TSV files of times for suffix tree traversal at
# different context lengths. Takes the FASTA to index as input.

# Die on errors
set -e

# Grab the FASTA to index.
FASTA="${1}"

# Grab a temp directory
INDEX=`mktemp -d`

for CONTEXT in {1..30}
do
    # Run on the FASTA file for each context.
    echo "#Running for ${CONTEXT} bases"
    # Grab the line for the second traversal and extract the first number: wall clock time in seconds.
    TIME_TAKEN=`./createIndex.sh ${INDEX} ${FASTA} --context ${CONTEXT} 2>&1 | grep "[timer - Level Index Construction]" | sed 's/[^0-9]*\([0-9]\+\.[0-9]\+\).*/\1/'`
    
    # Print them out as a two-column TSV
    printf "${CONTEXT}\t${TIME_TAKEN}\n"

done

# Get rid of the temp directory.
rm -Rf ${INDEX}


