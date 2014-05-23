#!/usr/bin/env bash
# traversalScaling.sh: Produce TSV files of times for suffix tree traversal at
# different context lengths. Takes the FASTA to index as input.

# Die on errors
set -e

# Grab the FASTA to index.
FASTA="${1}"

# Make a directory to aggregate output from background subshells.
OUT_DIR=`mktemp -d`

for CONTEXT in {1..30}
do
    (
        # Grab a temp directory
        INDEX=`mktemp -d`
        
        # Pick a temp file to send our output to.
        OUT_FILE="${OUT_DIR}/${CONTEXT}.tsv"
        
        # Make a scratch file
        SCRATCH_FILE=`mktemp`
    
        # Run on the FASTA file for each context.
        createIndex/createIndex ${INDEX} ${FASTA} --context ${CONTEXT} &> ${SCRATCH_FILE}
        
        # Grab the line for the second traversal and extract the first number: wall clock time in seconds.
        TIME_TAKEN=`cat ${SCRATCH_FILE} | grep "timer - Level Index Construction" | sed 's/[^0-9]*\([0-9]\+\.[0-9]\+\).*/\1/'`
        
        # Print them out as a two-column TSV
        printf "${CONTEXT}\t${TIME_TAKEN}\n" >> "$OUT_FILE"
        
        echo "# Context ${CONTEXT} done in ${TIME_TAKEN}"
        
        # Get rid of the temp directory.
        rm -Rf ${INDEX}
        
        # And the scratch file
        rm ${SCRATCH_FILE}
    ) &
    
done

# Wait for all the subshells to finish.
wait

# Aggregate all the output
cat ${OUT_DIR}/*.tsv

# Remove aggregated results.
rm -Rf ${OUT_DIR}


