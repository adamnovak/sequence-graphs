#!/usr/bin/env bash

set -e

# compareCamelAndCactus.sh: look at the HAL files from Camel and Cactus and see
# how they compare.

# What regions do we look at?
REGIONS="SMA MHC LRC_KIR"

function camelFilename()
{
    # Return the filename of the Camel HAL for the given upper-case region
    local REGION=$1
    
    echo "altRegions/${REGION}/graph.hal"
}

function cactusFilename()
{
    # Return the filename of the Cactus HAL for the given upper-case region
    local REGION=$1
    
    echo "cactus/${REGION,,}_work/${REGION,,}_star.hal"
}


function mkdtemp()
{
    # See <http://unix.stackexchange.com/a/84980>. Should work on OS X.
    mktemp -d 2>/dev/null || mktemp -d -t 'temp'
}


# We will keep series-labeled PR and coverage data here
MAIN_PR="precisionRecall.tsv"
MAIN_COVERAGE="coverage.tsv"

# Clear them out
rm -f ${MAIN_PR}
rm -f ${MAIN_COVERAGE}

# Make directories for our single-line files
PR_DIR=`mkdtemp`
COVERAGE_DIR=`mkdtemp`

for REGION in ${REGIONS}
do
    # We will use these files to hold the direct output of evaluateAlignment.
    TEMP_PR="${PR_DIR}/${REGION}.Camel"
    TEMP_COVERAGE="${COVERAGE_DIR}/${REGION}.Camel"
    
    # Get results for Camel
    mhc/evaluateAlignment.py \
        `camelFilename ${REGION}` \
        --truth altRegions/${REGION}/GRCAlignment.maf \
        --beds altRegions/${REGION}/genes/*/*.bed \
        --coverage_file ${TEMP_COVERAGE} \
        --precision_recall_file ${TEMP_PR} &
        
    # Set up to start Cactus instead
    
    TEMP_PR="${PR_DIR}/${REGION}.Cactus"
    TEMP_COVERAGE="${COVERAGE_DIR}/${REGION}.Cactus"
    
    # Get results for Cactus
    mhc/evaluateAlignment.py \
        `cactusFilename ${REGION}` \
        --truth altRegions/${REGION}/GRCAlignment.maf \
        --beds altRegions/${REGION}/genes/*/*.bed \
        --coverage_file ${TEMP_COVERAGE} \
        --precision_recall_file ${TEMP_PR} &
        
done

# Wait for all the evaluations to run
wait

# Tag the coverages by generator and save them
cat ${COVERAGE_DIR}/*.Camel | sed 's/^/Camel\t/' >> ${MAIN_COVERAGE}
cat ${COVERAGE_DIR}/*.Cactus | sed 's/^/Cactus\t/' >> ${MAIN_COVERAGE}

# Do the same for the precision/recall data
cat ${PR_DIR}/*.Camel | sed 's/^/Camel\t/' >> ${MAIN_PR}
cat ${PR_DIR}/*.Cactus | sed 's/^/Cactus\t/' >> ${MAIN_PR}

# Clean up carefuylly, without the possibility of an rm -Rf ""
rm ${COVERAGE_DIR}/*.Camel ${COVERAGE_DIR}/*.Cactus 
rmdir ${COVERAGE_DIR}
rm ${PR_DIR}/*.Camel ${PR_DIR}/*.Cactus
rmdir ${PR_DIR}







