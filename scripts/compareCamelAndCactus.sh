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
MAIN_GENES="geneCategories.tsv"

# Clear them out
rm -f ${MAIN_PR}
rm -f ${MAIN_COVERAGE}
rm -f ${MAIN_GENES}

# Make directories for our individual files
PR_DIR=`mkdtemp`
COVERAGE_DIR=`mkdtemp`
GENES_DIR=`mkdtemp`

for REGION in ${REGIONS}
do
    # We will use these files to hold the direct output of evaluateAlignment.
    TEMP_PR="${PR_DIR}/${REGION}.Camel"
    TEMP_COVERAGE="${COVERAGE_DIR}/${REGION}.Camel"
    TEMP_GENES="${GENES_DIR}/${REGION}.Camel"
    
    # Get results for Camel
    scripts/evaluateAlignment.py \
        `camelFilename ${REGION}` \
        --truth altRegions/${REGION}/GRCAlignment.maf \
        --beds altRegions/${REGION}/genes/*/*.bed \
        --coverage_file ${TEMP_COVERAGE} \
        --precision_recall_file ${TEMP_PR} \
        --gene_category_file ${TEMP_GENES} \
        --tag Camel ${REGION} &
        
    # Set up to start Cactus instead
    
    TEMP_PR="${PR_DIR}/${REGION}.Cactus"
    TEMP_COVERAGE="${COVERAGE_DIR}/${REGION}.Cactus"
    TEMP_GENES="${GENES_DIR}/${REGION}.Cactus"
    
    # Get results for Cactus
    scripts/evaluateAlignment.py \
        `cactusFilename ${REGION}` \
        --truth altRegions/${REGION}/GRCAlignment.maf \
        --beds altRegions/${REGION}/genes/*/*.bed \
        --coverage_file ${TEMP_COVERAGE} \
        --precision_recall_file ${TEMP_PR} \
        --gene_category_file ${TEMP_GENES} \
        --tag Cactus ${REGION} &
        
done

# Wait for all the evaluations to run
wait

for REGION in ${REGIONS}
do

    # Tag the coverages by generator and region and save them
    cat ${COVERAGE_DIR}/${REGION}.Camel >> ${MAIN_COVERAGE}
    cat ${COVERAGE_DIR}/${REGION}.Cactus >> ${MAIN_COVERAGE}

    # Do the same for the precision/recall data
    cat ${PR_DIR}/${REGION}.Camel >> ${MAIN_PR}
    cat ${PR_DIR}/${REGION}.Cactus >> ${MAIN_PR}
    
    # And the genes
    cat ${GENES_DIR}/${REGION}.Camel >> ${MAIN_GENES}
    cat ${GENES_DIR}/${REGION}.Cactus >> ${MAIN_GENES}
    
done
    

# Clean up carefuylly, without the possibility of an rm -Rf ""
rm ${COVERAGE_DIR}/*.Camel ${COVERAGE_DIR}/*.Cactus 
rmdir ${COVERAGE_DIR}
rm ${PR_DIR}/*.Camel ${PR_DIR}/*.Cactus
rmdir ${PR_DIR}
rm ${GENES_DIR}/*.Camel ${GENES_DIR}/*.Cactus
rmdir ${GENES_DIR}






