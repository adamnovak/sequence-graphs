#!/usr/bin/env bash

# continuousDeployment.sh: Take the Cactus and Camel HAL files, run them through
# hal2sg, and put them where Maciek can put a server up.

# TODO: make this script re-do the servers, maybe by keying into the AWS commands?

set -e

# What regions do we look at?
REGIONS="SMA MHC LRC_KIR BRCA1 BRCA2"

function camelFilename()
{
    # Return the filename of the Camel HAL for the given upper-case region
    local REGION=$1
    
    # Here is where we pick the scheme to publish
    echo "output/out3/${REGION}/zipCmis2Min20Ed2Range100.hal"
}

function cactusFilename()
{
    # Return the filename of the Cactus HAL for the given upper-case region
    local REGION=$1
    
    # We pull these from Glenn
    echo "/hive/users/hickey/ga4gh/cactus/output/${REGION,,}_star.hal"
}

# Where should we save our hal2sg output?
OUT_DIR="output/deploy"

for REGION in ${REGIONS}
do
    # For each region

    # Where does this region's output go?
    REGION_OUT_DIR="${OUT_DIR}/${REGION}"
    
    # What HAL should Cactus have made?
    CACTUS_HAL=$(cactusFilename ${REGION})
    if [ -e ${CACTUS_HAL} ]
    then
        # It exists.
        echo "Found Cactus results for ${REGION}"
        
        # Make a cactus directory
        mkdir -p "${REGION_OUT_DIR}/cactus"
        
        # Fill it up with the generated database. TODO: We're supposed to have
        # --refGenome ref and --noAncestors, but apparently those aren't allowed
        # for some reason?
        hal2sg "${CACTUS_HAL}" "${REGION_OUT_DIR}/cactus/database.fa" \
            "${REGION_OUT_DIR}/cactus/database.sql"
            
    else
        echo "No Cactus results available for ${REGION}"
    fi
    
    # What HAL should Camel have made?
    CAMEL_HAL=$(camelFilename ${REGION})
    if [ -e ${CAMEL_HAL} ]
    then
        # It exists
        echo "Found Camel results for ${REGION}"
        
        # Make a cactus directory
        mkdir -p "${REGION_OUT_DIR}/camel"
        
        # Fill it up with the generated database
        hal2sg "${CAMEL_HAL}" "${REGION_OUT_DIR}/camel/database.fa" \
            "${REGION_OUT_DIR}/camel/database.sql"
    
    else
        echo "No Camel results available for ${REGION}"
    fi

done

echo "Databases ready for deployment in ${OUT_DIR}"

