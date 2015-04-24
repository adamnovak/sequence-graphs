#!/usr/bin/env bash
set -e

# runCactus.sh: Run Cactus on all the regions.

for REGION in SMA MHC LRC_KIR; do
    
    # Lower-case it
    REGION_LOWER="${REGION,,}"
    
    # Grab the FASTAs
    FASTA_FILES=`ls /hive/users/anovak/sgdev/altRegions/${REGION}/*.fa | grep -v graph`
    
    # Decide on where to put the star tree file
    STAR_FILENAME="${REGION_LOWER}_star.txt"
    
    # Decide on what work directory to use
    WORK_DIR="${REGION_LOWER}_work"
    
    # And where to save the HAL
    HAL="${WORK_DIR}/${REGION_LOWER}_star.hal"
    
    # Start the tree line
    printf "(" >${STAR_FILENAME}
    
    # We are the first thing in the parens
    IS_FIRST=1
    
    for FASTA_FILE in ${FASTA_FILES}; do
        # For each FASTA
        
        if [ "$IS_FIRST" == "1" ]; then
            # Next time we aren't first
            IS_FIRST=0
        else
            # Do the seperator
            printf ", " >>${STAR_FILENAME}
        fi
    
        # We need the base name. See
        # <http://stackoverflow.com/a/2664746/402891>
        
        FILE_NAME="${FASTA_FILE##*/}"
        BASE_NAME="${FILE_NAME%.*}"
        
        printf "${BASE_NAME}:0.001" >>${STAR_FILENAME}
        
    done
    
    printf ");\n\n" >>${STAR_FILENAME}
    
    for FASTA_FILE in ${FASTA_FILES}; do
        FILE_NAME="${FASTA_FILE##*/}"
        BASE_NAME="${FILE_NAME%.*}"
        
        printf "${BASE_NAME} ${FASTA_FILE}\n" >>${STAR_FILENAME}
    done

    /cluster/home/hickey/genomes/progressiveCactus/bin/runProgressiveCactus.sh --batchSystem parasol --bigBatchSystem singleMachine --defaultMemory 8589934593 --bigMemoryThreshold 8589934592 --bigMaxMemory 893353197568 --bigMaxCpus 25 --maxThreads 25 --parasolCommand='/cluster/home/jcarmstr/bin/parasol -host=ku' --retryCount 3 ${STAR_FILENAME} ${WORK_DIR} ${HAL} --logInfo

done
