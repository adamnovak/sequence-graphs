#!/usr/bin/env bash
# runRegion.sh: Create a graph, HAL, and assembly hub for an alt region
# Usage: runRegion.sh <directory with the FASTAs>
# Uses ref.fa and GI*.fa from that directory
# Output will all be placed in that directory as well

# Stop on error
set -e

# Where should we look for input and save output?
DIR="${1}"

if [ ! -d "${DIR}" ]; then
    # We got a directory that isn't real
    echo "Directory ${DIR} is bad!"
    exit 1
fi

# Make the graph
./createIndex/createIndex "${DIR}/index" --mapType zip --context 25 --minEditBound 5 "${DIR}/ref.fa" "${DIR}"/GI*.fa --alignment "${DIR}/graph.c2h" --alignmentFasta "${DIR}/graph.fa" --mapStats "${DIR}/stats.txt" --lastGraph "${DIR}/graph.lastgraph"

# Make the tree
TREE="($(ls ${DIR}/GI*.fa | xargs -n 1 basename | tr '\n' ',' | sed 's/.fa//g')ref)rootSeq;"
echo "Making HAL with tree ${TREE}"

# Make the HAL
rm -f "${DIR}/graph.hal"
halAppendCactusSubtree --inMemory "${DIR}/graph.c2h" "${DIR}/graph.fa" "${TREE}" "${DIR}/graph.hal"

# Make the assembly hub
rm -Rf "${DIR}/tree" "${DIR}/hub"
hal2assemblyHub.py "${DIR}/graph.hal" "${DIR}/hub" --jobTree "${DIR}/tree" --maxThreads 32 --cpHalFileToOut --noUcscNames


