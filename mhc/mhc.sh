#!/usr/bin/env bash
set -e

# mhc.sh: Download some MHC haplotypes and index them.

# Download the MHCs
./getMHCs.sh

echo "Indexing..."

# Make sure to put refmhc.fa first so it becomes the reference against which
# everything maps.

time ../createIndex/createIndex --context 100 --scheme greedy --alignment mhc.c2h --alignmentFasta mhc.c2h.fasta --degrees degrees.txt mhc-index refmhc.fa GI*.fa

echo "Making tree..."

# Make a star tree off of rootSeq
FASTA_SEQ_NAMES=$(cat *.fa | grep ">" | cut -f 1 -d " " | sed 's/>//' | 
    tr '\n' ',')
# Drop the training comma. See <http://www.cyberciti.biz/faq/bash-remove-last-
# character-from-string-line-word/>
FASTA_SEQ_NAMES=${FASTA_SEQ_NAMES%?}

# Make the actual tree
TREE="(${FASTA_SEQ_NAMES})rootSeq;"

echo ${TREE}

echo "Creating HAL..."

rm -f mhc.hal
halAppendCactusSubtree mhc.c2h mhc.c2h.fasta "${TREE}" mhc.hal

echo "Creating AssemblyHub..."

rm -Rf hub
rm -Rf tree
hal2assemblyHub.py --jobTree tree mhc.hal hub --hub=mhc --shortLabel="All MHCs" --lod --bedDirs genes --tabBed
    
