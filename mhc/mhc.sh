#!/usr/bin/env bash
set -e

# mhc.sh: Download some MHC haplotypes and index them.

for GI_NUMBER in 568335879 568335954 568335976 568335986 568335989 568335992 568335994
do
    if [ ! -e GI${GI_NUMBER}.fa ]
    then
        # Download each of the hg38 MHC alts
        wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${GI_NUMBER}&rettype=fasta" \
            -O GI${GI_NUMBER}.fa
            
        # Set their FASTA headers to something nice and clean with no special characters.
        # TODO: Support multiple sequences.
        sed -i "s/>.*/>GI${GI_NUMBER}/" GI${GI_NUMBER}.fa
            
    fi
done



if [ ! -e refmhc.fa ]
then
    # Get the FASTA for hg38 chr6 (AKA GI568815592) 28,510,120-33,480,577
    wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=568815592&strand=1&seq_start=28510120&seq_stop=33480577&rettype=fasta&retmode=text" \
        -O refmhc.fa
        
    # Set the FASTA header to something nice and clean with no special characters.
    sed -i "s/>.*/>refmhc/" refmhc.fa
    
fi

echo "Indexing..."

# Make sure to put refmhc.fa first so it becomes the reference against which
# everything maps.

valgrind --tool=callgrind ../createIndex/createIndex --context 100 --scheme greedy --alignment mhc.c2h --alignmentFasta mhc.c2h.fasta --degrees degrees.txt mhc-index refmhc.fa GI*.fa

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

#echo "Creating AssemblyHub..."

#rm -Rf hub
#rm -Rf tree
#hal2assemblyHub.py --jobTree tree mhc.hal hub --hub=mhc --shortLabel="All MHCs" --lod --bedDirs genes --tabBed
    
