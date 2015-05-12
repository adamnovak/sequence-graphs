#!/usr/bin/env bash
# mergeGenome.sh: Merge a whole set of chromosomes into an idnex
set -e

DIR=indexGenome3
DIR2=indexGenome2

cd ${DIR}

# What's the index to accumulate in?
BASENAME=hg19

FILES=""

for CHROM in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
    FILES="${FILES} ../${DIR2}/${CHROM}.haplotype"
done

echo merge_rlcsa -10 ${BASENAME} ${FILES}

merge_rlcsa -10 ${BASENAME} ${FILES}


        

