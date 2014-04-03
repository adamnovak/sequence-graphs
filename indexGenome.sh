#!/usr/bin/env bash
# indexGenome.sh: download and time the indexing of the human genome.
set -e

DIR=indexGenome

mkdir -p ${DIR}
cd ${DIR}

# What's the index to accumulate in?
BASENAME=hg19

rm -f ${BASENAME}.*

for CHROM in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
    # Download each chromosome
    wget http://hgdownload.cse.ucsc.edu/goldenpath/${BASENAME}/chromosomes/${CHROM}.fa.gz
    gunzip ${CHROM}.fa.gz
    FASTA=${CHROM}.fa
    
    HAPLOTYPE=${CHROM}.haplotype
    
    # Do forward strand
    echo "${CHROM} forward strand..."
    cat $FASTA | grep -v ">" | tr acgtn ACGTN | tr -d '\n' > ${HAPLOTYPE}
    printf '\0' >> ${HAPLOTYPE}
    
    # Do reverse complement
    echo "${CHROM} reverse strand..."
    cat $FASTA | grep -v ">" | tr acgtn ACGTN | tr ACGT TGCA | tr -d '\n' | rev | tr -d '\n' >> ${HAPLOTYPE}
    printf '\0' >> ${HAPLOTYPE}
    
done

ls *.haplotype > list.txt
parallel_build list.txt ${BASENAME} 10



        

