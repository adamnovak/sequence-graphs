#!/usr/bin/env bash
# indexGenome.sh: download and time the indexing of the human genome.
set -e

mkdir indexGenome
cd indexGenome

for CHROM in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
    # Download each chromosome
    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/${CHROM}.fa.gz
    gunzip ${CHROM}.fa.gz
    FASTA=${CHROM}.fa
    
    
    #FASTA="/pod/podstore/data/reference/${CHROM}.fa"
    
    # Do forward strand
    echo "${CHROM} forward strand..."
    cat $FASTA | grep -v ">" | tr acgtn ACGTN | tr -d '\n' >> hg19.haplotype
    printf '\0' >> hg19.haplotype
    
    # Do reverse complement
    echo "${CHROM} reverse strand..."
    cat $FASTA | grep -v ">" | rev | tr acgtn ACGTN | tr ACGT TGCA | tr -d '\n' >> hg19.haplotype
    printf '\0' >> hg19.haplotype
done

time build_rlcsa hg19.haplotype 10


        

