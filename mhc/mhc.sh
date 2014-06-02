#!/usr/bin/env bash
set -e

# mhc.sh: Download some MHC haplotypes and index them.

for GI_NUMBER in 568335879 568335954 568335976 568335986 568335989 568335992 568335994
do
    if [ ! -e GI${GI_NUMBER}.fa ]
    then
        # Download each of the hg38 MHC alts
        wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${GI_NUMBER}&rettype=fasta" \
            -o GI${GI_NUMBER}.fa
    fi
done



if [ ! -e refmhc.fa ]
then
    # Get the FASTA for hg38 chr6:28,510,120-33,480,577 AKA GI568815592
    wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=568815592&strand=1&seq_start=28510120&seq_stop=33480577&rettype=fasta&retmode=text" \
        -o refmhc.fa
    
fi

echo "Indexing..."

time ../createIndex.sh --quiet --context 20 mhc-index *.fa
    
