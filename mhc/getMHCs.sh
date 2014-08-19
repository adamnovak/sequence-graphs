#!/usr/bin/env bash
set -e

# getMHCs.sh: Download some MHC haplotypes.

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
