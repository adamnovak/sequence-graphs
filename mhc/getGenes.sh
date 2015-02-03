#!/usr/bin/env bash
set -e

# getGenes.sh: pull genes to use for gene analysis from hgsql. You need to have
# hgsql set up to hit against the UCSC browser database for this to be useful.

mkdir -p genes/refmhc

CONTIG=chr6
START=28510119
END=33480576
# Go get genes between start and end on refmhc. Don't forget to filter out
# things that have dummy CDS starts at tx start, since those are pseudogenes.
hgsql -e "SELECT \"refmhc\", hg38.knownGene.txStart - ${START}, hg38.knownGene.txEnd - ${START}, hg38.kgXref.geneSymbol, 0, hg38.knownGene.strand FROM hg38.knownGene LEFT OUTER JOIN hg38.kgXref ON hg38.knownGene.name = hg38.kgXref.kgID WHERE hg38.knownGene.txStart != hg38.knownGene.cdsStart AND hg38.knownGene.chrom = \"${CONTIG}\" AND hg38.knownGene.txStart >= ${START} AND hg38.knownGene.txEnd < ${END};" | sed 1d > genes/refmhc/refmhc.bed

# Just get all the genes on the other contigs
UCSC_NAMES=(chr6_GL000250v2_alt chr6_GL000251v2_alt chr6_GL000252v2_alt chr6_GL000253v2_alt chr6_GL000254v2_alt chr6_GL000255v2_alt chr6_GL000256v2_alt)
GI_NAMES=(GI568335879 GI568335954 GI568335976 GI568335986 GI568335989 GI568335992 GI568335994)

for CONTIG_NUMBER in {0..6}
do

# Make a directory for genes for this genome
mkdir -p genes/${GI_NAMES[CONTIG_NUMBER]}

# Then write the gene BED
hgsql -e "SELECT \"${GI_NAMES[CONTIG_NUMBER]}\", hg38.knownGene.txStart, hg38.knownGene.txEnd, hg38.kgXref.geneSymbol, 0, hg38.knownGene.strand FROM hg38.knownGene LEFT OUTER JOIN hg38.kgXref ON hg38.knownGene.name = hg38.kgXref.kgID WHERE hg38.knownGene.txStart != hg38.knownGene.cdsStart AND hg38.knownGene.chrom = \"${UCSC_NAMES[CONTIG_NUMBER]}\";" | sed 1d > genes/${GI_NAMES[CONTIG_NUMBER]}/${GI_NAMES[CONTIG_NUMBER]}.bed

done
