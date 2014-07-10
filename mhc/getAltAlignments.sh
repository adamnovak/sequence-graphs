#!/usr/bin/env bash
# getAltAlignments.sh: Download all the alt alignments and convert them to PSL 
# with the GI sequence names.

set -e

# We need to know the accession.version-style names of all the MHC alts
ACCESSION_VERSIONS[1]=GL000250.2
ACCESSION_VERSIONS[2]=GL000251.2
ACCESSION_VERSIONS[3]=GL000252.2
ACCESSION_VERSIONS[4]=GL000253.2
ACCESSION_VERSIONS[5]=GL000254.2
ACCESSION_VERSIONS[6]=GL000255.2
ACCESSION_VERSIONS[7]=GL000256.2

# Also the GI IDs that we use internally in the MHC alignment
GI_IDS[1]=GI568335879
GI_IDS[2]=GI568335954
GI_IDS[3]=GI568335976
GI_IDS[4]=GI568335986
GI_IDS[5]=GI568335989
GI_IDS[6]=GI568335992
GI_IDS[7]=GI568335994

# Also the UCSC chrom.sizes names
UCSC_NAMES[1]=chr6_GL000250v2_alt
UCSC_NAMES[2]=chr6_GL000251v2_alt
UCSC_NAMES[3]=chr6_GL000252v2_alt
UCSC_NAMES[4]=chr6_GL000253v2_alt
UCSC_NAMES[5]=chr6_GL000254v2_alt
UCSC_NAMES[6]=chr6_GL000255v2_alt
UCSC_NAMES[7]=chr6_GL000256v2_alt


for ALT_INDEX in {1..7}
do

    # Download the alignment in GFF3 format
    wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/ALT_REF_LOCI_${ALT_INDEX}/alt_scaffolds/alignments/${ACCESSION_VERSIONS[ALT_INDEX]}_CM000668.2.gff -O ${GI_IDS[ALT_INDEX]}.gff3
    
    # Rename the reference to what we call it (chr6)
    sed -i "s/CM000668\.2/chr6/g" ${GI_IDS[ALT_INDEX]}.gff3
    
    # Escape the period in the accession.version name. See
    # <http://stackoverflow.com/a/2705678/402891> for the sed pattern.
    ACCESSION_PATTERN=$(echo ${ACCESSION_VERSIONS[ALT_INDEX]} | sed -e 's/[]\/$*.^|[]/\\&/g')
    
    # Rename the alt to the UCSC name (chr6_GL000250v2_alt or whatever)
    sed -i "s/${ACCESSION_PATTERN}/${UCSC_NAMES[ALT_INDEX]}/g" ${GI_IDS[ALT_INDEX]}.gff3
    
    # Convert it to PSL
    gff3ToPsl hg38.chrom.sizes ${GI_IDS[ALT_INDEX]}.gff3 ${GI_IDS[ALT_INDEX]}.psl
    
    # Fix up the PSL to have our GIwhatever name instead of the UCSC name for
    # the alt.
    sed -i "s/${UCSC_NAMES[ALT_INDEX]}/${GI_IDS[ALT_INDEX]}/g" ${GI_IDS[ALT_INDEX]}.psl
    
done



