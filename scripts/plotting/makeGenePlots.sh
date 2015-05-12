#!/usr/bin/env bash
set -e
shopt -s extglob

# makeGenePlots.sh: make plots of gene mapping info from compareReads.py. Run in
# the output directory of the run.

rm -Rf counts/
mkdir counts/
rm -f coverage.tsv
rm -f wrongness-fraction.tsv
for FILE in `ls *.counts`
do
    # Sort things out by category

    # Grab the scheme
    SCHEME=${FILE%.*}
    
    # We're going to count up how many mapped and unmapped bases it has.
    TOTAL_MAPPED=0
    TOTAL_UNMAPPED=0
    
    # And how many bases are in genes
    TOTAL_GENE2=0
    
    while read LINE || [[ -n $LINE ]]
    do
        
        # Split on spaces
        PARTS=(${LINE})
        
        # Parse out the category and count for each line
        CATEGORY="${PARTS[0]}"
        COUNT="${PARTS[1]}"
        
        # Save the count and scheme name to the file for the category
        printf "${SCHEME}\t${COUNT}\n" >> counts/${CATEGORY}.tsv
        
        if [[ ${CATEGORY} == gene2* ]]
        then
            # This one was supposed to be in a gene. Add it so we can see what
            # fraction of gene bases are 2wrong later.
            TOTAL_GENE2=$((${TOTAL_GENE2} + ${COUNT}))
        fi
        
        if [[ ${CATEGORY} == *2unmapped ]]
        then
            # This one was unmapped. Add it in.
            TOTAL_UNMAPPED=$((${TOTAL_UNMAPPED} + ${COUNT}))
        elif [[ ${CATEGORY} == "!*" ]]
        then
            # Skip this one
            echo "Metadata category ${CATEGORY}"
        else
            # This one was mapped. Add it in.
            TOTAL_MAPPED=$((${TOTAL_MAPPED} + ${COUNT}))
        fi
        
        if [[ ${CATEGORY} == "gene2wrong" ]]
        then
            # Save the count for when we get to gene2unmapped, the last category
            # for the genome, so we can compute a 
            # portion-of-gene-bases-wrongly-mapped value.
            GENE2WRONG="${COUNT}"
        fi
        
        if [[ ${CATEGORY} == "gene2unmapped" ]]
        then
            # This is the last category from each genome. Save stats per genome.
            # TODO: this is a hack.
            
            # Calculate coverage
            COVERAGE=$(echo "${TOTAL_MAPPED} / ( ${TOTAL_MAPPED} + ${TOTAL_UNMAPPED} )" | bc -l)
    
            # Save coverage point
            printf "${SCHEME}\t${COVERAGE}\n" >> coverage.tsv
            
            # And do the same for portion of gene bases wrongly mapped
            WRONGNESS=$(echo "${GENE2WRONG} / ${TOTAL_GENE2}" | bc -l)
            printf "${SCHEME}\t${WRONGNESS}\n" >> wrongness-fraction.tsv
            
            # Restart for next scheme
            TOTAL_MAPPED=0
            TOTAL_UNMAPPED=0
            TOTAL_GENE2=0
            
        fi
        
    done < ${FILE}
    
    
    
done

# While we're at it we totaled up per scheme coverage, so plot that.
boxplot.py coverage.tsv --x_label "Scheme" --x_sideways --y_label "Coverage" --max 1.01 --title "Read Mapping Coverage" --save coverage.png

# Also the gene2wrong/gene2* wrongness ratios
boxplot.py wrongness-fraction.tsv --x_label "Scheme" --x_sideways --y_label "Fraction of Gene Bases Mapped to Wrong Gene" --min 0 --title "Read Mapping Gene Wrongness" --save wrongness.png

for FILE in `ls counts/*.tsv`
do

    # Now make the actual charts

    # Pull out the filename without directory
    FILENAME="${FILE##*/}"
    # And the category name from that
    CATEGORY="${FILENAME%.*}"

    boxplot.py ${FILE} --x_label "Scheme" --x_sideways --y_label "Count" --title "Occurrences of ${CATEGORY} mappings" --save counts/${CATEGORY}.counts.png
    
done

rm -Rf genes/
mkdir genes/
for FILE in `ls *.genes`
do
    # Sort things out by category

    # Grab the scheme
    SCHEME=${FILE%.*}
    
    while read LINE || [[ -n $LINE ]]
    do
        
        # Split on spaces
        PARTS=(${LINE})
        
        # Parse out the category and gene name
        CATEGORY="${PARTS[0]}"
        GENE="${PARTS[1]}"
        
        # Save one entry in that category
        # TODO: this is a huge hack
        printf "${SCHEME}\t1\n" >> genes/${CATEGORY}.tsv
        
    done < ${FILE}
    
done

for FILE in `ls genes/*.tsv`
do

    # Now make the actual charts

    # Pull out the filename without directory
    FILENAME="${FILE##*/}"
    # And the category name from that
    CATEGORY="${FILENAME%.*}"

    barchart.py ${FILE} --x_label "Scheme" --x_sideways --y_label "Gene Instances" --title "Genes with ${CATEGORY} mappings" --save genes/${CATEGORY}.genes.png
    
done

