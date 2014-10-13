#!/usr/bin/env bash
set -e
shopt -s extglob

# makeGenePlots.sh: make plots of gene mapping info from compareGenes.py. Run in
# the output directory of the run.

rm -Rf counts/
mkdir counts/
for FILE in `ls *.counts`
do
    # Sort things out by category

    # Grab the scheme
    SCHEME=${FILE%.*}
    
    while read LINE || [[ -n $LINE ]]
    do
        
        # Split on spaces
        PARTS=(${LINE})
        
        # Parse out the category and count for each line
        CATEGORY="${PARTS[0]}"
        COUNT="${PARTS[1]}"
        
        # Save the count and scheme name to the file for the category
        printf "${SCHEME}\t${COUNT}\n" >> counts/${CATEGORY}.tsv
        
    done < ${FILE}
    
done

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

    barchart.py ${FILE} --x_label "Scheme" --x_sideways --y_label "Genes" --title "Genes with ${CATEGORY} mappings" --save genes/${CATEGORY}.genes.png
    
done

