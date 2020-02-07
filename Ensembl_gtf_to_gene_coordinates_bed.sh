#!/bin/bash

# sed -e '1,5d' removes first 5 lines, which are header.
# The tr just removes quotes and semicolons to clean things up a bit.
# awk '$3 == "gene"' is bc third field is entry type. Choose only gene entries.
# $14 selects the 14th field as separated by whitespace. Fields 9 and onward contain the gene ID, version, name, etc. For gene entries, the gene symbol should always be this field.
# The seen[$14]++ counts the number of times each gene symbol is seen, so that genes with different IDs but the same symbol can be given a unique identifier.
# This unique identifier is equivalent to R make.unique function with sep="_".

GTF=$1
BED=$2

sed -e '1,5d' $GTF | tr -d '\"' | tr -d ';' | awk '$3 == "gene"' | \
 awk '{ OFS="\t"}{seen[$14]++;if(seen[$14] == 1)print $1,$4,$5,$14;else print $1,$4,$5,$14"_"seen[$14] - 1}' > $BED
