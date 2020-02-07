#!/bin/bash

#Accept gene coordinates BED file, wide copy number matrices directory path, chromosome name, and output file path as command line inputs.

gene_coords=$1
copy_matrix_dir=$2
chr=$3
output_dir=$4

#Create a temp directory.

if [ ! -e copy_num_${chr}_tmp ]
then
mkdir copy_num_${chr}_tmp
fi

tmp_dir=copy_num_${chr}_tmp

#Get BED file containing annotated genes for the selected chromosome only.

grep -w $chr $gene_coords > $tmp_dir/annotated_genes.bed

#Get BED file for all intervals in DNA where there are at least 500 cells non-NA.
#This is somewhat of an arbitary threshold. But most intervals are going to have either all cells non-NA, all cells NA, or only a handful NA.
#So the exact cutoff doesn't really matter so much.
#AWK code to require 500+ cells NA adapted from here: https://unix.stackexchange.com/questions/552731/sum-occurrences-of-string-for-each-row-of-a-matrix-using-awk

sed -e '1,1d' $copy_matrix_dir/copy_number_chrom_${chr}.csv | awk -F "," '{OFS="\t"}{for(i=4;i<=NF;i++)if($i != "NA")c++;if(c >= 500)print $1,$2,$3;c=0}' | tr -d '\"' > $tmp_dir/mappable_DNA_intervals_coordinates.bed
sed -e '1,1d' $copy_matrix_dir/copy_number_chrom_${chr}.csv | cut -d ',' -f 4- | awk -F "," '{OFS="\t"}{for(i=4;i<=NF;i++)if($i != "NA")c++;if(c >= 500)print $0;c=0}' > $tmp_dir/mappable_DNA_intervals_values.csv

paste $tmp_dir/mappable_DNA_intervals_coordinates.bed $tmp_dir/mappable_DNA_intervals_values.csv > $tmp_dir/mappable_DNA_intervals_coordinates_and_values.bed

#Attach number of bases covered by DNA with a semicolon to the name field in annotated_genes.bed.
#This will also remove any genes in annotated_genes.bed that are actually located entirely in unmappable regions.

#Required to have bedtools installed and in the path here.

bedtools intersect -wo -a $tmp_dir/annotated_genes.bed -b $tmp_dir/mappable_DNA_intervals_coordinates_and_values.bed | \
 awk '{full_bed=$1":"$2":"$3":"$4;covered[full_bed]+=$NF} END {for(gene in covered)print gene"\t"covered[gene]}' | \
 awk '{ OFS="\t"}{split($1,a,":");print a[1],a[2],a[3],a[4]";"$NF}' > $tmp_dir/annotated_genes_incl_num_bases_covered_by_DNA_in_name_field.bed

#Next use intersectBed to get an intermediate format including the gene name per interval and the weight for the interval.
#In the intersectBed output, field 4 will contain the gene name and # covered bases in DNA, separated by a semicolon.
#Field 8 contains all the actual values in the interval in the DNA.
#Field 9 contains the number of bases overlap between the gene and the interval.

#The weight should be the number of bases overlap between the gene and the interval over the total covered bases for the gene in the DNA.

#Also add the appropriate cell names to the header of this file.

head -n 1 $copy_matrix_dir/copy_number_chrom_${chr}.csv | cut -d ',' -f 4- | awk -F "," '{ OFS=","}{printf("%s,%s,","\"""gene""\"","\"""interval_weight""\"");print $0}' > $tmp_dir/mappable_DNA_intervals_gene_names_and_weights_instead_of_coords.csv

bedtools intersect -wo -a $tmp_dir/annotated_genes_incl_num_bases_covered_by_DNA_in_name_field.bed -b $tmp_dir/mappable_DNA_intervals_coordinates_and_values.bed | \
 awk '{ OFS=","}{split($4,a,";");gene=a[1];total_bases=a[2];covered_bases=$9;values=$8;print "\""gene"\"",covered_bases/total_bases,values}' >> $tmp_dir/mappable_DNA_intervals_gene_names_and_weights_instead_of_coords.csv

#Create output directory. 

if [ ! -e $output_dir ]
then
mkdir $output_dir
fi

#Finally, use an R script to get weighted average per gene per cell across all intervals.
#This R script still needs to be written.

Rscript wide_to_genes.R $chr $output_dir

#Clean up by removing temp directory.

rm -r $tmp_dir
