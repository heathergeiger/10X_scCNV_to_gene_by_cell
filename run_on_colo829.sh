#!/bin/bash

curl -C - -O http://cf.10xgenomics.com/samples/cell-dna/1.0.0/colo829_G1_1k/colo829_G1_1k_node_cnv_calls.bed
mkdir wide_matrices
Rscript long_to_wide.R colo829_G1_1k_node_cnv_calls.bed wide_matrices

mkdir gene_matrices
for chr in {1..22}
do
./wide_to_genes.sh GRCh37_CellRanger_RNA_genes.bed wide_matrices $chr gene_matrices
echo Chromosome $chr gene x cell complete
done

for chr in {1..22}
do
cat gene_matrices/copy_number_genes_chrom_${chr}.csv
done | awk '!seen[$0]++' | sort -t ',' -k1,1 > gene_by_cell_matrix.csv

gzip gene_by_cell_matrix.csv
