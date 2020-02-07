Scripts to run on 10X single-cell CNV output by CellRanger DNA pipeline.

Result will be a gene x cell matrix.

Included is a BED file CellRanger_RNA_genes.bed with gene coordinates for the gene annotation from the CellRanger single-cell RNA pipeline, GRCh38.

Can run script Ensembl_gtf_to_gene_coordinates_bed.sh to get a similar file for different annotation.

Next, create a directory to store chromosome-level wide matrices (20kb intervals x n cells) in CSV format.

	mkdir wide_matrices

Run an R script to get wide matrices per chromosome output to this directory given a BED file and the directory name.
Here we will run on the colo829 BED file downloaded from 10X website.

	Rscript long_to_wide.R colo829_G1_1k_node_cnv_calls.bed wide_matrices

Run a BASH script to match up gene coordinates (in BED file) to 20kb intervals, and get copy number per gene using a weighted average across intervals.
This script accepts gene coordinates BED file, wide copy number matrices directory path, chromosome name, and output file path as command line inputs.
So need to run on each chromosome in a for loop.
Also need BEDtools installed and in the PATH before this.

	mkdir gene_matrices
	for chr in {1..22};do ./wide_to_genes.sh CellRanger_RNA_genes.bed wide_matrices $chr gene_matrices;echo Chromosome $chr gene x cell complete;done

Gene x cell matrices will now be available per chromosome for chr 1-22 under gene_matrices directory.
Clean up by removing wide matrices, which we do not need anymore.

	rm -r wide_matrices

Concatenate across all chromosome-level gene x cell matrices, then can remove those files as well.

	for chr in {1..22};do cat gene_matrices/copy_number_genes_chrom_${chr}.csv;done | awk '!seen[$0]++' | sort -t ',' -k1,1 > gene_by_cell_matrix.csv #awk command gets unique values without sorting.
	gzip gene_by_cell_matrix.csv
	rm -r gene_matrices
