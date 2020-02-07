#Load libraries.

library(dplyr)

#Accept chromosome name and output directory as command line inputs.

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]
output_dir <- args[2]

#Read in and format data.

tmp_dir <- paste0("copy_num_",chr,"_tmp")

info <- read.csv(paste0(tmp_dir,"/mappable_DNA_intervals_gene_names_and_weights_instead_of_coords.csv"),
	header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)

genes_and_weights <- info[,1:2]
values <- info[,3:ncol(info)]
values[values > 6] <- 6 #Truncate max copy number for an interval to 6.

#For each row, multiply all values by the appropriate weight for the gene.

values <- sweep(values,1,genes_and_weights$interval_weight,FUN="*")

#For each column, sum values for each gene now that we have multiplied by weights.

info <- data.frame(gene = genes_and_weights$gene,values,check.names=FALSE,stringsAsFactors=FALSE)
info <- info %>% group_by(gene) %>% summarise_all(sum)
info <- data.frame(info,check.names=FALSE,stringsAsFactors=FALSE)
rownames(info) <- info[,1]
info <- info[,2:ncol(info)]

#Round to nearest 4 digits.

info <- round(info,digits=4)

#Print.

write.csv(info,file=paste0(output_dir,"/copy_number_genes_chrom_",chr,".csv"),quote=TRUE)
