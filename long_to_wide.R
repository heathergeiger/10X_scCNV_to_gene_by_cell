#Turn off scientific notation.
#Otherwise large coordinates get changed from e.g. 200,000,000 to 2e8.

options(scipen=999)

#Accept input file path and output directory as command line argument.
#Output directory should already exist.

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_dir <- args[2]

#Read in and format file, 

info <- read.table(input_file,header=FALSE,stringsAsFactors=FALSE)
colnames(info) <- c("chr","start","end","id","copy_num","confidence")

#Remove cell IDs that are not actually individual cells.

info <- info[info$id <= (max(info$id)/2),]

#Remove intervals with a confidence score less than 16.

info <- info[info$confidence >= 16,]

#Convert cell IDs to character now that that is complete.

info$id <- as.character(info$id)

#For each chromosome, convert from wide to long.

for(chr_for_subset in unique(info$chr))
{
	info_subset <- info[info$chr == chr_for_subset,]
	
	#Order by id, then by start and end coordinates (reverse order for end coordinates).
	#Don't think this should happen, but if a cell has both a larger region and sub-regions, we would want the sub-regions to overwrite the relevant parts of the larger region.

	info_subset <- info_subset[order(as.numeric(info_subset$id),info_subset$start,-1*info_subset$end),]

	#Get min start coordinate and max end coordinate.

	min_start <- min(info_subset$start)
	max_end <- max(info_subset$end)

	#Create a matrix of NAs for all 20kb intervals between min start coordinate and max end coordinate.

	cell_ids <- unique(info$id)
	cell_ids <- cell_ids[order(as.numeric(cell_ids))]	

	start_coords <- seq(from=min_start,to=(max_end - 20e3),by=20e3)
	end_coords <- start_coords + 20e3
	start_and_end_coords_as_character <- paste0(start_coords,"-",end_coords)

	info_subset_wide <- matrix(NA,nrow=length(start_coords),ncol=length(cell_ids),dimnames=list(start_and_end_coords_as_character,cell_ids))

	#Run through the matrix line-by-line and add the relevant data to the big matrix in place of the NAs.
	#Let's say on a given line, you have cell 101, interval 20,000 to 800,000.
	#You would add the appropriate copy number to the matrix for rows named like so:
	#20000-40000
	#40000-60000
	#..
	#780000-800000
	#For column "101".

	#This runs surprisingly fast given the number of lines.

	info_subset_starts <- info_subset$start
	info_subset_ends <- info_subset$end
	info_subset_copynum <- info_subset$copy_num
	info_subset_cells <- info_subset$id

	running_line_count=0

	for(i in 1:nrow(info_subset))
	{
		rownames_this_interval_starts <- seq(from=info_subset_starts[i],to=(info_subset_ends[i] - 20e3),by=20e3)
		rownames_this_interval <- paste0(rownames_this_interval_starts,"-",rownames_this_interval_starts + 20e3)
		info_subset_wide[rownames_this_interval,info_subset_cells[i]] <- info_subset_copynum[i]
		running_line_count = running_line_count + 1
		if((running_line_count %% 1000) == 0){print(paste0("Chromosome ",chr_for_subset," line ",running_line_count," complete."))}
	}

	#Output in CSV format.

	to_cbind <- data.frame(chrom = chr_for_subset,
		start = start_coords,
		end = end_coords,
		stringsAsFactors=FALSE,
		row.names=rownames(info_subset_wide))

	write.csv(cbind(to_cbind,info_subset_wide),
		file=paste0(output_dir,"/copy_number_chrom_",chr_for_subset,".csv"),
		row.names=FALSE,quote=TRUE)
}
