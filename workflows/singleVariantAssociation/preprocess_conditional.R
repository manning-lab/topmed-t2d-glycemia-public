# preprocess_conditional.R
# Description: preprocess the phenotype file for a conditional analysis on some snps, extract the genotype information and add to the phenotypes
# Inputs:
# gds.files : comma separated list of gds file paths, this list must contain gds files for every snp
# phenotype.file : phenotype file to be edited (file, CSV/TSV)
# sample.file : text file with list of sample ids to include, one per line
# snps : comma separated list with the form <chr_number>:<position>,<chr_number>:<position> (string)

# make sure that the right packages are installed
packages <- c("data.table","SeqArray", "SeqVarTools")
to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install,repos='http://cran.us.r-project.org')

# Load packages
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
gds.files <- unlist(strsplit(input_args[1],","))
phenotype.file <- input_args[2]
id.col <- input_args[3]
sample.file <- input_args[4]
snps <- unlist(strsplit(input_args[5],","))
label <- input_args[6]

#### testing inputs
# gds.files <- c("/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/chunk1.freeze4.chrALL.pass.gtonly.minDP10.genotypes.gds")
# phenotype.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_AFEU_WesselJ_25AUG2017_T2D.csv"
# id.col <- "topmedid"
# sample.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_AFEU_WesselJ_25AUG2017_T2D.txt"
# snps <- unlist(strsplit("1:10485,1:10488",","))
####

# get the snps into a dataframe
snp.df <- data.frame()
for (s in snps){
  snp.df <- rbind(snp.df, as.numeric(unlist(strsplit(s,":"))))
}
colnames(snp.df) <- c("chr","pos")

# First get the sample ids
sample.ids <- unique(readLines(sample.file))

# Assign an array to hold the dosage info for each snp
dosage <- list()
alt_ref <- list()

# Loop through the gds files to find each snp, store whether we have found them
for (f in gds.files){
	# open gds file
	f.data <- seqOpen(f)

	# subset by the samples 
	seqSetFilter(f.data, sample.id=sample.ids)

	# Collect some info on the variants so we can see which are in this gds file
	f.var <- data.frame(id=seqGetData(f.data,"variant.id"), chr=seqGetData(f.data,"chromosome"), pos=seqGetData(f.data,"position"))
	
	# subset to only those variants in our snp set
	f.snps <- merge(snp.df, f.var, by.x = c("chr","pos"), by.y = c("chr","pos"))
	
	# if no snps are in this gds file, continue to next file
	if (length(f.snps) == 0) {
	  seqClose(f.data)
	  next
	}
	
	# filter to only the snps in this gds
	seqSetFilter(f.data, variant.id=f.snps$id, action="intersect")
	
	# get the dosage information for those snps
	dosage[[length(dosage)+1]] <- as.data.frame(altDosage(f.data))
	colnames(dosage[[length(dosage)]]) <- paste(f.snps$chr, f.snps$pos, sep=":")
	
	# record the alt and ref alleles
	alt_ref[[length(alt_ref)+1]] <- cbind(as.character(unlist(alt(f.data))), as.character(ref(f.data)))
	
	seqClose(f.data)
}

# combine the dosage info for each snp
dosage = do.call(cbind, dosage)

# Check to make sure that we have all the snps represented, throw and error if now
if (NCOL(dosage) < length(snps)) {
  stop("Not all of your snps were in the GDS files provided, exiting.")
}

# Fix dosage so that sample ids are in their own column
dosage$sample.id <- row.names(dosage)

# combine alt and ref info for each snp
alt_ref = do.call(rbind, alt_ref)
colnames(alt_ref) <- c("alt","ref")
row.names(alt_ref) <- colnames(dosage)

# Load the phenotype file 
phenotype.data <- fread(phenotype.file, data.table=F)

# subset by sample ids
phenotype.data <- phenotype.data[phenotype.data$topmedid %in% sample.ids,]
phenotype.data <- merge(phenotype.data, dosage, by.x = id.col, by.y = "sample.id")

# save new phenotype file
fwrite(phenotype.data, file = paste(label, ".csv", sep=""), sep=",")

# save the alleles
fwrite(alt_ref, file = paste(label, "_alleles.txt", sep=""), sep=" ")