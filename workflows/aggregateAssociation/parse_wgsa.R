## usage: trim annotation file to only the categories that you would like
## [1] = wgsa annotation file
## [2] = desired output file (csv format)
## [3] = desired categories (comma separated list)
## [4] = which categories need to be split (might have multiple entries, ie: VEP_ensembl_consequence)

## R --vanilla --args /Users/tmajaria/Documents/projects/topmed/data/varshney/annotation/freezes_2a_3a_4.snp_indel.annotated.general20170422.subset.gz anno.subset.full.csv #chr,pos,ref,alt,VEP_ensembl_Gene_ID,VEP_ensembl_Consequence VEP_ensembl_Gene_ID,VEP_ensembl_Consequence > parse_wgsa.R

args <- commandArgs(trailingOnly=T)
anno.file <- args[1]
out.file <- paste(args[2],".csv",sep="")
# target_columns <- unlist(strsplit(args[3],","))
# columns_to_split <- unlist(strsplit(args[4],","))

# anno.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/annotation/subset/freezes_2a_3a_4.snp_indel.annotated.general20170422.subset.gz.chr7"
# out.file <- "anno.subset.full.csv"
target_columns <- unlist(strsplit("#chr,pos,ref,alt,VEP_ensembl_Gene_ID,VEP_ensembl_Consequence",","))
columns_to_split <- unlist(strsplit("VEP_ensembl_Gene_ID,VEP_ensembl_Consequence",","))

devtools::install_github("UW-GAC/wgsaparsr@2.0.1.4")
library(wgsaparsr)
all_fields <- get_fields(anno.file)
print(all_fields)

# target_columns <- c("#chr", "pos", "ref", "alt", "VEP_ensembl_Gene_ID","VEP_ensembl_Consequence")
# columns_to_split <- c("VEP_ensembl_Gene_ID","VEP_ensembl_Consequence")

parse_to_file(source = anno.file, 
              destination = out.file, 
              desired_columns = target_columns, 
              to_split = columns_to_split, 
              chunk_size = 10000) 

####