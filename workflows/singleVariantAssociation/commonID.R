# Script: commonID.R
# Author: Jasen Jackson, Broad Institute, 2017
# This script takes as input: a genotype file and a 
# phenotype file and returns a list of common IDs. 

# Example use: R --vanilla --args GoT2D.chr22.gds GoT2D.phenotype.ped GoT2DsampleID GoT2D_T2D < commonID.R

library(SeqArray)
library(data.table)

input_args <- commandArgs(trailingOnly=T)
gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
ped <- input_args[2] #"GoT2D.phenotype.ped"
id.col <- input_args[3] # GoT2DSampleID
label <- input_args[4] # Name of common ID text and .RData files

print(paste("GDS file input:",gds))
print(paste("ped file input:",ped))
print(paste("id.column.name:",id.col))
print(paste("label:",label))

# get phenotype IDs
# ped_file <- read.table(ped, header = TRUE, as.is = FALSE)
ped_file <- fread(ped,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
rownames(ped_file) <- ped_file[,id.col] #$GoT2DsampleID

## get genotype IDs
genofile <- seqOpen(gds)
ids <- seqGetData(genofile, "sample.id") #read.gdsn(index.gdsn(genofile,"sample.id"))

#save text file of common IDs
both <- intersect(ids, ped_file[,id.col]) 
write(both, paste(label, ".commonIDs.txt", sep=""), sep="\n")
save(both, file=paste(label,".commonIDs.RData",sep=""))

