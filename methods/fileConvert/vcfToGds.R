# R --vanilla --args GoT2D.chr22.2017714.vcf biallelic.only

library(SeqArray)
input_args <- commandArgs(trailingOnly=T)
vcf <- input_args[1] 
method <- "biallelic.only"

# spl <- unlist(strsplit(vcf,"/",perl=TRUE))
# spl <- spl[length(spl)]
# spl <- unlist(strsplit(spl[length(spl)], "\\.",perl=TRUE))
# gds <- unlist(paste(paste(spl[1:length(spl)-1], collapse = '.'),".gds",sep=""))
# gds <- gds[length(gds)]
gds <- gsub('.vcf','.gds',vcf)


outfile <- seqVCF2GDS(vcf, gds, storage.option="LZMA_RA", verbose=FALSE)

fileConn<-file("output.txt")
writeLines(gds, fileConn)
close(fileConn)
