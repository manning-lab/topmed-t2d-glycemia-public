# Adapted from:
# 	Author: topmed analysis pipeline, smgogarten
# 	Link: https://github.com/smgogarten/analysis_pipeline/blob/master/R/vcf2gds.R

library(SeqArray)
library(tools)


input_args <- commandArgs(trailingOnly=T)
infile <- input_args[1]
gds_flag <- input_args[2]
vcf_flag <- input_args[3]

filetype <- file_ext(infile)

#outfiles <- c()

if (vcf_flag == "true"){
	vcf_out <- gsub(filetype,'vcf.gz',infile)
	vcf_out_tbi <- gsub(filetype,'vcf.gz.tbi',infile)
	#outfiles <- c(outfiles,vcf_out,vcf_out_tbi)
	vcf_gz_com <- paste("bcftools convert ", infile, " --output-type z --output ",vcf_out,sep="")
	vcf_tbi_com <- paste("tabix -p vcf ", vcf_out, sep="")

	system(vcf_gz_com, ignore.stdout = FALSE, ignore.stderr = FALSE)
	system(vcf_tbi_com, ignore.stdout = TRUE, ignore.stderr = TRUE)
}

if (gds_flag == "true"){
	gds_out <- gsub(filetype,'gds',infile)
	#outfiles <- c(outfiles,gds_out)
	if (filetype == "bcf") {
		infile <- pipe(paste("bcftools view", infile), "rt")
	}
	seqVCF2GDS(infile, gds_out, storage.option="LZMA_RA", verbose=FALSE)
	if (filetype == "bcf") {
		close(infile)
	}
}

print(gds_out)
print(vcf_out)
print(vcf_out_tbi)

fileConn<-file("output_gds.txt")
writeLines(paste(gds_out, fileConn)
close(fileConn)

fileConn<-file("output_vcf.txt")
writeLines(vcf_out, fileConn)
close(fileConn)

fileConn<-file("output_tbi.txt")
writeLines(vcf_out_tbi, fileConn)
close(fileConn)









       