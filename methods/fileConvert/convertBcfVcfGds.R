# Adapted from:
# 	Author: topmed analysis pipeline, smgogarten
# 	Link: https://github.com/smgogarten/analysis_pipeline/blob/master/R/vcf2gds.R

library(SeqArray)
library(tools)


input_args <- commandArgs(trailingOnly=T)
infile <- input_args[1]
out_type <- input_args[2]
filetype <- file_ext(infile)

if (identical(out_type, "gds")){
	gds_out <- input_args[3]
	if (filetype == "bcf") {
		infile <- pipe(paste("bcftools view", infile), "rt")
	}
	seqVCF2GDS(infile, gds_out, storage.option="LZMA_RA", verbose=FALSE)
	if (filetype == "bcf") {
		close(infile)
	}

} else {
	vcf_out <- input_args[3]
	tbi_out <- input_args[4]
	vcf_gz_com <- paste("bcftools convert ", infile, " --output-type z --output ",vcf_out,sep="")
	vcf_tbi_com <- paste("tabix -p vcf ", vcf_out, sep="")

	system(vcf_gz_com, ignore.stdout = FALSE, ignore.stderr = FALSE)
	system(vcf_tbi_com, ignore.stdout = TRUE, ignore.stderr = TRUE)

}