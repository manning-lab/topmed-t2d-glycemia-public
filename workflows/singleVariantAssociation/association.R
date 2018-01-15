# association.R

# Load packages
library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)

# Parse input arguments
input_args <- commandArgs(trailingOnly=T)
gds.file <- input_args[1] 
null.file <- input_args[2]
label <- input_args[3]
test <- input_args[4]
mac <- as.numeric(input_args[5])

# Load nullfile
load(null.file)

# Open gds file
gds.data <- seqOpen(gds.file)

# Filter by desired MAC
seqSetFilter(gds.data,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)
gds.freq <- seqAlleleFreq(gds.data, .progress=TRUE)
gds.maf <- pmin(gds.freq, 1-gds.freq)
gds.mac.filt <- 2 * gds.maf * (1-gds.maf) * length(nullmod$scanID) >= mac

# If no snps remain, return empty
if(sum(gds.mac.filt)==0) {
	print("No SNPs pass MAC filter. Finished Association Step")
	assoc <- NA

# Else move on to association testing
} else {

	# Filter to snps with mac greater than threshold
	seqSetFilter(gds.data, variant.sel=gds.mac.filt, action="intersect", verbose=TRUE)

	# Organize data for output
	pos <- seqGetData(gds.data,"position")
	allele <- seqGetData(gds.data, "allele")
	snps.pos <- cbind(pos,allele)
	
	# Print the number of snps were working with
	print("Filtered SNPs")
	print(dim(snps.pos))

	# Genotype data to the correct format
	gds.geno.data <- SeqVarData(gds.data)

	# Run association test
	assoc <- assocTestMM(genoData = gds.geno.data, nullMMobj = nullmod, test = test)
	print("Finished Association Step")
	print(dim(assoc))
	assoc <- cbind(snps.pos, assoc)

}
## save assoc object
save(assoc, file=paste(label, ".assoc.RData", sep=""))
