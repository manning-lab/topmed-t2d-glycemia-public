# GENESIS-test-v3.R .. accepts pre-configured null-model as input

# INPUT:
# (1) a genotype (gds) file
# (2) a label describing the workflow (for naming the assoc results),
# (3) the outcome type (dichot/cont),  
# (4) the null model (.RData file)

# OUTPUT:
# An assoc.RData file

# command line example
# R --vanilla --args \
#  GoT2D.chr22.biallelic.gds \
#  dichotomous \
#  GoT2D.nullModel.RData

library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
label <- input_args[6]
outcomeType <- input_args[8]
nullmod <- input_args[4]
minMAC <- 30 # hard coded

print(paste("GDS file:",gds))
print(paste("label:",label))
print(paste("outcome type:",outcomeType))
print(paste("null model:", nullmod))

## load null model
load(nullmod)

## load gds file
geno <- seqOpen(gds)

# filter by MAF
seqSetFilter(geno,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)

ref.freq <- seqAlleleFreq(geno, .progress=TRUE)
maf <- pmin(ref.freq, 1-ref.freq)
maf.filt <- 2 * maf * (1-maf) * length(nullmod$scanID) >= minMAC
print(table(maf.filt))

if(sum(maf.filt)==0) {
print("No SNPs pass MAC filter. Finished Association Step")
assoc <- NA

} else {
seqSetFilter(geno, variant.sel=maf.filt, action="intersect", verbose=TRUE)

## add position and rsID
pos <- seqGetData(geno,"position")
#chr.pos <- seqGetData(geno, "$chrom_pos") This line didn't seem to work after filtering geno
allele <- seqGetData(geno, "allele")
snps.pos <- cbind(pos,allele)
print("Filtered SNPs")
print(dim(snps.pos))

genoData <- SeqVarData(geno)
	
	
## dichotomous tests
if(outcomeType=="dichotomous" ) {
	assoc <- assocTestMM(genoData = genoData, nullMMobj = nullmod, test = "Score")
}

## continuous tests
if(outcomeType=="continuous" ) {
	assoc <- assocTestMM(genoData = genoData, nullMMobj = nullmod, test = "Wald")
}
print("Finished Association Step")
print(dim(assoc))
assoc <- cbind(snps.pos, assoc)

}
## save assoc object
save(assoc, file=paste(label, ".assoc.RData", sep=""))
