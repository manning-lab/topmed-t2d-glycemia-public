library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
null.file <- input_args[2]
label <- input_args[3]
test <- input_args[4]
minMAC <- 10 # hard coded

#### load null model
load(null.file)

#### run association test
geno <- seqOpen(gds)

#### get the right samples
seqSetFilter(geno,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)

#### filter to only passing sites
var_filt_flag <- data.table(varid=seqGetData(geno,"variant.id"), flag=seqGetData(geno,"annotation/filter"))
keep_var <- var_filt_flag$varid[var_filt_flag$flag == "PASS"]
leave_var <- var_filt_flag$varid[var_filt_flag$flag != "PASS"]
seqSetFilter(geno,variant.id=keep_var, action="intersect", verbose=TRUE)

#### filter by MAF
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
  
  ## dichotomous 
  assoc <- assocTestMM(genoData = genoData, nullMMobj = nullmod, test = test)
  
  print("Finished Association Step")
  print(dim(assoc))
  assoc <- cbind(snps.pos, assoc)
  
}
## save assoc object
save(assoc, file=paste(label, ".assoc.RData", sep=""))
