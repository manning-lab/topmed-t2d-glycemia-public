library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
sample.file <- input_args[2]
label <- input_args[3]
test <- input_args[4]
pval <- input_args[5]
group.file <- input_args[6]
null.file <- input_args[7]

# load common ID and kinship files
sample.ids <- unique(readLines(sample.file))

## load gds file using GWASTools/GENESIS
geno <- seqOpen(gds)
seqSetFilter(geno,sample.id=sample.ids)
genoData <- SeqVarData(geno)

## load groups
load(group.file)
groups <- groups[!duplicated(names(groups))]

## load null model
load(null.file)

#### run association test
if(test=="SKAT"){
  assoc <- assocTestSeq(genoData, nullmod, groups, test=test, pval.method=pval, weight.beta = c(1,25))
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
}

seqClose(geno)





