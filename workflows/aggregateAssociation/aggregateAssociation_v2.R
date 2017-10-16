library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
ped <- input_args[2] #"GoT2D.phenotype.ped"
commonID_file <- input_args[3]
id.column.name <- input_args[4]
label <- input_args[5]
outcome <- input_args[6]
outcomeType <- input_args[7]
test <- input_args[8]
pval <- input_args[9]
group.file <- input_args[10]
null.file <- input_args[11]

# load ped_file and re-order columns
ped_file <- read.table(ped, header = TRUE, as.is = FALSE)
rownames(ped_file) <- ped_file[,id.column.name]

# load common ID and kinship files
load(commonID_file)

## trim phenotype file 
ped_file.trimmed <- ped_file[ped_file[,id.column.name] %in% both, ]
ped_file.trimmed <- ped_file.trimmed[both,]


## load gds file using GWASTools/GENESIS
geno <- seqOpen(gds)
genoData <- SeqVarData(geno)

## load groups
load(group.file)

## load null model
load(null.file)

#### run association test
if(test=="SKAT"){
  assoc <- assocTestSeq(genoData, nullmod, groups, test=test, pval.method=pval, weight.beta = c(1,25))
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
}

seqClose(geno)





