library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
group.file <- input_args[2]
null.file <- input_args[3]
label <- input_args[4]
test <- input_args[5]
pval <- input_args[6]

if ( length(input_args) > 6 ){
	weights <- as.numeric(unlist(strsplit(input_args[7],",")))
} else {
	weights <- c(1,25)
}

if( !(test %in% c("SKAT", "Burden") ){
	stop("wrong test input")
} else if ( test == "SKAT" & !(pval %in% c("davies", "kuonen","liu") ) {
	stop("wrong pval input")
} else if ( test == "Burden" & !(pval %in% c("Score",  "Wald", "Firth") ) {
	stop("wrong pval input")
} else {
	print("inputs look fine")
}

## load gds file using GWASTools/GENESIS
geno <- seqOpen(gds)
genoData <- SeqVarData(geno)

## load groups
load(group.file)
groups <- groups[!duplicated(names(groups))]

## load null model
load(null.file)

#### run association test
if(test == "SKAT"){
  	assoc <- assocTestSeq(genoData, nullmod, groups, test=test, pval.method=pval, weight.beta = weights)
  	save(assoc, file=paste(label, ".assoc.RData", sep=""))
} else if (test == "Burden") {
	assoc <- assocTestSeq(genoData, nullmod, groups, test=test, burden.test=pval, weight.beta = weights)
	save(assoc, file=paste(label, ".assoc.RData", sep=""))
}

print("Everything seems to have worked fine")
seqClose(geno)