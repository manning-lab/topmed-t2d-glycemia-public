library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
label <- input_args[2]
test <- input_args[3]
pval <- input_args[4]
group.file <- input_args[5]
null.file <- input_args[6]

## load gds file using GWASTools/GENESIS
geno <- seqOpen(gds)
genoData <- SeqVarData(geno)

## load groups
group_ext <- tail(unlist(strsplit(basename(group.file),'\\.')),n=1)
if (group_ext == 'RData'){
  load(group.file)  
  
} else if (group_ext == 'tsv') {
  group.raw <- fread(group.file,sep="\t", data.table=F)
  var.df <- data.frame(id = seqGetData(geno, "variant.id"), pos = seqGetData(geno, "position"), ref = refChar(geno), alt = altChar(geno))
  var.df <- var.df[var.df$pos %in% group.raw$position,]
  group.raw <- group.raw[group.raw$position %in% var.df$pos,]
  group.var <- merge(group.raw, var.df, by.x=c('position','ref','alt'), by.y=c('pos','ref','alt'))
  groups <- list()
  
  for (gid in unique(groups.var$group_id)){
    groups[[gid]] <- groups.var[groups.var$group_id == gid,]
  }
  
} else if (group_ext == 'csv') {
  group.raw <- fread(group.file,sep=",", data.table=F)
  var.df <- data.frame(id = seqGetData(geno, "variant.id"), pos = seqGetData(geno, "position"), ref = refChar(geno), alt = altChar(geno))
  var.df <- var.df[var.df$pos %in% group.raw$position,]
  group.raw <- group.raw[group.raw$position %in% var.df$pos,]
  group.var <- merge(group.raw, var.df, by.x=c('position','ref','alt'), by.y=c('pos','ref','alt'))
  groups <- list()
  
  for (gid in unique(groups.var$group_id)){
    groups[[gid]] <- groups.var[groups.var$group_id == gid,]
  }
  
} else {
  stop("Group file does not have the required extension")
}

groups <- groups[!duplicated(names(groups))]

## load null model
load(null.file)

#### run association test
if(test=="SKAT"){
  assoc <- assocTestSeq(genoData, nullmod, groups, test=test, pval.method=pval, weight.beta = c(1,25))
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
}

seqClose(geno)





