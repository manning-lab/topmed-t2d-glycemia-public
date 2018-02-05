# aggAssociation.R

# Load packages
library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)

##### testing inputs
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/results/varshney/genomewide_v1/group_check/freeze4.chr21.pass.gtonly.minDP10.genotypes.gds"
# null.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/null_models/Model15_Sing_F4_sampAFEU_kinsGRM_ALL_cov.age.sex.study_pcNONE.Rda"
# label <- "group_test"
# test <- "SKAT"
# pval <- "kuonen"
# group.file <- "/Users/tmajaria/Documents/projects/topmed/results/varshney/genomewide_v2/groups/groups_16.RData"
# weights <- as.numeric(unlist(strsplit("1,50",",")))
#####

# Parse input args
input_args <- commandArgs(trailingOnly=T)
gds.file <- input_args[1] 
null.file <- input_args[2]
group.file <- input_args[3]
label <- input_args[4]
test <- input_args[5]
pval <- input_args[6]
weights <- as.numeric(unlist(strsplit(input_args[7],",")))


# Load nullfile
load(null.file)

# Open gds file
gds.data <- seqOpen(gds.file)

# Filter to samples in null model
seqSetFilter(gds.data,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)

# Genotype data to the correct format
gds.geno.data <- SeqVarData(gds.data)

## load groups
group_ext <- tail(unlist(strsplit(basename(group.file),'\\.')),n=1)
if (group_ext == 'RData'){
  # load if RData file
  load(group.file)  

} else if (group_ext == 'tsv') {
  # load with data table
  group.raw <- fread(group.file, data.table=F)
  var.df <- data.frame(id = seqGetData(gds.data, "variant.id"), pos = seqGetData(gds.data, "position"), ref = refChar(gds.data), alt = altChar(gds.data))
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

groups = groups[!duplicated(names(groups))]

#### run association test
if(test=="SKAT"){
  assoc <- assocTestSeq(gds.geno.data, nullmod, groups, test=test, pval.method=pval, weight.beta = weights)
  assoc$results = assoc$results[order(assoc$results$pval_0),]
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
} else {
  fwrite(list(), file=paste(label, ".assoc.RData", sep=""))
}

seqClose(gds.data)





