# INPUT:
# (1) a genotype (gds) file
# (2) a phenotype (ped) file
# (3) a GRM kinship.objects.RData file (from kinshipMatrix.R file)
# (4) a commonIDs .RData file
# (5) a chromatin state file (individual tissue)
# (6) an annotation file (for filtering variants)
# (7) annotation type to filter on (category in annotation)
# (8) annotation value to filter on (string)
# (9) the name of the column containing the IDs in the ped file
# (10) a label describing the workflow (for naming the assoc results),
# (11) the outcome under study
# (12) the outcome type (dichot/cont), 
# (13) chromatin states to include (, separated list, must match chromatin state file)
# (14) source file for group operations
# (15) type of aggregate test to perform (skat, burden)
# (16) pvalue type for burden test (score, wald, firth) or (kuonen, davies, liu)
# (16) covariates to use, separated by "," and no space. Example: "age,sexTrait,BMI"   

# OUTPUT:
# An assoc.RData file

# command line example
# R --vanilla --args \
#  GoT2D.chr22.biallelic.gds \
#  GoT2D.phenotype.ped \
#  GoT2D_T2D.kinship.objects.RData \ 
#  GoT2D_T2D.commonIDs.RData \
#  GoT2DsampleID \
#  GoT2D_T2D \
#  T2D \ 
#  dichotomous \
#  age sex BMI 
# 
library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
ped <- input_args[2] #"GoT2D.phenotype.ped"
# kinship_file <- input_args[3] # GoT2DSampleID
commonID_file <- input_args[4]
# chr_st.file <- inputs_args[5]
# anno.file <- input_args[6]
# anno.type <- input_args[7]
# anno.value <- input_args[8]
id.column.name <- input_args[9]
label <- input_args[10]
outcome <- input_args[11]
outcomeType <- input_args[12]
# minMAC <- 30 # hard coded
# chr_st.names <- unlist(strsplit(input_args[13],split=","))
# source.file <- input_args[14]
test <- input_args[15]
pval <- input_args[16]
# gene.file <- input_args[17]
# covariates <- NULL

# print(paste("GDS file:",gds))
print(paste("Ped file:",ped))
print(paste("Kinship file:",kinship_file))
print(paste("Common IDs file:",commonID_file))
print(paste("ID column name:",id.column.name))
print(paste("label:",label))
print(paste("outcome:",outcome))
print(paste("outcomeType:",outcomeType))


if(length(input_args)>17) {
  covariates <- strsplit(input_args[17],split=",")[[1]]
  print(paste("covariates:",paste(covariates,collapse=" ")))
}


### testing inputs ###

# gds <- "/Users/tmajaria/Documents/projects/topmed/code/testing_inputs/singleVariantFull/freeze4.chr21.pass.gtonly.minDP10.genotypes.gds"
# ped <- "/Users/tmajaria/Documents/projects/topmed/code/testing_inputs/singleVariantFull/Pooled_T2D_Traits_FHS_JHS_CFS_Amish_freeze4_pcair_pcs.ped"
# kinship_file <- "/Users/tmajaria/Documents/projects/topmed/code/testing_inputs/singleVariantFull/freeze4_round2_pcrelate.gds"
# commonID_file <- "/Users/tmajaria/Documents/projects/topmed/code/testing_inputs/singleVariantFull/freeze4.chr21.pass.gtonly.minDP10.genotypes.commonIDs.RData"
# chr_st.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/data/local_cs_states/Islets.chromatinStates.bed"
# anno.file <- ""

# id.column.name <- "TOPMEDID"
# label <- "group_test"
# outcome <- "T2D"
# outcomeType <- "dichotomous"
# minMAC <- 30 # hard coded
# chr_st.names <- c("1_Active_TSS")
# source.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/code/get_agg_units.R"
# covariates <- strsplit("T2D_age,sex",split=",")[[1]]
# print(paste("covariates:",paste(covariates,collapse=" ")))
# test <- "burden"



# load ped_file and re-order columns
ped_file <- read.table(ped, header = TRUE, as.is = FALSE)
rownames(ped_file) <- ped_file[,id.column.name]
head(ped_file)

# load common ID and kinship files
load(commonID_file)

library(gdsfmt)
geno.pcrelate <- openfn.gds(kinship_file)
# load(paste(kinship_file)) #ped_file.trimmed, PCRgrm, filteredSNPs, both

## trim phenotype file 
ped_file.trimmed <- ped_file[ped_file[,id.column.name] %in% both, ]
ped_file.trimmed <- ped_file.trimmed[both,]


## load gds file using GWASTools/GENESIS
#geno <- GdsGenotypeReader(filename=gds)
#genoData <- GenotypeData(geno)

# convert phenotype dataframe to ScanAnnotationDataFrame object
if(!is.null(covariates)) {
  df <- cbind(setNames(data.frame(both,ped_file.trimmed[,outcome]), c("scanID",outcome)), covs = ped_file.trimmed[,covariates])
  
  # change covariates = [age, sex, BMI] to covariates = [covs.age, covs.sex, covs.BMI]
  covariates <- paste("covs.", covariates, sep="")
  
  # double check sex column;
  if("covs.sex" %in% covariates) {
    print("Sex is one of the covariates")
    if(sum(c(1,2) %in% names(table(df[,"covs.sex"])))==2) {
      df$sex <- ifelse(df[,"covs.sex"]==2,"F","M")
      print("Converted sex from 1/2 to M/F")
      print(table(sex=df$sex,covs.sex=df$covs.sex,useNA="always"))
    } else {
      df$sex <- df$covs.sex
      print(table(sex=df$sex,covs.sex=df$covs.sex,useNA="always"))
    }
    covariates <- c(setdiff(covariates,"covs.sex"),"sex")
    print(paste("covariates:",paste(covariates,collapse=" ")))	
  }
} else {
  df <- setNames(data.frame(both,ped_file.trimmed[,outcome]), c("scanID",outcome))
  
}

scanAnnot <- ScanAnnotationDataFrame(df)
str(scanAnnot)

#### Null model - Future Versions of this script can be extended to input null models

# get GRM in correct format
PCRgrm <- pcrelateMakeGRM(geno.pcrelate, scan.include = both, scaleKin = 2)

if(outcomeType=="dichotomous" & is.null(covariates)) {
  nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = binomial, covMatList = PCRgrm)
}
if(outcomeType=="dichotomous" & !is.null(covariates)) {
  nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = binomial, covMatList = PCRgrm, covars = covariates)
}
if(outcomeType=="continuous" & !is.null(covariates)) {
  nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = gaussian, covMatList = PCRgrm, covars = covariates)
}

if(outcomeType=="continuous" & is.null(covariates)) {
  nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = gaussian, covMatList = PCRgrm)
}


#### run association test
geno <- seqOpen(gds)

# filter by MAF
# seqSetFilter(geno,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)

# ref.freq <- seqAlleleFreq(geno, .progress=TRUE)
# maf <- pmin(ref.freq, 1-ref.freq)
# maf.filt <- 2 * maf * (1-maf) * length(nullmod$scanID) >= minMAC
# print(table(maf.filt))

# if(sum(maf.filt)==0) {
#   print("No SNPs pass MAC filter. Finished Association Step")
#   assoc <- NA
#   
# } else {
  # seqSetFilter(geno, variant.sel=maf.filt, action="intersect", verbose=TRUE)
  
  ## add position and rsID
  # pos <- seqGetData(geno,"position")
  #chr.pos <- seqGetData(geno, "$chrom_pos") This line didn't seem to work after filtering geno
  # allele <- seqGetData(geno, "allele")
  # snps.pos <- cbind(pos,allele)
  # print("Filtered SNPs")
  # print(dim(snps.pos))
  
  genoData <- SeqVarData(geno)
  
  # get the groups
  
  # variant.groups <- groupByGene(geno.gds,chr_st.file,anno.file,gene.file,chr_st.names,anno.value,minmaf)
  load("/Users/tmajaria/Documents/projects/topmed/results/varshney/v1/groups_v1.RData")
  
  groups <- groups[!sapply(groups, is.null)]
 
  ## dichotomous 
  if(outcomeType=="dichotomous" ) {
    if (test=="burden"){
      assoc <- assocTestSeq(genoData, nullmod, groups, test=test, burden.test=pval)
    } else if (test == "skat"){
      assoc <- assocTestSeq(genoData, nullmod, groups, test=test, pval.method=pval)
    } else {
      assoc <- c()
      print("sorry, that didnt work")
    }
    # assoc  <- assocTestMM(genoData = genoData, nullMMobj = nullmod, test = "Score")
  }
  
  ## continuous 
  if(outcomeType=="continuous" ) {
    if (test=="Burden"){
      assoc <- assocTestSeq(genoData, nullmod, groups, test=test, burden.test=pval)
    } else if (test == "SKAT"){
      assoc <- assocTestSeq(genoData, nullmod, groups, test=test, pval.method=pval,weight.beta <- c(1,25))
    } else {
      assoc <- c()
      print("sorry, that didnt work")
    }
  }
  print("Finished Association Step")
  print(dim(assoc))
  # assoc <- cbind(snps.pos, assoc)
  
# }
## save assoc object
  seqClose(gds)
save(assoc, file=paste(label, ".assoc.RData", sep=""))





