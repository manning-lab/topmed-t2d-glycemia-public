# This is a ported over copy of the DNANexus app genesis_v0.7 written by Jen Brody 

#== Args 
args<-commandArgs(TRUE)
#===mandatory parameters
phenotype.file <- args[1]
outcome.name <- args[2]
outcome.type <-  args[3]
covariate.string <- args[4]
genotype.files <- args[5]
output.file <- args[6]

#==optional parameters
kinship.matrix <- args[7]
pheno.id <- args[8]


# added these to JSON
test.stat <-  args[9] # Score, Wald, Firth
# test.type  <-  args[10] # Burden, Single, SKAT
conditional <- args[10]
het_vars <-  args[11]



# GLOBAL VARIABLES
collapsing.tests <- c("SKAT",  "Burden")
# test.type.vals <- c("Single","SKAT", "Burden")
test.stat.vals <- c("Score", "Wald", "Firth")

# FROM pipelineFunctions.R

library(plyr)

# This function provides a basic check that the phenotype file is a suitable
# for seqMeta.  It checks that the required column names are in the phenotype
# data frame. 
#
checkPhenotype <- function(p, outcome, covariates, id.col=NULL, gender.col=NULL) {
  if (!is.null(id.col)) {
    if (anyDuplicated(p[ , id.col])) {
      stop("Duplicated phenotype ids.")
    }
  }
  
  missing.covariates <- !(covariates %in% colnames(p))
  if (any(missing.covariates)) {
    msg <- paste("Covariates:", covariates[missing.covariates], "not found in phenotype file.\n", sep=" ")
    print(colnames(p))
    print(covariates %in% colnames(p))
    print(covariates[covariates %in% colnames(p)])
    stop(msg)
  } 
  return(invisible(NULL)) 
}



# This function reduces a data set to only the variables used in a model
# removing subjects with missing data.  Also, it makes the row names of
# the resulting data fram the subject identifier
#
# JAB addition: subsets to complete cases (i.e. no NAs in outcome or covariates)
#
# p: a data frame containing the variables in the model
#
# formula: a character vector which can be coered to an object of class 
#          "formula" (with as.formula): a symbolic description of the model to
#          be fitted. The details of model specification are given under 
#          'Details' in the "lm" help file.
#
# id: (optional) colunm name identifier of the subjects
#
# gender: (optional) colunm name identifier for the gender classification of
#         the subjects.
#
# returns: data frame with only the columns specified in the formula and with
#          the (optional) row names as the subject identifier.
#

reducePheno <- function(pheno.data, 
                        outcome, 
                        covariates = NULL, 
                        hetvars = NULL, 
                        id=NULL, 
                        gender=NULL) {
  checkPhenotype(pheno.data, outcome, covariates, id.col=id, gender.col=gender)   
  if (!is.null(id)) {
    rownames(pheno.data) <- pheno.data[ ,id]
  }
  
  all.terms <- unique(c(outcome, covariates, hetvars, gender))
  cat('all terms',print(all.terms),'\n')
  pheno.data <- as.data.frame(pheno.data) 
  pheno <- na.omit(pheno.data[, all.terms, drop=F])
  return(pheno)
}

# Calculate MAF
#
# dose: matrix with dosages rows individuals, columns are variants


Maf <- function(dose){
       aaf <- colMeans(dose,na.rm=T)/2
       return(min(1-aaf, aaf))
}


split.by.comma <- function(cur.string){
    cur.string <- gsub('"', '', cur.string)
    out <- unlist(strsplit(cur.string, ","))
    if (length(out) == 0){
        out = NULL
    }
    return(out)
}
    


filterByMAF <- function(gds, sample.id=NULL, mac.min=NA, maf.min=NA, verbose=TRUE) {
    if ((!is.na(mac.min) & mac.min > 1) |
        (!is.na(maf.min) & maf.min > 0)) {
        if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
        seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
        ref.freq <- seqAlleleFreq(gds)
        maf <- pmin(ref.freq, 1-ref.freq)
        if (!is.na(mac.min)) {
            maf.filt <- 2 * maf * (1-maf) * length(sample.id) >= mac.min
            if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAC >=", mac.min))
        } else {
            maf.filt <- maf >= maf.min
            if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAF >=", maf.min))
        }
        seqSetFilter(gds, variant.sel=maf.filt, action="intersect", verbose=verbose)
    }
}


GetFamilyDistribution <- function(response.type) {
               if (tolower(response.type) == "continuous"){
                      family = "gaussian"
               } else if (tolower(response.type) == "dichotomous"){
                      family = "binomial"
               } else {
                      msg = paste("Don't know how to deal with response type", response.type)
                      stop(msg)
               }
               return(family)
           }

GetKinshipMatrix <- function(kinship.matrix){
  cat('Loading Kinship Matrix:',kinship.matrix,'\n')
  if(grepl('Rda',kinship.matrix,ignore.case=TRUE)){
    kmatr = get(load(kinship.matrix))
  }
  else{
    kmatr = as.matrix(read.csv(kinship.matrix,as.is=T,check.names=F,row.names=1))
  }

  cat('Loaded Kinship NROW:',NROW(kmatr),' NCOL:',NCOL(kmatr),'\n')
  kmatr
}

cat('output.file',output.file,'\n')
cat('kinship.matrix',kinship.matrix,'\n')
cat('test.stat',test.stat,'\n')
# cat('test.type',test.type,'\n')
cat('outcome.type',outcome.type,'\n')
cat('het_vars',het_vars,'\n')
cat('conditional',conditional,'\n')

if(conditional != 'NA'){
  cpos = strsplit(conditional,':')[[1]][2]
}else{
  cpos = FALSE
}
cat('conditional',conditional,'\t pos',cpos,'\n')

if (!(test.stat %in% test.stat.vals)){
     msg = paste("The requested test statistic:", test.stat, "is not available (Use Firth, Score, Wald!")
     stop(msg)
}


# if(!test.type %in% test.type.vals){
    # stop("Test type must be one of ",paste(test.type.vals,sep=','))
# }


suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GWASTools))
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(gdsfmt))
suppressMessages(library(bdsmatrix))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))
#suppressMessages(library(parallel))
## Setup
source("/home/dnanexus/pipelineFunctions.R")
covariates <- split.by.comma(covariate.string)  


## phenotype 
phenotype.data <- read.csv(phenotype.file, header=TRUE, as.is=TRUE)
if(NCOL(phenotype.data) < 2){
    
     msg = paste("Is the phenotype file a CSV?  Too few columns from read.csv()")
     warning(msg)
}


cat('Input pheno N=',nrow(phenotype.data),'\n')
if(het_vars != 'NA'){
    cat('prep pheno with het vars')
    pheno <- reducePheno(phenotype.data, outcome.name, covariates=covariates,hetvars=het_vars, id=pheno.id)
}else{
    cat('prep pheno without het vars\n')
    het_vars = NA
    pheno <- reducePheno(phenotype.data, outcome.name, covariates=covariates, id=pheno.id)
}
cat('Output pheno N=',nrow(pheno),'\n')

## Report dropped individuals
dropped.ids.selector <- !(phenotype.data[[pheno.id]] %in% row.names(pheno))
dropped.ids <- phenotype.data[[pheno.id]][dropped.ids.selector] 
if (NROW(dropped.ids) != 0 ) {
  cat("Dropped because of incomplete cases:", length(dropped.ids) )
}

# For GDS files
f <- seqOpen(genotype.files)
sample.ids <- seqGetData(f, "sample.id")
all.terms <- unique(c(outcome.name, covariates, het_vars))
pheno <- pheno[row.names(pheno) %in% sample.ids,na.omit(all.terms),drop=F]
cat('Output pheno after mergeing with Genos N=',nrow(pheno),'\n')
if(nrow(pheno) == 0){
    msg = paste("Phenotype ID column doesn't match IDs in GDS")
    stop(msg)
}

full.sample.ids <- sample.ids 

#subset to phenotyped samples
seqSetFilter(f,sample.id = row.names(pheno))

# order pheno to the GDS subject order
sample.ids <- seqGetData(f, "sample.id")
pheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]



## Conditional analaysis
if(cpos >0){
    
    pos = seqGetData(f, "position")
  cat('Conditioning on ',conditional,'...\n')
  cidx = which(pos == as.numeric(cpos))
  if(any(cidx)){
    
    seqSetFilter(f,variant.sel=cidx, sample.id = row.names(pheno),verbose=FALSE)
    pheno$csnp = altDosage(f,  use.names=FALSE)
  }else{
    stop('Can not find snp ',conditional,' with position ',cpos,' to condition on in data file')
  }
  
  dropConditionalCases = NROW(pheno)-NROW(pheno[complete.cases(pheno),])
  if(dropConditionalCases > 0){
    cat('Warning: Dropping ',dropConditionalCases,' samples due to missing conditional genotype calls\n')
  }

  pheno = pheno[complete.cases(pheno),]
  
  covariates[length(covariates) + 1] <- 'csnp'
}




## Load KINSHIP matrix
## Kinship doesn't contain all samples
kmatr = GetKinshipMatrix(kinship.matrix)
pheno = pheno[row.names(pheno) %in% row.names(kmatr),,drop=F]
kmatr = kmatr[row.names(kmatr) %in% row.names(pheno),colnames(kmatr) %in% row.names(pheno)]
cat('Output pheno in Kinship N=',nrow(pheno),'\n')
kmatr = kmatr[match(row.names(pheno),row.names(kmatr)),match(row.names(pheno),colnames(kmatr))]
if(nrow(pheno) == 0){
    msg = paste("Phenotype ID column doesn't match IDs in Kinship Matrix")
    stop(msg)
}


# Get sample ids to check order 
seqSetFilter(f,sample.id = row.names(pheno))
sample.ids <- seqGetData(f, "sample.id")

if (!(identical(sample.ids,row.names(pheno)) && identical(row.names(kmatr),row.names(pheno)))){
        stop("Something is off problem with re-ordering")
}

seqClose(f)

sample.data <- data.frame(scanID = row.names(pheno),  
                    pheno, 
                    stringsAsFactors=F)
scan.annotated.frame <- ScanAnnotationDataFrame(sample.data)
modified.pheno = pheno[full.sample.ids,,drop=FALSE]
row.names(modified.pheno) <- full.sample.ids

sample.data.for.annotated <- data.frame(sample.id = full.sample.ids,
                                        modified.pheno,
                                        stringsAsFactors=F)
rm(modified.pheno)

annotated.frame <- AnnotatedDataFrame(sample.data.for.annotated)

###################
## NULL MODEL
##################
# Should depend on response type

cat('start fit....\n')
kmatr = as.matrix(kmatr)
if (test.stat == 'Firth'){
  cat('WARNING: Firth test does NOT use kinship information - unrelated only')
  nullmod <- fitNullReg(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type))

}else if (!is.na(het_vars)){
  cat('Fitting model with heterogeneous variances')
  nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     group.var = het_vars,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kmatr)
}else{
  cat('Fitting model ')
  nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kmatr)
}
save(nullmod,annotated.frame,file=paste(output.file,".RDa",sep="")
