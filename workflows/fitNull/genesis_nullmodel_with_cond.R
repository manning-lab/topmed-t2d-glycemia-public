args<-commandArgs(TRUE)
phenotype.file <- args[1]
outcome.name <- args[2]
outcome.type <-  args[3]
covariate.string <- args[4]
sample.file = args[5]
label <- args[6]
kinship.matrix <- args[7]
pheno.id <- args[8]
conditional <- args[9]
genotype.files <- args[10]

# genotype.files <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr10.pass_and_fail.gtonly.minDP10.SUBSET.1000000.gds"
# phenotype.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_Glycemic_Traits_freeze5b_TDM_12062017_FG.ped"
# kinship.matrix <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/grm/freeze.5b.auto.pass.gtonly.minDP10.mmap.grm.DP.fixed.0.001.matrix.cor.Rda"
# sample.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_Glycemic_Traits_freeze5b_TDM_12062017_FG_sampleids.txt"
# pheno.id <- "TOPMEDID"
# label <- "cond_test"
# outcome.name <- "FastingGlucose"
# outcome.type <- "Continuous"
# conditional <- "A:6687759"
# covariate.string <- "sex,age_FG,STUDY_ANCESTRY"

library(plyr)
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
  return(list(pheno,all.terms))
}

split.by.comma <- function(cur.string){
  cur.string <- gsub('"', '', cur.string)
  out <- unlist(strsplit(cur.string, ","))
  if (length(out) == 0){
    out = NULL
  }
  return(out)
}

GetFamilyDistribution <- function(response.type) {
  if (response.type == "continuous"){
    family = "gaussian"
  } else if (response.type == "dichotomous"){
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

if(conditional != 'NA'){
  cpos = strsplit(conditional,':')[[1]][2]
}else{
  cpos = FALSE
}
cat('conditional',conditional,'\t pos',cpos,'\n')

suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GWASTools))
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(gdsfmt))
suppressMessages(library(bdsmatrix))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))

covariates <- split.by.comma(covariate.string)  

## phenotype 
phenotype.data <- fread(phenotype.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
if (NCOL(phenotype.data) < 2){
  phenotype.data <- fread(phenotype.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
}

phenotype.data <- phenotype.data[!duplicated(phenotype.data[,1]),]
pa <- reducePheno(phenotype.data, outcome.name, covariates=covariates, id=pheno.id)
pheno <- pa[[1]]
all.terms <- pa[[2]]

dropped.ids.selector <- !(phenotype.data[[pheno.id]] %in% row.names(pheno))
dropped.ids <- phenotype.data[[pheno.id]][dropped.ids.selector] 
if (NROW(dropped.ids) != 0 ) {
  cat("Dropped because of incomplete cases:", length(dropped.ids) )
}

sample.ids <- unique(readLines(sample.file))
pheno <- pheno[row.names(pheno) %in% sample.ids,na.omit(all.terms,drop=F)]

if(nrow(pheno) == 0){
  msg = paste("Phenotype ID column doesn't match IDs in GDS")
  stop(msg)
}

if(cpos >0){
  cat('Conditioning on ',conditional,'...\n')
  f <- seqOpen(genotype.files)
  pos = seqGetData(f, "position")
  variant.ids = seqGetData(f, "variant.id")
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

  seqClose(f)
}


## Load KINSHIP matrix
## Kinship doesn't contain all samples
kmatr = GetKinshipMatrix(kinship.matrix)
pheno = pheno[row.names(pheno) %in% row.names(kmatr),,drop=F]
kmatr = kmatr[row.names(kmatr) %in% row.names(pheno),colnames(kmatr) %in% row.names(pheno)]
kmatr = kmatr[match(row.names(pheno),row.names(kmatr)),match(row.names(pheno),colnames(kmatr))]
if(nrow(pheno) == 0){
  msg = paste("Phenotype ID column doesn't match IDs in Kinship Matrix")
  stop(msg)
}



sample.data <- data.frame(scanID = row.names(pheno),  
                          pheno, 
                          stringsAsFactors=F)
scan.annotated.frame <- ScanAnnotationDataFrame(sample.data)
modified.pheno = pheno[sample.ids,,drop=FALSE]
row.names(modified.pheno) <- sample.ids

sample.data.for.annotated <- data.frame(sample.id = sample.ids,
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
cat('Fitting model ')
nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kmatr)

save(nullmod,annotated.frame,file=paste(label,"_null.RDa",sep=""))

library("ggplot2")
if (tolower(outcome.name) == "t2d"){
  
  cov_to_remove <- c("study_ancestry","STUDY_ANCESTRY","sex")
  quant_covars <- covariates[! covariates %in% cov_to_remove]
  
  pheno <- pheno[!is.na(pheno[,outcome.name]),]
  case = pheno[pheno[,outcome.name] == 2 | pheno[,outcome.name] == 1,]
  control = pheno[pheno[,outcome.name] == 0,]

  male.val <- "M"
  female.val <- "F"

  m = pheno[pheno$sex==male.val,]
  f = pheno[pheno$sex==female.val,]

  if (NCOL(m)==0 && NCOL(f)==0){
    male.val <- 1
    female.val <- 2
  }

  m_case = case[case$sex==male.val,]
  m_control = control[control$sex==male.val,]
  f_case = case[case$sex==female.val,]
  f_control = control[control$sex==female.val,]

  all_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  case_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  ctrl_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  m_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  f_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  mcase_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  mctrl_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  fcase_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  fctrl_stats <- data.frame(matrix(NA, nrow = 3, ncol = length(quant_covars)))
  
  colnames(all_stats) <- quant_covars
  colnames(case_stats) <- quant_covars
  colnames(ctrl_stats) <- quant_covars
  colnames(m_stats) <- quant_covars
  colnames(f_stats) <- quant_covars
  colnames(mcase_stats) <- quant_covars
  colnames(mctrl_stats) <- quant_covars
  colnames(fcase_stats) <- quant_covars
  colnames(fctrl_stats) <- quant_covars

  row.names(all_stats) <- c("mean","median","sd")
  row.names(case_stats) <- c("mean","median","sd")
  row.names(ctrl_stats) <- c("mean","median","sd")
  row.names(m_stats) <- c("mean","median","sd")
  row.names(f_stats) <- c("mean","median","sd")
  row.names(mctrl_stats) <- c("mean","median","sd")
  row.names(mcase_stats) <- c("mean","median","sd")
  row.names(fcase_stats) <- c("mean","median","sd")
  row.names(fctrl_stats) <- c("mean","median","sd")

  if ("STUDY_ANCESTRY" %in% covariates){
    study_ancestry_val <- "STUDY_ANCESTRY"
  } else {
    study_ancestry_val <- "study_ancestry"
  }

  pdf(paste(label,"_plots.pdf",sep=""),width=11)
  layout(matrix(seq(1,6*(length(quant_covars)-1)),nrow=length(quant_covars)-1,ncol=6,byrow=T))
  for (i in quant_covars){
    print(i)
    all_stats[1,i] <- mean(pheno[,i])
    all_stats[2,i] <- median(pheno[,i])
    all_stats[3,i] <- sd(pheno[,i])
    
    case_stats[1,i] <- mean(case[,i])
    case_stats[2,i] <- median(case[,i])
    case_stats[3,i] <- sd(case[,i])
    
    ctrl_stats[1,i] <- mean(control[,i])
    ctrl_stats[2,i] <- median(control[,i])
    ctrl_stats[3,i] <- sd(control[,i])
    
    m_stats[1,i] <- mean(m[,i])
    m_stats[2,i] <- median(m[,i])
    m_stats[3,i] <- sd(m[,i])
    
    f_stats[1,i] <- mean(f[,i])
    f_stats[2,i] <- median(f[,i])
    f_stats[3,i] <- sd(f[,i])
    
    mcase_stats[1,i] <- mean(m_case[,i])
    mcase_stats[2,i] <- median(m_case[,i])
    mcase_stats[3,i] <- sd(m_case[,i])
    
    mctrl_stats[1,i] <- mean(m_control[,i])
    mctrl_stats[2,i] <- median(m_control[,i])
    mctrl_stats[3,i] <- sd(m_control[,i])
    
    fcase_stats[1,i] <- mean(f_case[,i])
    fcase_stats[2,i] <- median(f_case[,i])
    fcase_stats[3,i] <- sd(f_case[,i])
    
    fctrl_stats[1,i] <- mean(f_control[,i])
    fctrl_stats[2,i] <- median(f_control[,i])
    fctrl_stats[3,i] <- sd(f_control[,i])

    plot <- ggplot(pheno, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_violin() + labs(title = "All samples"))
    
    
    plot <- ggplot(case, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_violin() + labs(title = "Case samples"))
    
    
    plot <- ggplot(control, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_violin() + labs(title = "Control samples"))
    
    
    plot <- ggplot(pheno, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_violin(aes(fill=factor(sex))) + labs(title = "All samples by sex"))
    
    
    plot <- ggplot(case, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_violin(aes(fill=factor(sex))) + labs(title = "Case samples by sex"))
    
    
    plot <- ggplot(control, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_violin(aes(fill=factor(sex))) + labs(title = "Control samples by sex"))
    
  }

  dev.off()

  all_stat_frames <- list(all_stats,case_stats,ctrl_stats,m_stats,f_stats,mcase_stats,mctrl_stats,fcase_stats,fctrl_stats)
  names(all_stat_frames) <- c("all","case","control","all_male","all_female","male_case","male_control","female_case","female_control")

  for (i in seq(1,length(all_stat_frames))){
    fwrite(all_stat_frames[[names(all_stat_frames)[i]]], file=paste(label,names(all_stat_frames)[i],"stats.csv",sep="_"))
  }

} else {
  pdf(paste(label,"_plots.pdf",sep=""),width=11)
  dev.off()
  fwrite(data.frame(na="NA"), file=paste(label,"_stats.csv",sep=""))
}
