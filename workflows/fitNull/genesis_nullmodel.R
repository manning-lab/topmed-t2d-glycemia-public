input_args <- commandArgs(trailingOnly=T)
genotype.file <- input_args[1]
phenotype.file <- input_args[2]
outcome.name <- input_args[3]
outcome.type <-  input_args[4]
covariate.string <- input_args[5]
sample.file <- input_args[6]
label <- input_args[7]
kinship.matrix <- input_args[8]
id.col <- input_args[9]

######
# genotype.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr10.pass_and_fail.gtonly.minDP10.gds"
# phenotype.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_AFEU_WesselJ_25AUG2017_T2D.csv"
# kinship.matrix <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/grm/freeze4.autopass.gtonly.minDP10.mmap.grm.fixed.001.matrix.Rda"
# sample.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_AFEU_WesselJ_25AUG2017_T2D_sampleids.txt"
# id.col <- "topmedid"
# label <- "dnanexus_test"
# outcome.name <- "t2d_ctrl"
# outcome.type <- "dichotomous"
# covariate.string <- "sex,last_exam_age,study"
######

suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))

covariates <- unlist(strsplit(covariate.string,","))

## phenotype 
phenotype.data <- fread(phenotype.file,header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)

if (outcome.type == "continuous"){
  phenotype.data[,outcome.name] <- as.numeric(phenotype.data[,outcome.name])
}

for (cur_cov in covariates){
  if (cur_cov %in% c("age_FG", "age", "FG", "FI","PC","age_T2D","age_FI","last_exam_age")){
    phenotype.data[,cur_cov] <- as.numeric(phenotype.data[,cur_cov])
  }
}

phenotype.data <- phenotype.data[!duplicated(phenotype.data[,id.col]),]
all_vals <- c(id.col,outcome.name,covariates)
phenotype.slim <- phenotype.data[,all_vals]

sample.ids <- unique(readLines(sample.file))
phenotype.slim <- phenotype.slim[phenotype.slim[,id.col] %in% sample.ids,na.omit(all_vals,drop=F)]

genotype.data <- seqOpen(genotype.file)
genotype.ids <- seqGetData(genotype.data, "sample.id")

phenotype.slim <- phenotype.slim[phenotype.slim[,id.col] %in% genotype.ids,]

seqSetFilter(genotype.data,sample.id=phenotype.slim[,id.col])
genotype.ids <- seqGetData(genotype.data, "sample.id")
phenotype.slim <- phenotype.slim[match(genotype.ids,phenotype.slim[,id.col]),,drop=F]

seqClose(genotype.data)

kinship <- get(load(kinship.matrix))


kinship <- kinship[row.names(kinship) %in% phenotype.slim[,id.col],colnames(kinship) %in% phenotype.slim[,id.col]]

kinship <- kinship[match(phenotype.slim[,id.col],row.names(kmatr)),match(phenotype.slim[,id.col],colnames(kmatr))]


sample.data <- data.frame(scanID = phenotype.slim[,id.col],  
                          phenotype.slim, 
                          stringsAsFactors=F)
scan.annotated.frame <- ScanAnnotationDataFrame(sample.data)
row.names(phenotype.slim) <- genotype.ids

sample.data.for.annotated <- data.frame(sample.id = genotype.ids,
                                        phenotype.slim,
                                        stringsAsFactors=F)

annotated.frame <- AnnotatedDataFrame(sample.data.for.annotated)

###################
## NULL MODEL
##################
# Should depend on response type

cat('start fit....\n')
kinship = as.matrix(kinship)
cat('Fitting model ')
nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kinship)

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
