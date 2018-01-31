
##inverse normal function
inversenormal <- function(x) {
  # inverse normal if you have missing data
  return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
}

##make residuals transform function

make.residuals.transform <- function(outcome,covariates,pheno,by.sex=F,transform=identity,controls.only=F,label="_adj",covariates.cases=NA,return.raw=T,casesTF=NA,femaleTF=NA) {
  
  if(is.na(covariates.cases)) covariates.cases=covariates
  
  residuals <- rep(NA,nrow(pheno)) 
  
  if(by.sex==F) {
    
    ## case specific
    residuals[which(!casesTF)] <- do.lm(formula.2=as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF)
    
    if(controls.only==F) {
      residuals[which(casesTF)] <- do.lm(formula.2=as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=casesTF)
    }  
    
    residual.name <- paste(outcome,label,sep="")
    
    ## apply transformation to residual
    residuals_invn <- transform(residuals)
    
    if(return.raw) {
      to.return <- cbind(pheno[,outcome],residuals)
      
      colnames(to.return) <- c(outcome,residual.name)
      
    } else {
      to.return <- cbind(residuals, residuals_invn)
      colnames(to.return) <- c(residual.name, paste(outcome,label,"_invn",sep=""))
    }
    
    
    return(to.return)
    
    
  } else {
    ## sex-specific residuals and transformations  
    residuals.males <- rep(NA,nrow(pheno))
    residuals.females <- rep(NA,nrow(pheno))
    if(controls.only==F){
      print("controls,males")
      residuals[!casesTF & !femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & !femaleTF)
      residuals.males[!casesTF & !femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & !femaleTF)
      
      print("controls,females")
      residuals[!casesTF & femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & femaleTF)
      residuals.females[!casesTF & femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & femaleTF)
      
      print("cases,males")
      residuals[casesTF & !femaleTF] <- do.lm(as.formula(paste(outcome,covariates.cases,sep="~")),pheno.2=pheno,subset.2=casesTF & !femaleTF)
      residuals.males[casesTF & !femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=casesTF & !femaleTF)
      
      print("cases,females")
      residuals[casesTF & femaleTF] <- do.lm(as.formula(paste(outcome,covariates.cases,sep="~")),pheno.2=pheno,subset.2=casesTF & femaleTF)
      residuals.females[casesTF & femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=casesTF & femaleTF)
      
    }else{
      
      print("controls,males")
      residuals[!casesTF & !femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & !femaleTF)
      residuals.males[!casesTF & !femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & !femaleTF)
      
      print("controls,females")
      residuals[!casesTF & femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & femaleTF)
      residuals.females[!casesTF & femaleTF] <- do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & femaleTF)
    }
    ## residuals, not inverse normalized
    
    residual.name <- paste(outcome,label,"_by_sex",sep="")
    
    to.return1 <- cbind(residuals,residuals.males,residuals.females)
    colnames(to.return1) <- c(residual.name,paste(outcome,label,c("_males","_females"),sep=""))
    
    ## apply transformation
    residual.name <- paste(outcome,label, "_invn", sep="")
	
	##does not transform:
    residuals_invn <- transform(residuals)
    residuals.females_invn <- transform(residuals.females)
    residuals.males_invn  <- transform(residuals.males)
    to.return2 <- cbind(residuals_invn, residuals.females_invn, residuals.males_invn)
	######tried 3 different methods for calling transform with the inversenormal function, not have worked. Placed transform in the for loop instead for the time being.  
 
    
	colnames(to.return2) <- c(residual.name,paste(outcome,label,c("_invn_males","_invn_females"),sep=""))
    
    to.return <-cbind(to.return1, to.return2)
    return(to.return)
  }
  
}

##do.lm function
do.lm <- function(formula.2, pheno.2, subset.2) {
  print(dim(pheno.2[which(subset.2),]))
  print(summary(lm(formula.2,data=pheno.2[which(subset.2),],na.action="na.exclude")))
  residuals(lm(formula.2,data=pheno.2[which(subset.2),],na.action="na.exclude"))
}

##prepare files
in.d <- subset(in.d, FastingInsulin >0) ## to remove FastingInsulin of 0. 
in.d$logFastingInsulin <- log(in.d$FastingInsulin)

pooled.traits <- in.d
#pooled.traits <- in.d[,c("Individual_ID","Family_ID","Mother_ID","Father_ID","sex","T2D","STUDY_ANCESTRY","age_FI","T2D_FI","BMI_FI","FastingInsulin","logFastingInsulin","FastingGlucose")]
cohort.traits <- in.d
#cohort.traits <- in.d[,c("Individual_ID","Family_ID","Mother_ID","Father_ID","sex","T2D","STUDY_ANCESTRY","age_FI", "T2D_FI","BMI_FI","FastingInsulin","logFastingInsulin","FastingGlucose")]

pooled.traits$logFastingInsulin_adj <- NA
pooled.traits$logFastingInsulin_adj_invn <- NA

in.d.cases.FI <- !in.d$T2D_FI %in% c(0,1)

in.d.females <- in.d$sex==2
table(in.d.females)

cohorts <- names(table(in.d$STUDY_ANCESTRY))

####WHI only includes one sex, which breaks the for loop as written. Made separate loop just for WHI cohorts to run make.residuals function by.sex=F, but treat as by.sex=T when pooling data #####
remove <- c("WHI_AA", "WHI_AS", "WHI_EU", "WHI_HA")
WHI <- cohorts[cohorts %in% remove]
cohorts <- cohorts[!cohorts %in% remove]

sink("make.residuals.txt")

for(cohort in cohorts) {
  print(paste("COHORT:",cohort))
  
  in.d[which(in.d$STUDY_ANCESTRY==cohort),]
  
  print(paste("FASTING INSULIN"))
  tmp.FI <- make.residuals.transform("logFastingInsulin","age_FI+BMI_FI", controls.only=T, by.sex=T,
                                     in.d[which(in.d$STUDY_ANCESTRY==cohort),],
                                     casesTF=in.d.cases.FI[which(in.d$STUDY_ANCESTRY==cohort)],
                                     femaleTF=in.d.females[which(in.d$STUDY_ANCESTRY==cohort)])
  
  ###for transform
  tmp.FI <- transform(tmp.FI, logFastingInsulin_adj_invn = inversenormal(logFastingInsulin_adj_by_sex))
  tmp.FI <- transform(tmp.FI, logFastingInsulin_adj_invn_males = inversenormal(logFastingInsulin_adj_males))
  tmp.FI <- transform(tmp.FI, logFastingInsulin_adj_invn_females= inversenormal(logFastingInsulin_adj_females))
  
  pooled.traits$logFastingInsulin_adj[which(in.d$STUDY_ANCESTRY==cohort)] <- tmp.FI[,"logFastingInsulin_adj_by_sex"]
  pooled.traits$logFastingInsulin_adj_invn[which(in.d$STUDY_ANCESTRY==cohort)] <- tmp.FI[,"logFastingInsulin_adj_invn"]
  
  
  tmp.cohort.traits.FI <- cbind(logFastingInsulin_adj=rep(NA,nrow(pooled.traits)),
                                logFastingInsulin_adj_invn=rep(NA,nrow(pooled.traits)))
  
  tmp.cohort.traits.FI[which(in.d$STUDY_ANCESTRY==cohort),"logFastingInsulin_adj"] <- tmp.FI[,"logFastingInsulin_adj_by_sex"]
  tmp.cohort.traits.FI[which(in.d$STUDY_ANCESTRY==cohort),"logFastingInsulin_adj_invn"] <- tmp.FI[,"logFastingInsulin_adj_invn"]
  colnames(tmp.cohort.traits.FI) <- paste(colnames(tmp.cohort.traits.FI),"_",cohort,sep="")
  

  cohort.traits <- cbind(cohort.traits,tmp.cohort.traits.FI)
  
}

for(cohort in WHI) {
  print(paste("COHORT:",cohort))
  
  in.d[which(in.d$STUDY_ANCESTRY==cohort),]
  
  print(paste("FASTING INSULIN"))
  tmp.FI <- make.residuals.transform("logFastingInsulin","age_FI+BMI_FI", controls.only=T, by.sex=F,
                                     in.d[which(in.d$STUDY_ANCESTRY==cohort),],
                                     casesTF=in.d.cases.FI[which(in.d$STUDY_ANCESTRY==cohort)],
                                     femaleTF=in.d.females[which(in.d$STUDY_ANCESTRY==cohort)])
  
  ##for transform
  tmp.FI <- transform(tmp.FI, logFastingInsulin_adj_invn = inversenormal(logFastingInsulin_adj))

  
  pooled.traits$logFastingInsulin_adj[which(in.d$STUDY_ANCESTRY==cohort)] <- tmp.FI[,"logFastingInsulin_adj"]
  pooled.traits$logFastingInsulin_adj_invn[which(in.d$STUDY_ANCESTRY==cohort)] <- tmp.FI[,"logFastingInsulin_adj_invn"]
  
  
  tmp.cohort.traits.FI <- cbind(logFastingInsulin_adj=rep(NA,nrow(pooled.traits)),
                                logFastingInsulin_adj_invn=rep(NA,nrow(pooled.traits)))
  
  tmp.cohort.traits.FI[which(in.d$STUDY_ANCESTRY==cohort),"logFastingInsulin_adj"] <- tmp.FI[,"logFastingInsulin_adj"]
  tmp.cohort.traits.FI[which(in.d$STUDY_ANCESTRY==cohort),"logFastingInsulin_adj_invn"] <- tmp.FI[,"logFastingInsulin_adj_invn"]
  colnames(tmp.cohort.traits.FI) <- paste(colnames(tmp.cohort.traits.FI),"_",cohort,sep="")
  

  cohort.traits <- cbind(cohort.traits,tmp.cohort.traits.FI)
  
}
sink()