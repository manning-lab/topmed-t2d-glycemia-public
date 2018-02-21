require(ggplot2)

compare.residuals <- function(trait,CCnotsep,CCsep,CCsep.transform,CCsep.transform.adjcohort,cohortCaseControl,givenGender,countryCaseControl,rawtrait,rawtrait.label) {
  
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y)
  }
  plot(1:5,1:5,type="n",axes=FALSE,xlab="",ylab="")
  legend("bottomleft",pch=19,col=rainbow(length(levels(as.factor(cohortCaseControl)))),legend=levels(as.factor(cohortCaseControl)))
  
  pairs(data.frame(CCnotsep=CCnotsep[,trait],
                   CCsep=CCsep[,trait],
                   CCsep.transform=CCsep.transform[,trait],
                   CCsep.transform.adjcohort=CCsep.transform.adjcohort[,trait]),
        col=rainbow(13)[as.factor(phenotypes$cohortCaseControl)],main=trait,panel=panel.smooth,diag.panel=panel.hist)
   
  
  make.plots.cohort2(rawtrait,givenGender,cohortCaseControl,mainlabel=rawtrait.label)
  make.plots.cohort2(CCnotsep[,trait],givenGender,countryCaseControl,mainlabel=paste(trait,"CCnotsep"))
  make.plots.cohort2(CCsep[,trait],givenGender,cohortCaseControl,mainlabel=paste(trait,"CCsep"))
  make.plots.cohort2(CCsep[,trait],givenGender,countryCaseControl,mainlabel=paste(trait,"CCsep"))
  make.plots.cohort2(CCsep.transform[,trait],givenGender,cohortCaseControl,mainlabel=paste(trait,"CCsep.transform"))
  make.plots.cohort2(CCsep.transform[,trait],givenGender,countryCaseControl,mainlabel=paste(trait,"CCsep.transform"))
  make.plots.cohort2(CCsep.transform.adjcohort[,trait],givenGender,cohortCaseControl,mainlabel=paste(trait,"CCsep.transform.adjcohort"))
  make.plots.cohort2(CCsep.transform.adjcohort[,trait],givenGender,countryCaseControl,mainlabel=paste(trait,"CCsep.transform.adjcohort"))
}

make.plots.cohort2 <-  function(trait,sexinfo="givenGender",grouping,mainlabel="") {

  qplotStatement <- qplot(trait, facets=grouping~., fill=grouping, colour=sexinfo, 
                            data=data.frame(trait,sexinfo,grouping),main=mainlabel) + theme_bw() 
                            
  #print(qplotStatement)
  print(qplotStatement)
}


make.plots <-  function(columnName,phenotypes=phenotypes,sexinfo="givenGender") {
  print("In make plots")
  qplotStatement <- paste(
    "qplot(",
    columnName,
    ", facets=countryCaseControl~., fill=countryCaseControl, colour=",sexinfo,", data=phenotypes[ !is.na(phenotypes[,",
    columnName,
    "]),]) + theme_bw() + opts(axis.text.x= theme_text(angle= 90,hjust=1))",
    sep=""
  )
  #print(qplotStatement)
  print(eval(parse(text=qplotStatement)))
}


inversenormal <- function(x) {
  # inverse normal if you have missing data
  return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
}

ztrans <- function(x) {
  mu = mean(x, na.rm=TRUE)
  sigma = sd(x, na.rm=TRUE)
  return((x - mu)/sigma)
}


do.lm <- function(formula.2, pheno.2, subset.2) {
  print(dim(pheno.2[which(subset.2),]))
  print(summary(lm(formula.2,data=pheno.2[which(subset.2),],na.action="na.exclude")))
  residuals(lm(formula.2,data=pheno.2[which(subset.2),],na.action="na.exclude"))
}

make.residuals.combined <- function(outcome,covariates,pheno,by.sex=F,transform=identity,controls.only=F,label="_adj",covariates.cases=NA,return.raw=T,casesTF=NA,femaleTF=NA) {
  if(by.sex==F) {
    residuals <- rep(NA,nrow(pheno))
  
    if(controls.only) {
      residuals[!casesTF] <- transform(do.lm(formula.2=as.formula(paste(outcome,covariates,sep="~")),
                                   pheno.2=pheno,subset.2=!casesTF))
      
    } else {
      residuals <- transform(do.lm(formula.2=as.formula(paste(outcome,covariates,sep="~")),
                                  pheno.2=pheno,subset.2=rep(TRUE,nrow(pheno))))
    }
    
    residual.name <- paste(outcome,label,sep="")

    
    if(return.raw) {
      to.return <- cbind(pheno[,outcome],residuals)
      colnames(to.return) <- c(outcome,residual.name)
    } else {
      to.return <- cbind(residuals)
      colnames(to.return) <- residual.name 
    }    
    return(to.return)
    
  } else if(by.sex==T) {
    ## sex-specific residuals and transformations  
    residuals.males <- rep(NA,nrow(pheno))
    residuals.females <- rep(NA,nrow(pheno))
    
    print("controls,males")
    residuals.males[ !femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2= !femaleTF))
    
    print("controls,females")
    residuals.females[femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=femaleTF))
        
    to.return <- cbind(residuals.males,residuals.females)
    colnames(to.return) <- c(paste(outcome,label,c("_males","_females"),sep=""))
    return(to.return)
    
  }
}
  
make.residuals <- function(outcome,covariates,pheno,by.sex=F,transform=identity,controls.only=F,label="_adj",covariates.cases=NA,return.raw=T,casesTF=NA,femaleTF=NA) {
 
  if(is.na(covariates.cases)) covariates.cases=covariates
  
  residuals <- rep(NA,nrow(pheno))
  
  if(by.sex==F) {
    ## case specific
    residuals[which(!casesTF)] <- transform(do.lm(formula.2=as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF))
    
    if(controls.only==F) {
      residuals[which(casesTF)] <- transform(do.lm(formula.2=as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=casesTF))
    }  
  
    residual.name <- paste(outcome,label,sep="")
    
    
    if(return.raw) {
      to.return <- cbind(pheno[,outcome],residuals)
      
      colnames(to.return) <- c(outcome,residual.name)
      
    } else {
      to.return <- cbind(residuals)
      colnames(to.return) <- residual.name 
    }  
    
    
    return(to.return)
    
    
    
  } else {
    ## sex-specific residuals and transformations  
    residuals.males <- rep(NA,nrow(pheno))
    residuals.females <- rep(NA,nrow(pheno))

    
    print("controls,males")
    residuals[!casesTF & !femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & !femaleTF))
    residuals.males[!casesTF & !femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & !femaleTF))

    print("controls,females")
    residuals[!casesTF & femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & femaleTF))
    residuals.females[!casesTF & femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=!casesTF & femaleTF))
    
    if(controls.only==F) {
      print("cases,males")
      residuals[casesTF & !femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates.cases,sep="~")),pheno.2=pheno,subset.2=casesTF & !femaleTF))
      residuals.males[casesTF & !femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=casesTF & !femaleTF))
    
      print("cases,females")
      residuals[casesTF & femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates.cases,sep="~")),pheno.2=pheno,subset.2=casesTF & femaleTF))
      residuals.females[casesTF & femaleTF] <- transform(do.lm(as.formula(paste(outcome,covariates,sep="~")),pheno.2=pheno,subset.2=casesTF & femaleTF))
    }
    
    residual.name <- paste(outcome,label,"_by_sex",sep="")

    to.return <- cbind(residuals,residuals.males,residuals.females)
    colnames(to.return) <- c(residual.name,paste(outcome,label,c("_males","_females"),sep=""))
    return(to.return)
    }

}

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
    residuals <- transform(residuals)
    
    if(return.raw) {
      to.return <- cbind(pheno[,outcome],residuals)
      
      colnames(to.return) <- c(outcome,residual.name)
      
    } else {
      to.return <- cbind(residuals)
      colnames(to.return) <- residual.name 
    }  
    
    
    return(to.return)
    
    
    
  } else {
    ## sex-specific residuals and transformations  
    residuals.males <- rep(NA,nrow(pheno))
    residuals.females <- rep(NA,nrow(pheno))
    
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
    
    residual.name <- paste(outcome,label,"_by_sex",sep="")
    
    ## apply transformation
    residuals <- transform(residuals)
    residuals.females <- transform(residuals.females)
    residuals.males <- transform(residuals.males)
    
    to.return <- cbind(residuals,residuals.males,residuals.females)
    colnames(to.return) <- c(residual.name,paste(outcome,label,c("_males","_females"),sep=""))
    return(to.return)
  }
  
}



library(e1071)
get.stats <- function(colname) {
  col <- phenotypes.transformed[,colname]
  f.names <- c("mean","sd","min","max","skewness","kurtosis")
  to.return <- c()
  if(is.numeric(col)) {
    to.return <- c(column=colname,N=format(sum(!is.na(col))),
                   sapply(f.names,function(f.name){f<-get(f.name);return(f.name=f(col,na.rm=T))}))
  } else {
    to.return <- c(column=colname,N=format(sum(!is.na(col))),rep(NA,6))
  } 
  names(to.return) <- c("colname","N.nonmissing",f.names)
  return(to.return)
}