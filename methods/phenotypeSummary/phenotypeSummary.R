input_args <- commandArgs(trailingOnly=T)
ped.file <- input_args[1]
outcome <- input_args[2]
covars <- unlist(strsplit(input_args[3],","))
label <- input_args[4]
cohort_column <- input_args[5]

# ped.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_AFEU_WesselJ_25AUG2017_T2D.csv"
# outcome <- "t2d_ctrl"
# covars <- unlist(strsplit("last_exam_age,sex,study",","))
# label <- "testing"
# cohort_column <- "study"

# Load ped.data data
library(data.table)
ped.data <- fread(ped.file,header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
ped.data = ped.data[,c(outcome,covars)]

# determine outcome type
is.continuous <- ifelse(length(unique(ped.data[,outcome])) < 4 ,F,T)

# get number of outcomes
if (!(is.continuous)){
  outcome.vals <- sort(unique(ped.data[,outcome]))
}

# remove cohort column and any study based columns if its in covars
covars = covars[!(covars %in% c(cohort_column,"STUDY_ANCESTRY","study_ancestry","topmed_project"))]



# determine types of covars
val.type <- apply(ped.data[,covars], 2, function(x) ifelse(length(unique(x)) <= 2 ,"categorical","continuous"))

stats <- list()
# row_names <- c()
for (study in unique(ped.data[,cohort_column])){
  row_names <- c(paste(study,"total",sep=" "))
  for (v in outcome.vals){
    row_names <- c(row_names,paste(study,":",outcome,"=",v,sep=" "))
  }
  # row_names <- c(row_names,paste(study,"total",sep=" "),paste(study,"controls",sep=" "),paste(study,"cases",sep=" "))
  ped.cur = ped.data[ped.data[,cohort_column] == study,]
  total <- length(ped.cur[,1])
  samp <- c(total)
  for (v in outcome.vals){
    samp <- c(samp, length(ped.cur[ped.cur[,outcome] == v,1]))
  }
  
  
  all_dat <- data.frame(V1 = row_names, Samples=samp)
  # total_dat <- c(length(ped.cur[,1]))
  # case_dat <- c(length(ped.cur[ped.cur[,outcome] == 1,1]))
  # ctrl_dat <- c(length(ped.cur[ped.cur[,outcome] == 0,1]))
  for (c in covars[val.type == "continuous"]){
  #   if (c == cohort_column){
  #     next
  #   }
    # total_dat <- c(total_dat, mean(ped.cur[,c], na.rm = T), median(ped.cur[,c], na.rm = T), sd(ped.cur[,c], na.rm = T), min(ped.cur[,c], na.rm = T), max(ped.cur[,c], na.rm = T))
    # case_dat <- c(case_dat, mean(ped.cur[ped.cur[,outcome] == 1,c], na.rm = T), median(ped.cur[ped.cur[,outcome] == 1,c], na.rm = T), sd(ped.cur[ped.cur[,outcome] == 1,c], na.rm = T), min(ped.cur[ped.cur[,outcome] == 1,c], na.rm = T), max(ped.cur[ped.cur[,outcome] == 1,c], na.rm = T))
    # ctrl_dat <- c(ctrl_dat, mean(ped.cur[ped.cur[,outcome] == 0,c], na.rm = T), median(ped.cur[ped.cur[,outcome] == 0,c], na.rm = T), sd(ped.cur[ped.cur[,outcome] == 0,c], na.rm = T), min(ped.cur[ped.cur[,outcome] == 0,c], na.rm = T), max(ped.cur[ped.cur[,outcome] == 0,c], na.rm = T))
    means <- c(mean(ped.cur[,c]))
    medians <- c(median(ped.cur[,c]))
    sds <- c(sd(ped.cur[,c]))
    mins <- c(min(ped.cur[,c]))
    maxes <- c(max(ped.cur[,c]))
    for (v in outcome.vals){
      ped.new <- ped.cur[ped.cur[,outcome] == v,c]
      means = c(means,mean(ped.new))
      medians = c(medians, median(ped.new))
      sds = c(sds, sd(ped.new))
      mins = c(mins, min(ped.new))
    }
    
    all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
    
  }
  
  
  ############## need to do below this point
  for (c in covars[val.type == "categorical"]){
    if (c == outcome){
      next
    }
    
    
    curvals <- unique(ped.cur[,c])
    for (cv in curvals){
      
      pd.c <- ped.cur[ped.cur[,c] == cv,]
      percents <- c(length(pd.c[,1])/total)
      for (v in outcome.vals){
        ped.new <- pd.c[pd.c[,outcome] == v,c]
        percents <- c(percents,length(ped.new)/length(ped.cur[ped.cur[,outcome] == v,1]))
      }
      all_dat <- cbind(all_dat,percents)  
    }
    
    
    
    # uvals <- unique(ped.cur[,c])
    # total_dat <- c(total_dat, sum(ped.cur[,c] == uvals[1])/length(ped.cur[,1]),sum(ped.cur[,c] == uvals[2])/length(ped.cur[,1]))
    # case_dat <- c(case_dat, sum(ped.cur[ped.cur[,outcome] == 1,c] == uvals[1])/length(ped.cur[ped.cur[,outcome] == 1,c]),sum(ped.cur[ped.cur[,outcome] == 1,c] == uvals[2])/length(ped.cur[ped.cur[,outcome] == 1,c]))
    # ctrl_dat <- c(ctrl_dat, sum(ped.cur[ped.cur[,outcome] == 0,c] == uvals[1])/length(ped.cur[ped.cur[,outcome] == 0,c]),sum(ped.cur[ped.cur[,outcome] == 0,c] == uvals[2])/length(ped.cur[ped.cur[,outcome] == 0,c]))
  }
  
  # all_dat <- cbind(all_dat,percents)
  stats[[length(stats)+1]] <- all_dat
}

cont_vals <- covars[val.type == "continuous"]
cont_cols <- c()
for (c in cont_vals){
  if (c == cohort_column){
    next
  }
  cont_cols <- c(cont_cols, paste(c,"mean",sep=" "), paste(c,"median",sep=" "), paste(c,"STD",sep=" "), paste(c,"min",sep=" "), paste(c,"max",sep=" "))
}
dich_vals <- covars[val.type == "categorical"]
dich_cols <- c()
for (c in dich_vals){
  if (c == outcome){
    next
  }
  uvals <- unique(ped.cur[,c])
  dich_cols <- c(dich_cols, paste(c,"%",uvals[1],sep=" "),paste(c,"%",uvals[2],sep=" "))
}

col_names <- c("Samples",cont_cols,dich_cols)
stats.df <- as.data.frame(do.call(rbind,stats))#, row.names = row_names)
row.names(stats.df) <- stats.df$V1
stats.df <- stats.df[,2:length(stats.df[1,])]
colnames(stats.df) <- col_names

fwrite(stats.df,file=paste(label,"_stats.csv",sep=""),sep=",")

if (!(is.continuous)){
  library(ggplot2)
  
  quant_covars <- covars[val.type == "continuous"]
  quant_covars = quant_covars[quant_covars != cohort_column]
  
  ped.data <- ped.data[!is.na(ped.data[,outcome]),]
  case = ped.data[ped.data[,outcome] == 2 | ped.data[,outcome] == 1,]
  control = ped.data[ped.data[,outcome] == 0,]
  
  male.val <- "M"
  female.val <- "F"
  
  m = ped.data[ped.data$sex==male.val,]
  f = ped.data[ped.data$sex==female.val,]
  
  if (NCOL(m)==0 && NCOL(f)==0){
    male.val <- 1
    female.val <- 2
  }
  
  m_case = case[case$sex==male.val,]
  m_control = control[control$sex==male.val,]
  f_case = case[case$sex==female.val,]
  f_control = control[control$sex==female.val,]
  
  
  
  pdf(paste(label,"_plots.pdf",sep=""),width=11)
  if (length(quant_covars) > 1){
    layout(matrix(seq(1,6*(length(quant_covars)-1)),nrow=length(quant_covars)-1,ncol=6,byrow=T))
  }else{
    layout(matrix(seq(1,6*(length(quant_covars))),nrow=length(quant_covars),ncol=6,byrow=T))
  }
  for (i in quant_covars){
    print(i)
    # 
    # plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    # print(plot + geom_violin() + labs(title = "All samples"))
    # 
    # 
    # plot <- ggplot(case, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    # print(plot + geom_violin() + labs(title = "Case samples"))
    # 
    # 
    # plot <- ggplot(control, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    # print(plot + geom_violin() + labs(title = "Control samples"))
    # 
    # plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    # print(plot + geom_violin(aes(fill=factor(get(outcome)))) + labs(title = "All samples"))
    # 
    # 
    # plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    # print(plot + geom_violin(aes(fill=factor(sex))) + labs(title = "All samples by sex"))
    # 
    # 
    # plot <- ggplot(case, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    # print(plot + geom_violin(aes(fill=factor(sex))) + labs(title = "Case samples by sex"))
    # 
    # 
    # plot <- ggplot(control, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    # print(plot + geom_violin(aes(fill=factor(sex))) + labs(title = "Control samples by sex"))
    
    
    
    plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_boxplot() + labs(title = "All samples"))
    
    
    plot <- ggplot(case, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_boxplot() + labs(title = "Case samples"))
    
    
    plot <- ggplot(control, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_boxplot() + labs(title = "Control samples"))
    
    plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_boxplot(aes(fill=factor(get(outcome)))) + labs(title = "All samples"))
    
    
    plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_boxplot(aes(fill=factor(sex))) + labs(title = "All samples by sex"))
    
    
    plot <- ggplot(case, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_boxplot(aes(fill=factor(sex))) + labs(title = "Case samples by sex"))
    
    
    plot <- ggplot(control, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
    print(plot + geom_boxplot(aes(fill=factor(sex))) + labs(title = "Control samples by sex"))
  }
  
  dev.off()
} else {
  pdf(paste(label,"_plots.pdf",sep=""),width=11)
  fwrite(data.frame(na="NA"), file=paste(label,"_stats.csv",sep=""))
  dev.off()
}
