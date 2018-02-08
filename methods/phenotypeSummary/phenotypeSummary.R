input_args <- commandArgs(trailingOnly=T)
ped.file <- input_args[1]
outcome <- input_args[2]
covars <- unlist(strsplit(input_args[3],","))
conditional.string <- input_args[4]
ivar.string <- input_args[5]
group.var <- input_args[6]
label <- input_args[7]
cohort_column <- input_args[8]
sex_column <- input_args[9]

# add ivars
if (!(ivar.string == "NA")){
  covars <- c(covars,unlist(strsplit(ivar.string,",")))
}

# add het var
if (!(group.var == "NA")){
  covars <- c(covars,group.var)
}

# get the length of covars
ncovar <- length(covars)

# If this is conditional, combine with covariates
if (!(conditional.string == "NA")){
  conditional = unlist(strsplit(conditional.string,","))
  covars = c(covars,conditional)
}

# Load phenotype data
library(data.table)
ped.data <- fread(ped.file,header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
print(head(ped.data))
print(summary(ped.data))

print(cohort_column)
print(outcome)
print(covars)
ped.data = na.omit(as.data.frame(ped.data[,unique(c(cohort_column,sex_column,outcome,covars)),drop=F]))
print(summary(ped.data))

# Change phenotype names
if (!(conditional.string == "NA")){
  for (c in seq(1,length(conditional))){
     conditional[c] = sub(":","_",conditional[c])
     if(grepl("[[:digit:]]", substr(conditional[c], 1, 1))){
       conditional[c] <- paste("chr",conditional[c],sep="")
     }
     colnames(ped.data)[ncovar+1+c] <- conditional[c]
     covars[ncovar+c] <- conditional[c]
  }
}

# Determine outcome type
is.continuous <- ifelse(length(unique(ped.data[,outcome])) < 4 ,F,T)

# Get number of possible outcomes
if (!(is.continuous)){
  outcome.vals <- sort(unique(ped.data[,outcome]))
} else {
  covars = c(outcome,covars)
}

# Remove cohort column and any study based columns if its in covars
covars = covars[!(covars %in% c(cohort_column,"STUDY_ANCESTRY","study_ancestry","topmed_project"))]

# Determine types of covars
val.type <- apply(ped.data[,covars], 2, function(x) ifelse(length(unique(x)) <= 2 ,"categorical","continuous"))

# Gather each cohorts stats in a list to be combined later
stats <- list()

# Two different loops for catgorical vs continuous
if (!(is.continuous)){
  # Loop through cohorts
  for (study in unique(ped.data[,cohort_column])){
    
    # Generate the row names that we'll use
    row_names <- c(paste(study,"total",sep=" "))
    for (v in outcome.vals){
      row_names <- c(row_names,paste(study,":",outcome,"=",v,sep=" "))
    }
    
    # Subset the phenotype data down to only the study we're currently using
    ped.cur = ped.data[ped.data[,cohort_column] == study,]
    
    # Get the total number of samples in this cohort
    total <- length(ped.cur[,1])
    
    # Start collecting numbers of samples per outcome condition
    samp <- c(total)
    for (v in outcome.vals){
      samp <- c(samp, length(ped.cur[ped.cur[,outcome] == v,1]))
    }
    
    # This is where we will store all stats for this cohort
    all_dat <- data.frame(V1 = row_names, Samples=samp)
    
    # Loop through the continuous covariates to calculate mean, median, std, min, max
    for (c in covars[val.type == "continuous"]){
      if (c == cohort_column){
        next
      }
      
      # Store the values for the whole cohort first
      means <- c(mean(ped.cur[,c]))
      medians <- c(median(ped.cur[,c]))
      sds <- c(sd(ped.cur[,c]))
      mins <- c(min(ped.cur[,c]))
      maxes <- c(max(ped.cur[,c]))
      
      # Loop through each outcome condition and calculate values
      for (v in outcome.vals){
        ped.new <- ped.cur[ped.cur[,outcome] == v,c]
        means = c(means,mean(ped.new))
        medians = c(medians, median(ped.new))
        sds = c(sds, sd(ped.new))
        mins = c(mins, min(ped.new))
      }
      
      # Store the values back in the whole data frame
      all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
      
    }
    
    # Loop through the categorical covariates
    for (c in covars[val.type == "categorical"]){
      if (c == cohort_column){
        next
      }
      
      # Get all possible values of the current covar
      curvals <- unique(ped.data[,c])
      
      # Loop through the possible values of this covar
      for (cv in curvals){
        
        # Subset the phenotypes to only cases with this covar value
        pd.c <- ped.cur[ped.cur[,c] == cv,]
        
        # Total percent with this value for the covar
        percents <- c(length(pd.c[,1])/total)
        
        # Calculate percents of samples with this value over total with given condition
        for (v in outcome.vals){
          ped.new <- pd.c[pd.c[,outcome] == v,c]
          percents <- c(percents,length(ped.new)/length(ped.cur[ped.cur[,outcome] == v,1]))
        }
        
        # Store back in the total data frame
        all_dat <- cbind(all_dat,percents)  
      }
    }
    
    # Store the cohort specific stats in the list
    stats[[length(stats)+1]] <- all_dat
  }
  
  # Get totals for all cohorts
  row_names <- c("All cohorts")
  for (v in outcome.vals){
    row_names <- c(row_names,paste("All cohorts",":",outcome,"=",v,sep=" "))
  }
  
  # total samples
  total <- length(ped.data[,1])
  samp <- c(total)
  for (v in outcome.vals){
    samp <- c(samp, length(ped.data[ped.data[,outcome] == v,1]))
  }
  
  all_dat <- data.frame(V1 = row_names, Samples=samp)
  
  # Loop through the continuous covariates to calculate mean, median, std, min, max
  for (c in covars[val.type == "continuous"]){
    if (c == cohort_column){
      next
    }
    
    # Store the values for the whole cohort first
    means <- c(mean(ped.data[,c]))
    medians <- c(median(ped.data[,c]))
    sds <- c(sd(ped.data[,c]))
    mins <- c(min(ped.data[,c]))
    maxes <- c(max(ped.data[,c]))
    
    # Loop through each outcome condition and calculate values
    for (v in outcome.vals){
      ped.new <- ped.data[ped.data[,outcome] == v,c]
      means = c(means,mean(ped.new))
      medians = c(medians, median(ped.new))
      sds = c(sds, sd(ped.new))
      mins = c(mins, min(ped.new))
    }
    
    # Store the values back in the whole data frame
    all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
    
  }
  
  # Loop through the categorical covariates
  for (c in covars[val.type == "categorical"]){
    if (c == cohort_column){
      next
    }
    
    # Get all possible values of the current covar
    curvals <- unique(ped.data[,c])
    
    # Loop through the possible values of this covar
    for (cv in curvals){
      
      # Subset the phenotypes to only cases with this covar value
      pd.c <- ped.data[ped.data[,c] == cv,]
      
      # Total percent with this value for the covar
      percents <- c(length(pd.c[,1])/total)
      
      # Calculate percents of samples with this value over total with given condition
      for (v in outcome.vals){
        ped.new <- pd.c[pd.c[,outcome] == v,c]
        percents <- c(percents,length(ped.new)/length(ped.data[ped.data[,outcome] == v,1]))
      }
      
      # Store back in the total data frame
      all_dat <- cbind(all_dat,percents)  
    }
    stats[[length(stats)+1]] <- all_dat
  }
  
# If we have a continuous outcome
} else {
  for (study in unique(ped.data[,cohort_column])){
    
    # Subset the phenotype data down to only the study we're currently using
    ped.cur = ped.data[ped.data[,cohort_column] == study,]
    
    # Get the total number of samples in this cohort
    total <- length(ped.cur[,1])
    
    # This is where we will store all stats for this cohort
    all_dat <- data.frame(V1 = study, Samples=total)
    
    # Loop through the continuous covariates to calculate mean, median, std, min, max
    for (c in covars[val.type == "continuous"]){
      if (c == cohort_column){
        next
      }
      
      # Store the values for the whole cohort first
      means <- c(mean(ped.cur[,c]))
      medians <- c(median(ped.cur[,c]))
      sds <- c(sd(ped.cur[,c]))
      mins <- c(min(ped.cur[,c]))
      maxes <- c(max(ped.cur[,c]))
      
      # Store the values back in the whole data frame
      all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
      
    }
    
    # Loop through the categorical covariates
    for (c in covars[val.type == "categorical"]){
      if (c == cohort_column){
        next
      }
      
      # Get all possible values of the current covar
      curvals <- unique(ped.data[,c])
      
      # Loop through the possible values of this covar
      for (cv in curvals){
        
        # Subset the phenotypes to only cases with this covar value
        pd.c <- ped.cur[ped.cur[,c] == cv,]
        
        # Total percent with this value for the covar
        percents <- c(length(pd.c[,1])/total)
        
        # Store back in the total data frame
        all_dat <- cbind(all_dat,percents)  
      }
    }
    
    # Store the cohort specific stats in the list
    stats[[length(stats)+1]] <- all_dat
  }
  
  # Get the total number of samples in this cohort
  total <- length(ped.data[,1])
  
  # This is where we will store all stats for this cohort
  all_dat <- data.frame(V1 = "All cohorts", Samples=total)
  
  # Loop through the continuous covariates to calculate mean, median, std, min, max
  for (c in covars[val.type == "continuous"]){
    if (c == cohort_column){
      next
    }
    
    # Store the values for the whole cohort first
    means <- c(mean(ped.data[,c]))
    medians <- c(median(ped.data[,c]))
    sds <- c(sd(ped.data[,c]))
    mins <- c(min(ped.data[,c]))
    maxes <- c(max(ped.data[,c]))
    
    # Store the values back in the whole data frame
    all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
    
  }
  
  # Loop through the categorical covariates
  for (c in covars[val.type == "categorical"]){
    if (c == cohort_column){
      next
    }
    
    # Get all possible values of the current covar
    curvals <- unique(ped.data[,c])
    
    # Loop through the possible values of this covar
    for (cv in curvals){
      
      # Subset the phenotypes to only cases with this covar value
      pd.c <- ped.data[ped.data[,c] == cv,]
      
      # Total percent with this value for the covar
      percents <- c(length(pd.c[,1])/total)
      
      # Store back in the total data frame
      all_dat <- cbind(all_dat,percents)  
    }
  }
  
  # Store the cohort specific stats in the list
  stats[[length(stats)+1]] <- all_dat
}
# Determine the column and row names for the data frame
# Continuous covars go first
cont_vals <- covars[val.type == "continuous"]
cont_cols <- c()
for (c in cont_vals){
  if (c == cohort_column){
    next
  }
  cont_cols <- c(cont_cols, paste(c,"mean",sep=" "), paste(c,"median",sep=" "), paste(c,"STD",sep=" "), paste(c,"min",sep=" "), paste(c,"max",sep=" "))
}

# Next are categorical covars
dich_vals <- covars[val.type == "categorical"]
dich_cols <- c()
for (c in dich_vals){
  if (c == outcome){
    next
  }
  uvals <- unique(ped.cur[,c])
  dich_cols <- c(dich_cols, paste(c,"%",uvals[1],sep=" "),paste(c,"%",uvals[2],sep=" "))
}

# Add samples to headers and combine
col_names <- c("Samples",cont_cols,dich_cols)

# Combine list to data frame containing all cohorts
stats.df <- as.data.frame(do.call(rbind,stats))#, row.names = row_names)

# Move the row names to row.names of data frame
row.names(stats.df) <- stats.df$V1

# Remove row names column
stats.df <- stats.df[,2:length(stats.df[1,])]

# Assign the correct column names
colnames(stats.df) <- col_names

# Write the final stats out to a csv file
fwrite(stats.df,file=paste(label,"_stats.csv",sep=""),sep=",",row.names = T)

# Make some plots
# Load ggplot for plots
library(ggplot2)

# Gather continuous covariates
quant_covars <- covars[val.type == "continuous"]

# Make sure we dont include the cohort column
quant_covars = quant_covars[quant_covars != cohort_column]

# Initialize the pdf
pdf(paste(label,"_plots.pdf",sep=""),width=11)

# For each continuous covariate
for (i in quant_covars){
  print(i)
  
  # Plot all samples
  plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
  print(plot + geom_boxplot() + labs(title = "All samples",x="Cohort",y=i))
  
  if (!(is.continuous)){
  # Plot samples devided by outcome    
  plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1))
  print(plot + geom_boxplot(aes(fill=factor(get(outcome)))) + labs(title = "All samples by outcome",x="Cohort",y=i,fill=aes_string(outcome)))
  }
  # Plot samples devided by sex
  plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
  print(plot + geom_boxplot(aes(fill=factor(get(sex_column)))) + labs(title = "All samples by sex",x="Cohort",y=i,fill=sex_column))
}

dev.off()

