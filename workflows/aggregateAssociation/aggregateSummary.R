# This script takes the results from a GENESIS
# association test and output (1) the results
# in a CSV (2) the top results in a csv (3)
# a manhattan plot and (4) a QQ plot of the 
# p-values

## Example usage: R --vanilla --args GoT2D_T2D.assoc.RData Score.pval label T2D_covariates red black F F < summary.R

#install.packages("calibrate_1.7.2.tar.gz")
# install.packages("qqman",repos='http://cran.us.r-project.org')  #_0.14.tar.gz")
library(qqman)
library(data.table)

input_args <- commandArgs(trailingOnly=T)
label <- input_args[1]

assoc.files <- c() # list of input .assoc.RData files
assoc.compilation <- c() # matrix of association results
numAssocFiles <- (length(input_args) - 1)
all_assoc <- list()

for (i in 1:numAssocFiles) {
  print(i)
  assoc.files[i] <- input_args[i+1]
  load(assoc.files[i])
  res <- assoc$results
  all_assoc[[i]]<-assoc

  if (!is.na(res)[1]){
    print(dim(res))
    res <- res[!is.na(res[,"pval_0"]),]
    
    #add to assoc.compilation
    assoc.compilation <- rbind(assoc.compilation, res)
    
    if (i == 1) {
      write.table(res,paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
    } else {
      write.table(res,paste(label, ".assoc.csv", sep=""),col.names=FALSE,sep=",",row.names=F, append=TRUE)
    }	
  }
}

# load(assoc.file)

ppi <- 300
# results <- assoc$results
# results$chr <- rep(10,length(results[,1]))
pdf(paste(label,".qqplot.pdf",sep=""))
qq(res$pval_0)
dev.off()

l <- list()
for (i in seq(1,length(all_assoc))){
# for (i in seq(1,length(res[,1]))){
  l2 <- list()
    for (j in seq(1,length(all_assoc[[i]]$variantInfo))){
      l2[[length(l2)+1]] <- data.frame(P=rep(all_assoc[[i]]$results$pval_0[j],length(all_assoc[[i]]$variantInfo[[j]][,1])), BP=all_assoc[[i]]$variantInfo[[j]]$pos, CHR=all_assoc[[i]]$variantInfo[[j]]$chr)
  # l[[length(l)+1]] <- data.frame(P=rep(res$pval_0[i],length(groups[[i]]$variant.id)), BP=groups[[i]]$position, CHR=groups[[i]]$chromosome)
    }
  l <- unlist(list(l,l2),recursive=F)
}

df <- l[[1]]
for (i in seq(2,length(l))){
  df <- rbind(df,l[[i]])
}
df$CHR <- as.numeric(as.vector(df$CHR))
pdf(paste(label,".mhplot.pdf",sep=""))
manhattan(df,chr="CHR",bp="BP",p="P", main=label)
dev.off()

write.csv(results, paste(label, ".groupAssoc.csv", sep=""))

