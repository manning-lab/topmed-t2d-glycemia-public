# This script takes the results from a GENESIS
# association test and output (1) the results
# in a CSV (2) the top results in a csv (3)
# a manhattan plot and (4) a QQ plot of the 
# p-values

## Example usage: R --vanilla --args GoT2D_T2D.assoc.RData Score.pval label T2D_covariates red black F F < summary.R

#install.packages("calibrate_1.7.2.tar.gz")
install.packages("qqman",repos='http://cran.us.r-project.org')  #_0.14.tar.gz")
library(qqman)
library(data.table)

input_args <- commandArgs(trailingOnly=T)
assoc.file <- input_args[1]
label <- input_args[2]

load(assoc.file)

ppi <- 300
results <- assoc$results
results$chr <- rep(10,length(results[,1]))
pdf(paste(label,".qqplot.pdf",sep="")
qq(results$Score.pval)
dev.off()

l <- list()
for (i in seq(1,length(row.names(results)))){
  l[[length(l)+1]] <- data.frame(P=rep(results$Score.pval[i],length(assoc$variantInfo[[i]][,1])), BP=assoc$variantInfo[[i]]$pos, CHR=assoc$variantInfo[[i]]$chr)
}

df <- l[[1]]
for (i in seq(2,length(l))){
  df <- rbind(df,l[[i]])
}
df$CHR <- as.numeric(as.vector(df$CHR))
pdf(paste(label,".mhplot.pdf",sep="")
manhattan(df,chr="CHR",bp="BP",p="P", main=label)
dev.off()

write.csv(results, paste(label, ".groupAssoc.csv", sep=""))

