# This script takes the results from a GENESIS
# association test and output (1) the results
# in a CSV (2) the top results in a csv (3)
# a manhattan plot and (4) a QQ plot of the 
# p-values

## Example usage: R --vanilla --args GoT2D_T2D.assoc.RData Score.pval label T2D_covariates red black F F < summary.R

#install.packages("calibrate_1.7.2.tar.gz")
install.packages("qqman",repos='http://cran.us.r-project.org')  #_0.14.tar.gz")
library(qqman)

input_args <- commandArgs(trailingOnly=T)

pval <- input_args[1] # name of P-value column
label <- input_args[2] 
title <- input_args[3] # plot label (Q: How do I have user input separated by periods and then split for label?)

## load files and add to assoc.compilation

assoc.files <- c() # list of input .assoc.RData files
assoc.compilation <- c() # matrix of association results
numAssocFiles <- (length(input_args) - 3)
print(numAssocFiles)
print(input_args)

print(date())

for (i in 1:numAssocFiles) {
	print(i)
	assoc.files[i] <- input_args[i+3]

	print(assoc.files[i])
	load(assoc.files[i])
#	print(summary(assoc))
	#fix column names and types
#	names(assoc)[names(assoc)=="chr"] <- "CHR"
#	names(assoc)[names(assoc)=="rsid"] <- "SNP"
#	names(assoc)[names(assoc)=="pos"] <- "BP"
#	names(assoc)[names(assoc)==pval] <- "P"
#	assoc <-  transform(assoc, CHR = as.numeric(CHR))
#	assoc <-  transform(assoc, P = as.numeric(P))
#	assoc <-  transform(assoc, BP = as.numeric(BP))
	
	#remove rows with null values
	print(dim(assoc))
	assoc <- assoc[!is.na(assoc[,pval]),]
	print(dim(assoc))

	#add to assoc.compilation
#	assoc.compilation <- rbind(assoc.compilation, assoc)

	if (i == 1) {
		write.table(assoc,paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
	} else {
		write.table(assoc,paste(label, ".assoc.csv", sep=""),col.names=FALSE,sep=",",row.names=F, append=TRUE)
	}
		
}
print(date())
print(list.files())
install.packages("data.table",repos='http://cran.us.r-project.org')
library(data.table)
assoc.compilation <- fread(paste(label, ".assoc.csv", sep=""),sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)

print(date())

	print(dim(assoc.compilation))

print(summary(assoc.compilation))

assoc.compilation$chr <- as.numeric(as.character(assoc.compilation$chr))
assoc.compilation$pos <- as.numeric(as.character(assoc.compilation$pos))
assoc.compilation$P <- as.numeric(as.character(assoc.compilation[,pval]))

print(summary(assoc.compilation))

# save association results as CSV
#write.csv(assoc.compilation, paste(label, ".assoc.csv", sep=""))
write.csv(assoc.compilation[assoc.compilation[,pval] < 0.0001, ], paste(label, ".topassoc.csv", sep=""))

ppi <- 300

# plot results: Q-Q
png(paste(label, ".qqplot.png", sep=""),width = 7*ppi, height = 7*ppi, res=ppi, type="cairo")
qq(assoc.compilation$P)
dev.off()

# plot results
png(paste(label, ".mhplot.png", sep=""), width = 7*ppi, height = 7*ppi, res = ppi, type="cairo")
manhattan(assoc.compilation,chr="chr",bp="pos",p="P", main=title)
dev.off()


