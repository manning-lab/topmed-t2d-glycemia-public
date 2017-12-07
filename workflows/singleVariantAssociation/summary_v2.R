install.packages("qqman",repos='http://cran.us.r-project.org')  #_0.14.tar.gz")
library(qqman)

input_args <- commandArgs(trailingOnly=T)

pval <- input_args[1]
label <- input_args[2] 

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

  if (!is.na(assoc)[1]){
    print(dim(assoc))
    assoc <- assoc[!is.na(assoc[,pval]),]
    print(dim(assoc))
    
    if (i == 1) {
      write.table(assoc,paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
    } else {
      write.table(assoc,paste(label, ".assoc.csv", sep=""),col.names=FALSE,sep=",",row.names=F, append=TRUE)
    }	
  }
}

library(data.table)
assoc.compilation <- fread(paste(label, ".assoc.csv", sep=""),sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)

assoc.compilation$chr <- as.numeric(as.character(assoc.compilation$chr))
assoc.compilation$pos <- as.numeric(as.character(assoc.compilation$pos))
assoc.compilation$P <- as.numeric(as.character(assoc.compilation[,pval]))

write.csv(assoc.compilation[assoc.compilation[,pval] < 0.0001, ], paste(label, ".topassoc.csv", sep=""))

ppi <- 300

# plot results: Q-Q
png(paste(label, ".all.qqplot.png", sep=""),width = 7*ppi, height = 7*ppi, res=ppi, type="cairo")
qq(assoc.compilation$P,main="All variants")
dev.off()

png(paste(label, ".common.qqplot.png", sep=""),width = 7*ppi, height = 7*ppi, res=ppi, type="cairo")
qq(assoc.compilation$P[assoc.compilation$MAF>=0.05],main="Variants with MAF>=0.05")
dev.off()

png(paste(label, ".uncommon.qqplot.png", sep=""),width = 7*ppi, height = 7*ppi, res=ppi, type="cairo")
qq(assoc.compilation$P[assoc.compilation$MAF<=0.05],main="Variants with MAF<=0.05")
dev.off()

# plot results
png(paste(label, ".all.mhplot.png", sep=""), width = 7*ppi, height = 7*ppi, res = ppi, type="cairo")
manhattan(assoc.compilation,chr="chr",bp="pos",p="P", main="All variants")
dev.off()

png(paste(label, ".common.mhplot.png", sep=""), width = 7*ppi, height = 7*ppi, res = ppi, type="cairo")
manhattan(assoc.compilation[assoc.compilation$MAF>=0.05,],chr="chr",bp="pos",p="P", main="Variants with MAF>=0.05")
dev.off()

png(paste(label, ".uncommon.mhplot.png", sep=""), width = 7*ppi, height = 7*ppi, res = ppi, type="cairo")
manhattan(assoc.compilation[assoc.compilation$MAF<=0.05,],chr="chr",bp="pos",p="P", main="Variants with MAF<=0.05")
dev.off()
