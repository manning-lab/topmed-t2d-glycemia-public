# metal summary
# Check if required packages are installed (sourced from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them)
packages <- c("qqman","data.table")
to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install,repos='http://cran.us.r-project.org')

# Load packages
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
marker.column <- input_args[1]
freq.column <- input_args[2]
pval.column <- input_args[3]
cols.tokeep <- unlist(strsplit(input_args[4],","))
assoc.names <- unlist(strsplit(input_args[5],","))
out.pref <- input_args[6]
metal.file <- input_args[7]
assoc.files <- unlist(strsplit(input_args[8],","))

# if you provide less names than the number of assoc files, use arbitrary names
if (length(assoc.names) < length(assoc.files)){
  assoc.names = as.character(seq(1,length(assoc.files)))
}

# load metal results
metal.data <- fread(metal.file,data.table=F)

# loop through results
for (f in seq(1,length(assoc.files)){
  assoc.data <- fread(assoc.files[f],data.table=F)
  assoc.data = assoc.data[,c(marker.column, freq.column, pval.column, cols.tokeep)]
  colnames(f.data)[2:length(f.data[1,])] <- sub("^",paste(assoc.names[f],"_",sep=""),colnames(assoc.data)[2:length(f.data[1,])])
  metal.data <- merge(metal.data, f.data[,cols.tokeep], by.x=c("MarkerName"), by.y=c(marker.column))
}

# order based on meta pvalue
metal.data = metal.data[order(metal.data[,"P-value"]),]

# write results out to file
fwrite(metal.data, file = paste(out.pref,"_all.csv",sep=""), sep=",")

png(filename = paste(out.pref,"_all_plots.png",sep=""),width = 11, height = 11, units = "in", res=400, type = "cairo")
par(mfrow=c(3,3))
qq(metal.data[,"P-value"],main="All variants")
qq(metal.data[metal.data$MAF>0.05,"P-value"],main="Variants with MAF>0.05")
qq(metal.data[metal.data$MAF<=0.05,"P-value"],main="Variants with MAF<=0.05")
manhattan(metal.data,chr="chr",bp="pos",p="P-value", main="All variants")
manhattan(metal.data[metal.data$MAF>0.05,],chr="chr",bp="pos",p="P-value", main="Variants with MAF>=0.05")
manhattan(metal.data[metal.data$MAF<=0.05,],chr="chr",bp="pos",p="P-value", main="Variants with MAF<=0.05")
dev.off()
