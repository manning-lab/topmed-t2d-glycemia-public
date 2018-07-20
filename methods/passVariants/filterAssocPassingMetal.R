## for metal
args <- commandArgs(trailingOnly=T)
res.file <- args[1]
pass.file <- args[2]

library(data.table)

# load files
res.data <- fread(res.file, data.table = F, stringsAsFactors = F)
pass.data <- fread(pass.file, data.table = F, stringsAsFactors = F)

# fix markers
# chr10:112998590:C:T
pass.data$V2 <- gsub(":","-",pass.data$V2)

# keep passing variants
res.data <- res.data[res.data$MarkerName %in% pass.data$V2,]

# write new file
fwrite(res.data, file = res.file, sep = ",", col.names = T, quote = F)