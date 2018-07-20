
args <- commandArgs(trailingOnly=T)
res.file <- args[1]
pass.file <- args[2]

out.file <- basename(res.file)

library(data.table)

# load files
res.data <- fread(res.file, data.table = F, stringsAsFactors = F)
pass.data <- fread(pass.file, data.table = F, stringsAsFactors = F, header = F)$V1

# keep passing variants
res.data <- res.data[res.data$MarkerName %in% pass.data,]

# write new file
fwrite(res.data, file = out.file, sep = ",", col.names = T, quote = F)