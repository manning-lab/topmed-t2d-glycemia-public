## for metal
args <- commandArgs(trailingOnly=T)
res.file <- args[1]
pass.file <- args[2]

res.file <- "/Users/tmajaria/Documents/projects/public_workflows/fineMap/paintor/test_inputs/t2d_Model4_cov_age_sex_topmedproject_ALLT2D_selfreportAncestry.METAL.assoc.csv"
pass.file <- "/Users/tmajaria/Documents/projects/topmed/code/topmed-t2d-glycemia-public/methods/passVariants/test_inputs/freeze.5b.passing.variants.all.minDP10.csv"

out.file <- basename(res.file)

library(data.table)

# load files
res.data <- fread(res.file, data.table = F, stringsAsFactors = F)
pass.data <- fread(pass.file, data.table = F, stringsAsFactors = F, header = F)$V1

# fix markers
# chr10:112998590:C:T
pass.data <- gsub("-",":",pass.data)

# keep passing variants
res.data <- res.data[res.data$MarkerName %in% pass.data,]

# write new file
fwrite(res.data, file = out.file, sep = ",", col.names = T, quote = F)