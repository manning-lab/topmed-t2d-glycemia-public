library(SeqArray)

## testing inputs:
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/code/testing_inputs/singleVariantFull/freeze4.chr21.pass.gtonly.minDP10.genotypes.gds"
# gds.subset.ids <- c("NWD146248","NWD146267","NWD146274","NWD146380")
# gds.subset.label <- "test"
## real inputs
args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
gds.subset.ids <- unlist(strsplit(args[2],","))
gds.subset.file <- args[3]

# load the gds file
gds.data <- seqOpen(gds.file)

# subset data
seqSetFilter(gds.data, sample.id=gds.subset.ids)

# save the subset
seqExport(gds.data, gds.subset.file)

# close the gds file
seqClose(gds.data)

# clean up
cleanup.gds(gds.subset.file)