library(SeqArray)

## testing inputs:
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/code/testing_inputs/singleVariantFull/freeze4.chr21.pass.gtonly.minDP10.genotypes.gds"
# gds.subset.ids <- c("NWD146248","NWD146267","NWD146274","NWD146380")
# gds.subset.label <- "test"
## real inputs
args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
gds.sample.file <- args[2]
gds.variant.file <- args[3]
label <- args[4]

# get the right output file
gds.subset.file <- paste(label,"gds",sep=".")

# load the gds file
gds.data <- seqOpen(gds.file)

# filter samples
if(gds.sample.file != "none"){
	sample.ids <- unique(readLines(gds.sample.file))
	seqSetFilter(gds.data, sample.id=sample.ids)
}

# filter variants
if(gds.variant.file != "none"){
	variant.ids <- unique(readLines(gds.variant.file))
	seqSetFilter(gds.data, variant.id=variant.ids, action="intersect")
}

# save the subset
seqExport(gds.data, gds.subset.file)

# close the gds file
seqClose(gds.data)

# clean up
cleanup.gds(gds.subset.file)