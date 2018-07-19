
args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr10.pass_and_fail.gtonly.minDP10.gds"

library(SeqVarTools)
library(dplyr)
library(tidyr)
library(data.table)

.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             filter=seqGetData(gds.data, "annotation/filter"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    group_by_("variant.id") %>%
    as.data.frame()
}

gds.data <- seqOpen(gds.file)
df <- .expandAlleles(gds.data)
seqClose(gds.data)

chr <- paste0("chr",df$chromosome[1])
df <- df[df$filter == "PASS",]
df$MarkerName <- paste("chr", df$chromosome, "-", df$position, "-", df$ref, "-", df$alt, sep = "")
df <- df[, "MarkerName"]
df <- unique(df)

fwrite(list(df), file = paste0("freeze.5b.passing.variants.", chr, ".minDP10.csv"), sep = "\t", col.names = F, quote = F)

