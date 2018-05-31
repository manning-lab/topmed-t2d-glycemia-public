# Make ferrer based groups

args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
hub.file <- args[2]
reg.file <- args[3]
chr <- args[4]
reg.states <- c("Active_enhancers_I", "Active_enhancers_I_CTCF", "Active_promoters", "Active_promoters_CTCF")
enh.states <- c("Active_enhancers_I", "Active_enhancers_I_CTCF")
prom.states <- c("Active_promoters", "Active_promoters_CTCF")
minmaf <- 0.01

####################################################################################################################################
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr10.pass_and_fail.gtonly.minDP10.gds"
# hub.file <- "/Users/tmajaria/Documents/projects/topmed/data/ferrer/Islet_annotations_FerrerLab/Islet_enhancer_hubs.grch38.csv"
# reg.file <- "/Users/tmajaria/Documents/projects/topmed/data/ferrer/Islet_annotations_FerrerLab/Islet_regulome.parsed.grch38.csv"
# reg.states <- c("Active_enhancers_I", "Active_enhancers_I_CTCF", "Active_promoters", "Active_promoters_CTCF")
# enh.states <- c("Active_enhancers_I", "Active_enhancers_I_CTCF")
# prom.states <- c("Active_promoters", "Active_promoters_CTCF")
# chr = "chr10"
# minmaf <- 0.01
####################################################################################################################################


# load packages and functions
library(SeqVarTools)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(data.table)
library(SeqArray)
library(biomaRt)

.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             nAlleles=seqNumAllele(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    rename_(allele="alt") %>%
    group_by_("variant.id") %>%
    mutate_(allele.index=~1:n()) %>%
    as.data.frame()
}
#########################
# load the islet hubs
islet.hub <- fread(hub.file, data.table = F, stringsAsFactors = F)
islet.hub <- islet.hub[islet.hub$chr == chr,]
hub.unique <- unique(islet.hub$hub.id)

# make hub granges
hub.gr <- GRanges(seqnames = islet.hub$chr, ranges = IRanges(start = islet.hub$start, end = islet.hub$end), hub_id = islet.hub$hub.id, hub_genes = islet.hub$hub.genes)

# make bait granges
bait.gr <- GRanges(seqnames = islet.hub$bait.chr, ranges = IRanges(start = islet.hub$bait.start, end = islet.hub$bait.end), hub_id = islet.hub$hub.id, hub_genes = islet.hub$hub.genes)

# make enh granges
enh.hub.gr <- GRanges(seqnames = islet.hub$enh.chr, ranges = IRanges(start = islet.hub$enh.start, end = islet.hub$enh.end), hub_id = islet.hub$hub.id, hub_genes = islet.hub$hub.genes)
enh.hub.gr <- unique(enh.hub.gr)

# load islet regulome states
islet.reg <- fread(reg.file, data.table = F, stringsAsFactors = F)

# subset by chr
islet.reg <- islet.reg[islet.reg$chr == chr,]

# subset by the states we want
islet.reg <- islet.reg[islet.reg$state %in% reg.states,]

# make promoter granges
islet.reg.prom <- islet.reg[islet.reg$state %in% prom.states,]
prom.gr <- GRanges(seqnames = islet.reg.prom$chr, ranges = IRanges(start = islet.reg.prom$start, end = islet.reg.prom$end), state = islet.reg.prom$state)

# intersect bait regions with promoters
prom.ovp.ids <- findOverlaps(prom.gr, bait.gr)
prom.gr <- prom.gr[queryHits(prom.ovp.ids),]
prom.gr$hub_id <- bait.gr[subjectHits(prom.ovp.ids),]$hub_id
prom.gr <- unique(prom.gr)

# intersect enhancers with reg enhancers to get state
islet.reg.enh <- islet.reg[islet.reg$state %in% enh.states,]
reg.enh.gr <- GRanges(seqnames = islet.reg.enh$chr, ranges = IRanges(start = islet.reg.enh$start, end = islet.reg.enh$end), state = islet.reg.enh$state)

# intersect enh regions with enhancers
enh.ovp.ids <- findOverlaps(reg.enh.gr, enh.hub.gr)
enh.gr <- reg.enh.gr[queryHits(enh.ovp.ids),]
enh.gr$hub_id <- enh.hub.gr[subjectHits(enh.ovp.ids),]$hub_id
enh.gr <- unique(enh.gr)

# combine prom+enh 
all.regions <- as.data.frame(do.call(c, list(prom.gr, enh.gr)))
all.regions <- all.regions[,c(1,2,3,6,7)]
names(all.regions) <- c("chr","start","end","state","hub_id")
all.regions$group_id <- paste(all.regions$hub_id, all.regions$state, all.regions$chr, all.regions$start, all.regions$end, sep = ";")

## open gds file
gds.data <- seqOpen(gds.file)

## collecting the groups
groups <- list()

## loop through genes
for (gind in seq(1,nrow(all.regions))){
  # for (hind in seq(1,5)){
  # hind = 166
  
  print(paste(gind,"/",nrow(all.regions),sep=""))
  
  # set cur region
  cur.df <- all.regions[gind,]
  
  # set cur gene id
  cur.id <- cur.df$group_id[1]
  
  # to granges
  cur.gr <- GRanges(seqnames = cur.df$chr, ranges = IRanges(start = cur.df$start, end = cur.df$end), state = cur.df$state, hub_id = cur.df$hub_id, group_id = cur.df$group_id)
  
  ## subset gds
  seqSetFilter(gds.data,cur.gr)
  
  # skip if no var
  var.ids <- seqGetData(gds.data,"variant.id")
  num_var <- length(var.ids)
  if(num_var==0){
    next
  }
  
  ## get the variants to a better format
  gds.df <- .expandAlleles(gds.data)
  ref.freq <- seqAlleleFreq(gds.data)
  gds.df$maf <- pmin(ref.freq, 1-ref.freq)
  
  ## subset the variants to those with maf<1% and maf > 0
  gds.df <- gds.df[gds.df$maf<minmaf,]
  gds.df <- gds.df[gds.df$maf>0.0,]

  # check if we have no var again
  if(nrow(gds.df)==0){
    next
  }
  
  # fix chromosome
  gds.df$chromosome <- sub("^","chr",gds.df$chromosome)
  
  # add state
  gds.df$state <- cur.gr$state
  
  # add group_id
  gds.df$group_id <- cur.gr$group_id
  
  # subset to my columns
  gds.df <- gds.df[,c("group_id","position","chromosome","ref","allele","maf","state")]
  names(gds.df) <- c("group_id","position","chromosome","ref","alt","maf","state")
  groups[[cur.id]] <- gds.df
  
}

## close the geno file
seqClose(gds.data)

## remove any groups that dont have variants
groups <- groups[!sapply(groups, is.null)]
groups.final <- do.call(rbind,groups)

fwrite(groups.final, file=paste0("freeze5b_dp10.", chr,".mask9.csv"),sep=",")
