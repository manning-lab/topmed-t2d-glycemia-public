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

# load islet regulome states
islet.reg <- fread(reg.file, data.table = F, stringsAsFactors = F)

# subset by chr
islet.reg <- islet.reg[islet.reg$chr == chr,]

# subset by the states we want
islet.reg <- islet.reg[islet.reg$state %in% reg.states,]

# make into granges 
islet.reg.gr <- GRanges(seqnames = islet.reg$chr, ranges = IRanges(start = islet.reg$start, end = islet.reg$end), state = islet.reg$state)

## open gds file
gds.data <- seqOpen(gds.file)

## collecting the groups
groups <- list()

## loop through genes
for (hind in seq(1,length(hub.unique))){
  # for (hind in seq(1,5)){
  # hind = 166
  
  print(paste(hind,"/",length(hub.unique),sep=""))
  
  # set ptv and state flags
  enh.flag = F
  prom.flag = F
  
  # set cur gene id
  cur_hub.id <- hub.unique[hind]
  
  ## set cur gene
  cur_hub.df <- islet.hub[islet.hub$hub.id == cur_hub.id,]
  cur_hub.genes <- unique( unlist( lapply( unique( cur_hub.df$hub.genes ), function(x) unlist( strsplit( x , ";" ) ) ) ) )
  
  # to granges
  cur_hub.gr <- GRanges(seqnames = cur_hub.df$chr, ranges = IRanges(start = cur_hub.df$start, end = cur_hub.df$end), hub.id = cur_hub.df$hub.id, hub.genes = cur_hub.df$hub.genes)
  cur_hub.bait.gr <- GRanges(seqnames = cur_hub.df$bait.chr, ranges = IRanges(start = cur_hub.df$bait.start, end = cur_hub.df$bait.end), hub.id = cur_hub.df$hub.id, hub.genes = cur_hub.df$hub.genes)
  cur_hub.enh.gr <- GRanges(seqnames = cur_hub.df$enh.chr, ranges = IRanges(start = cur_hub.df$enh.start, end = cur_hub.df$enh.end), hub.id = cur_hub.df$hub.id, hub.genes = cur_hub.df$hub.genes)
  
  # get the enhancer regions in the hub, intersect with islet regulome data
  cur_hub.enh.ids <- findOverlaps(islet.reg.gr, cur_hub.enh.gr)
  cur_hub.enh.gr <- unique(islet.reg.gr[queryHits(cur_hub.enh.ids),])
  cur_hub.enh.gr <- cur_hub.enh.gr[cur_hub.enh.gr$state %in% enh.states,]
  
  if (NROW(cur_hub.enh.gr) > 0) {
    enh.flag = T
    mcols(cur_hub.enh.gr) <- data.frame(state = cur_hub.enh.gr$state, group_id = rep(cur_hub.id, NROW(cur_hub.enh.gr)))
  }
  
  # get the promoter regions in the hub, intersect with islet regulome data
  cur_hub.bait.ids <- findOverlaps(islet.reg.gr, cur_hub.bait.gr)
  cur_hub.bait.gr <- unique(islet.reg.gr[queryHits(cur_hub.bait.ids),])
  cur_hub.bait.gr <- cur_hub.bait.gr[cur_hub.bait.gr$state %in% prom.states,]
  
  if (NROW(cur_hub.bait.gr) > 0) {
    prom.flag = T
    mcols(cur_hub.bait.gr) <- data.frame(state = cur_hub.bait.gr$state, group_id = rep(cur_hub.id, NROW(cur_hub.bait.gr)))
  }
  
  # combine ptv, enhancer, and promoter ranges
  if (prom.flag & enh.flag){
    cur_hub.var_regions <- do.call(c,list(cur_hub.bait.gr,cur_hub.enh.gr))
  } else if (prom.flag){
    cur_hub.var_regions <- cur_hub.bait.gr
  } else if (enh.flag) {
    cur_hub.var_regions <- cur_hub.enh.gr
  } else {
    next
  }
  
  ## subset gds
  seqSetFilter(gds.data,cur_hub.var_regions)
  
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
  
  # get to granges
  gds.gr <- GRanges(seqnames=sub("^","chr",gds.df$chromosome), ranges=IRanges(start=gds.df$position,end=gds.df$position), position=gds.df$position, ref=gds.df$ref, allele = gds.df$allele, nAllelles = gds.df$nAlleles, allele.index = gds.df$allele.index, variant.id = gds.df$variant.id, maf = gds.df$maf)
  
  # add annotation
  gds.gr.ovp.id <- findOverlaps(gds.gr,cur_hub.var_regions)
  gds.gr.ovp <- gds.gr[queryHits(gds.gr.ovp.id),]
  gds.gr.ovp$state <- cur_hub.var_regions$state[subjectHits(gds.gr.ovp.id)]
  gds.gr.ovp$group_id <- cur_hub.var_regions$group_id[subjectHits(gds.gr.ovp.id)]
  
  gds.gr.df <- data.frame(variant.id = gds.gr.ovp$variant.id, position = gds.gr.ovp$position, chromosome = seqnames(gds.gr.ovp), ref = gds.gr.ovp$ref, allele = gds.gr.ovp$allele, nAlleles = gds.gr.ovp$nAllelles, allele.index = gds.gr.ovp$allele.index, maf = gds.gr.ovp$maf, state = gds.gr.ovp$state, group_id = gds.gr.ovp$group_id)
  gds.gr.df <- gds.gr.df[!duplicated(gds.gr.df$variant.id),]
  gds.gr.df <- gds.gr.df[,c("group_id","position","chromosome","ref","allele","maf","state")]
  names(gds.gr.df) <- c("group_id","position","chromosome","ref","alt","maf","state")
  gds.gr.df$genes <- paste(cur_hub.genes, collapse = ";")
  groups[[cur_hub.id]] <- gds.gr.df
  
}

## close the geno file
seqClose(gds.data)

## remove any groups that dont have variants
groups <- groups[!sapply(groups, is.null)]
groups.final <- do.call(rbind,groups)

fwrite(groups.final, file=paste0("freeze5b_dp10.", chr,".mask8.no.ptv.csv"),sep=",")
