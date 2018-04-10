args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
expr.file <- args[2]
exon.file <- args[3]
tfbs.file <- args[4]
states.file <- args[5]
ptv.file <- args[6]
genh.file <- args[7]
states.names <- c("1_Active_TSS","9_Active_enhancer_1","10_Active_enhancer_2")
genes.pad <-5000
minmaf <- 0.01
outpref <- args[8]

####################################################################################################################################
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr10.pass_and_fail.gtonly.minDP10.gds"
# expr.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/freeze5b/mask2/t2dreamdb_rnaseq_FPKM2_transcripts_03262018.csv"
# exon.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/freeze5b/mask2/t2dreamdb_rnaseq_FPKM2_exons_03262018.csv"
# tfbs.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/freeze5b/mask2/all_tfbs_chr10.hg38.bed"
# states.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/freeze5b/mask2/Islets.chromatinStates.hg38.v2.bed"
# ptv.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/freeze5b/mask2/freeze5b_dp10_ptv_mask1_chr10.tsv"
# genh.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/genehancer/gene_to_enhancer.RData"
# states.names <- c("1_Active_TSS","9_Active_enhancer_1","10_Active_enhancer_2")
# genes.pad <-5000
# chr = "chr10"
# minmaf <- 0.01
####################################################################################################################################



library(SeqVarTools)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(data.table)
library(SeqArray)
library(biomaRt)

getAggList <- function(gds, variants.orig){
  filtOrig <- seqGetFilter(gds)
  seqSetFilter(gds, variant.id=variants.orig)
  variants.new <- .expandAlleles(gds)
  group <- data.frame(variant.id=variants.new$variant.id, allele.index=variants.new$allele.index)
}
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

spl <- unlist(strsplit(gds.file,"\\."))
chr <- spl[startsWith(spl,"chr")]

# load genes
genes.pan <- fread(expr.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
genes.pan$chromosome_name <- sub("^", "chr", genes.pan$chromosome_name)
genes.pan <- genes.pan[genes.pan$chromosome_name == chr,]
genes.pan <- genes.pan[!duplicated(genes.pan$transcript_id),]
genes.pan <- genes.pan[!is.na(genes.pan$transcript_id),]
genes.gr <- GRanges(seqnames=genes.pan$chromosome_name,ranges=IRanges(start=genes.pan$transcript_start - genes.pad, end=genes.pan$transcript_end + genes.pad), trans_id=genes.pan$transcript_id, gene_id=genes.pan$ensembl_gene_id, symbol=genes.pan$hgnc_symbol)

# load exons
genes.exons <- fread(exon.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
genes.exons <- genes.exons[genes.exons$transcript_id %in% genes.gr$trans_id,]
genes.exons.gr <- GRanges(seqnames=rep(chr,NROW(genes.exons)),ranges=IRanges(start=genes.exons$exon_start, end=genes.exons$exon_end), transcript_id=genes.exons$transcript_id)

# load ptvs
ptv.raw <- fread(ptv.file,data.table=F)
ptv.gr <- GRanges(seqnames=ptv.raw$chromosome,ranges=IRanges(start=ptv.raw$position, end=ptv.raw$position+1),symbol=ptv.raw$group_id,trans_id=ptv.raw$transcript, ref = ptv.raw$ref, alt = ptv.raw$alt, vep_csq = ptv.raw$vep_csq)

## load chromatin states
states.raw <- fread(states.file,data.table=F, sep="\t")
states.raw = states.raw[states.raw$V1 == chr,]
states.subset <- states.raw[states.raw$V4 %in% states.names,] # 4216
states.gr <- GRanges(seqnames=states.subset$V1,ranges=IRanges(start=states.subset$V2,end=states.subset$V3),state=states.subset$V4)

## load tfbs data
tfbs.raw <- fread(tfbs.file, data.table=F)
# tfbs.raw$V1 <- sub("^","chr",tfbs.raw$V1)
tfbs.raw <- tfbs.raw[tfbs.raw$V1 == chr,]
tfbs.gr <- GRanges(seqnames=tfbs.raw$V1,ranges=IRanges(start=tfbs.raw$V2,end=tfbs.raw$V3),tf=tfbs.raw$V4)

# get tfbs in promoter
tfbs.prom.id <- findOverlaps(tfbs.gr,states.gr[states.gr$state == "1_Active_TSS",])
tfbs.prom.gr <- tfbs.gr[queryHits(tfbs.prom.id),]
tfbs.prom.gr$state <- states.gr[states.gr$state == "1_Active_TSS",]$state[subjectHits(tfbs.prom.id)]
tfbs.prom.gr$tfbs_state <- tfbs.prom.gr$state

# subset tfbs by enhancer
tfbs.enh.id <- findOverlaps(tfbs.gr,states.gr[states.gr$state %in% c("9_Active_enhancer_1","10_Active_enhancer_2"),])
tfbs.enh.gr <- tfbs.gr[queryHits(tfbs.enh.id),]
tfbs.enh.gr$state <- states.gr[states.gr$state %in% c("9_Active_enhancer_1","10_Active_enhancer_2"),]$state[subjectHits(tfbs.enh.id)]
tfbs.enh.gr$tfbs_state <- tfbs.enh.gr$state

# load genehnacer data
load(genh.file)

#enhancers to granges
genh.gr <- GRanges(seqnames = genh1$chrom, ranges = IRanges(start = genh1$start, end = genh1$end), enh_id = row.names(genh1))

#subset tfbs in enhancers to those in genehancers
tfbs.enh.genh.id <- findOverlaps(tfbs.enh.gr,genh.gr)
tfbs.enh.genh.gr <- tfbs.enh.gr[queryHits(tfbs.enh.genh.id),]
tfbs.enh.genh.gr$enh_id <- genh.gr$enh_id[subjectHits(tfbs.enh.genh.id)]

## open gds file
gds.data <- seqOpen(gds.file)

## collecting the groups
groups <- list()

## loop through genes
for (gind in seq(1,NROW(genes.gr))){
  # for (gind in seq(1,5)){
    ptv.flag <- F
    prom.flag <- F
    enh.flag <- F
  
  print(paste(gind,"/",length(genes.gr[,1]),sep=""))
  
  ## set cur trans
  cur_trans.gr <- genes.gr[gind,]
  cur_trans.symbol <- genes.gr$symbol[gind]
  cur_trans.trans_id <- genes.gr$trans_id[gind]
  
  # get exon coordinates
  cur_trans.exons.gr <- genes.exons.gr[genes.exons.gr$transcript_id == cur_trans.trans_id,]
  if (NROW(cur_trans.exons.gr) == 0){
    next
  }
  
  # get ptv variant in gene
  cur_trans.ptv.ovp <- findOverlaps(ptv.gr,cur_trans.exons.gr)
  cur_trans.ptv.gr <- ptv.gr[queryHits(cur_trans.ptv.ovp),]
  if (NROW(cur_trans.ptv.gr) > 0){
    ptv.flag = T
    cur_trans.ptv.gr$tfbs_state <- "NA"
    mcols(cur_trans.ptv.gr) <- data.frame(tfbs_state = cur_trans.ptv.gr$tfbs_state, ptv = cur_trans.ptv.gr$vep_csq)
  }
  
  # get ranges in prom-enhancer
  cur_trans.prom.ids <- findOverlaps(tfbs.prom.gr,cur_trans.gr)
  cur_trans.prom.gr <- tfbs.prom.gr[queryHits(cur_trans.prom.ids),]
  if (NROW(cur_trans.prom.gr) > 0){
    prom.flag = T
    cur_trans.prom.gr$ptv <- "NA"
    mcols(cur_trans.prom.gr) <- data.frame(tfbs_state=cur_trans.prom.gr$tfbs_state,ptv=cur_trans.prom.gr$ptv)
  }
  
  # get ranges for enahncers mapped to gene
  cur_trans.enh.ids <- unique(gene_to_enh[[cur_trans.gr$symbol]])
  cur_trans.enh.gr <- tfbs.enh.genh.gr[tfbs.enh.genh.gr$enh_id %in% cur_trans.enh.ids,]
  if (NROW(cur_trans.enh.gr) > 0){
    enh.flag = T
    cur_trans.enh.gr$ptv <- "NA"
    mcols(cur_trans.enh.gr) <- data.frame(tfbs_state=cur_trans.enh.gr$tfbs_state,ptv=cur_trans.enh.gr$ptv)
  }
  
  # combine ptv, prom, and enhancer
  if (ptv.flag & prom.flag & enh.flag){
    cur_trans.var_regions <- do.call(c,list(cur_trans.ptv.gr,cur_trans.prom.gr,cur_trans.enh.gr))
  } else if (ptv.flag & prom.flag){
    cur_trans.var_regions <- do.call(c,list(cur_trans.ptv.gr,cur_trans.prom.gr))
  } else if (ptv.flag & enh.flag){
    cur_trans.var_regions <- do.call(c,list(cur_trans.ptv.gr,cur_trans.enh.gr))
  } else if (ptv.flag) {
    cur_trans.var_regions <- cur_trans.ptv.gr
  } else if (prom.flag){
    cur_trans.var_regions <- cur_trans.prom.gr
  } else if (enh.flag) {
    cur_trans.var_regions <- cur_trans.enh.gr
  } else {
    next
  }
  
  ## subset gds
  seqSetFilter(gds.data,cur_trans.var_regions)
  
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
  gds.gr <- GRanges(seqnames=sub("^","chr",gds.df$chromosome), ranges=IRanges(start=gds.df$position,end=gds.df$position+1), position=gds.df$position, ref=gds.df$ref, allele = gds.df$allele, nAllelles = gds.df$nAlleles, allele.index = gds.df$allele.index, variant.id = gds.df$variant.id, maf = gds.df$maf)
  
  # add annotation
  gds.gr.ovp.id <- findOverlaps(gds.gr,cur_trans.var_regions)
  gds.gr.ovp <- gds.gr[queryHits(gds.gr.ovp.id),]
  gds.gr.ovp$ptv <- cur_trans.var_regions$ptv[subjectHits(gds.gr.ovp.id)]
  gds.gr.ovp$tfbs_state <- cur_trans.var_regions$tfbs_state[subjectHits(gds.gr.ovp.id)]
  
if(length(gds.gr.ovp$variant.id) == 0){
    next
}

  gds.gr.df <- data.frame(variant.id = gds.gr.ovp$variant.id, position = gds.gr.ovp$position, chromosome = seqnames(gds.gr.ovp), ref = gds.gr.ovp$ref, allele = gds.gr.ovp$allele, nAlleles = gds.gr.ovp$nAllelles, allele.index = gds.gr.ovp$allele.index, maf = gds.gr.ovp$maf, ptv = gds.gr.ovp$ptv, tfbs_state = gds.gr.ovp$tfbs_state, group_id = cur_trans.trans_id)
  gds.gr.df <- gds.gr.df[!duplicated(gds.gr.df$variant.id),]
  # ptv.final.id <- merge(gds.gr.df[gds.gr.df$ptv != "NA",],cur_trans.ptv.df, by.x = c("chromosome","position","ref","allele"), by.y = c("chromosome","position","ref","alt"))$variant.id
  
  # gds.gr.df <- rbind(gds.gr.df[gds.gr.df$ptv == "NA",],gds.gr.df[gds.gr.df$variant.id %in% ptv.final.id,])
  
  # gds.final.df <- data.frame()
  # for (v in unique(gds.gr.df$variant.id)){
  #   new.tfbs_states <- paste(unique(gds.gr.df[gds.gr.df$variant.id == v, "tfbs_state"]),collapse=",")
  #   new.ptv <- paste(unique(gds.gr.df[gds.gr.df$variant.id == v, "ptv"]),collapse=",")
  #   new.row <- gds.gr.df[gds.gr.df$variant.id == v,][1,]
  #   new.row$tfbs_state <- new.tfbs_states
  #   new.row$ptv <- new.ptv
  #   gds.final.df <- rbind(gds.final.df,new.row)
  # }
  # 
  
  gds.gr.df$annotation <- paste(gds.gr.df$ptv,gds.gr.df$tfbs_state,sep="_")
  gds.gr.df$annotation <- sub("NA_","",gds.gr.df$annotation)
  gds.gr.df$annotation <- sub("_NA","",gds.gr.df$annotation)
  gds.gr.df <- gds.gr.df[,c("group_id","variant.id","position","chromosome","ref","allele","nAlleles","allele.index","maf","annotation")]
  # groups[[cur_trans.symbol]] <- gds.final.df
  groups[[cur_trans.trans_id]] <- gds.gr.df
  
}

## close the geno file
seqClose(gds.data)

## remove any groups that dont have loci
groups <- groups[!sapply(groups, is.null)]
groups.final <- do.call(rbind,groups)


## save groups
# save(groups, "groups", file="/Users/tmajaria/Documents/projects/topmed/code/varshney/freeze5b/mask3/mask3.RData")
# fwrite(groups.final, file="/Users/tmajaria/Documents/projects/topmed/code/varshney/freeze5b/mask3/mask3.csv",sep=",")
save(groups, "groups", file=paste(outpref,'mask3.RData',sep='.'))
fwrite(groups.final, file=paste(outpref,'mask3.csv',sep='.'),sep=",")
