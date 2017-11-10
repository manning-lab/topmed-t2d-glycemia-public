args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
genes.pan.file <- args[2]
anno.file <- args[3]
tfbs.file <- args[4]
states.file <- args[5]
label <- args[6]


minmaf <- 0.01
genes.pad <-5000
anno.value <- c("splice_acceptor_variant","splice_donor_variant","splice_region_variant","stop_gained","stop_lost", "start_gained", "start_lost", "frameshift_variant")
states.names <- c("1_Active_TSS")
states.prom <- c("10_Active_enhancer_2","9_Active_enhancer_1")

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


## open gds file
gds.data <- seqOpen(gds.file)
chr_num <- seqGetData(gds.data,"chromosome")[1]
chr <- paste("chr",chr_num,sep="")

genes.pan <- fread(genes.pan.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
genes.pan.gr <- GRanges(seqnames=rep(chr_num,length(genes.pan$ensembl_gene_id)),ranges=IRanges(start=genes.pan$transcript_start-genes.pad, end=genes.pan$transcript_end+genes.pad),symbol=genes.pan$hgnc_symbol,trans_id=genes.pan$ensembl_transcript_id,id=genes.pan$ensembl_gene_id)
genes.pan.gr_nopad <- GRanges(seqnames=rep(chr_num,length(genes.pan$ensembl_gene_id)),ranges=IRanges(start=genes.pan$transcript_start, end=genes.pan$transcript_end),symbol=genes.pan$hgnc_symbol,trans_id=genes.pan$ensembl_transcript_id,id=genes.pan$ensembl_gene_id)
genes.pan.gr <- genes.pan.gr[!duplicated(genes.pan.gr$trans_id),]
genes.pan.gr_nopad <- genes.pan.gr_nopad[!duplicated(genes.pan.gr_nopad$trans_id),]

## load annotations
anno.raw <- fread(anno.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
a <- anno.raw[,anno.field] != "."
anno.sub_one <- anno.raw[a,]

rm(anno.raw)

## load chromatin states
states.unparsed <- fread(states.file,sep="\t",header=F,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
states.raw <- data.frame(name=states.unparsed[,4],start=states.unparsed[,2],stop=states.unparsed[,3], chr=states.unparsed[,1],tissue=rep('islets',length(states.unparsed[,1])),category=rep('chromatin_state',length(states.unparsed[,1])),type=states.unparsed[,4])
states.raw[,4] <- sub("chr","",states.raw[,4]) # 1329
states.subset <- states.raw[which(states.raw[,7] %in% states.names),] # 14363

states.gr <- GRanges(seqnames=states.subset$chr,ranges=IRanges(start=states.subset$start,end=states.subset$stop),state=states.subset$type)
states.gr <- states.gr[seqnames(states.gr) == chr_num,]

states.prom_sub <- states.raw[which(states.raw[,7] %in% states.prom),] # 14363
states.grprom <- GRanges(seqnames=states.prom_sub$chr,ranges=IRanges(start=states.prom_sub$start,end=states.prom_sub$stop),state=states.prom_sub$type)
states.grprom <- states.grprom[seqnames(states.grprom) == chr_num,]

rm(states.unparsed)
rm(states.raw)
rm(states.subset)

## load tfbs data
tfbs.unparsed <- fread(tfbs.file, sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
head(tfbs.unparsed)
# tfbs.gr_o <- GRanges(seqnames=sub("^","chr",tfbs.unparsed$V1),ranges=IRanges(start=tfbs.unparsed$V2-100,end=tfbs.unparsed$V3+100),tf=tfbs.unparsed$V4)
tfbs.gr_o <- GRanges(seqnames=sub("chr","",tfbs.unparsed$V1),ranges=IRanges(start=tfbs.unparsed$V2,end=tfbs.unparsed$V3),tf=tfbs.unparsed$V4)
rm(tfbs.unparsed)
tfbs.phits <- findOverlaps(tfbs.gr_o,states.grprom)
tfbs.gr <- tfbs.gr_o[queryHits(tfbs.phits),]
# tfbs.gr <- tfbs.gr[seqnames(tfbs.gr) == chr]

# genes.pan.gr <- genes.pan.gr[genes.pan.gr$id %in% anno.sub_one$VEP_ensembl_Gene_ID,]
# genes.pan.gr_nopad <- genes.pan.gr_nopad[genes.pan.gr_nopad$id %in% anno.sub_one$VEP_ensembl_Gene_ID,]

## collecting the groups
groups <- list()

## loop through genes
for (gind in seq(1,length(genes.pan.gr[,1]))){
  
  print(paste(gind,"/",length(genes.pan.gr[,1]),sep=""))
  
  ## set cur gene
  cur_gene.gr <- genes.pan.gr[gind,]
  cur_gene.gr_nopad <- genes.pan.gr_nopad[gind,]
  head(cur_gene.gr)
  ## subset gds
  seqSetFilter(gds.data,cur_gene.gr)
  
  num_var <- length(seqGetData(gds.data,"variant.id"))
  
  if(num_var==0){
    next
  }
  
  ## get the variants to a better format
  gds.df <- .expandAlleles(gds.data)
  ref.freq <- seqAlleleFreq(gds.data)
  gds.df$maf <- pmin(ref.freq, 1-ref.freq)
  
  ## subset the variants to those with maf<1% and maf > 0
  gds.df.umaf <- gds.df[gds.df$maf<minmaf,]
  gds.df.umaf <- gds.df.umaf[gds.df.umaf$maf>0,]
  gds.gr <- GRanges(seqnames=gds.df.umaf$chromosome, ranges=IRanges(start=gds.df.umaf$position,end=gds.df.umaf$position))
  names(gds.gr) <- gds.df.umaf$variant.id
  
  # print(length(gds.df.umaf$variant.id))
  # subset annotations to our range
  anno.subset <- anno.sub_one[anno.sub_one[,"VEP_ensembl_Gene_ID"] == cur_gene.gr$id,]
  # head(anno.subset)
  # anno.subset <- anno.sub_one
  # print(length(anno.subset[,1]))
  if (length(anno.subset[,1])>0){
    
    ## get only those annotations that have our consequences
    anno.subvals <- anno.subset[apply(as.matrix(anno.subset[,anno.field]), 1, function(x) any(sapply(anno.value, function(y) str_detect(unlist(str_split(x, pattern=",")),y)))),]
    anno.df <- data.frame(pos=anno.subvals$pos,ref=as.vector(anno.subvals$ref), alt=as.vector(anno.subvals$alt),annotation=as.vector(anno.subvals[,anno.field]))
    gds_match <- data.frame(pos=gds.df.umaf$position,ref=as.vector(gds.df.umaf$ref), alt=as.vector(gds.df.umaf$allele))
    gds_paste = do.call("paste", gds_match[,1:3])
    anno_paste = do.call("paste", anno.df[,1:3])
    anno_paste2 <- cbind(anno_paste,as.vector(anno.df$annotation))
    uv <- unique(anno_paste2[,1])
    anno_paste3 <- data.frame(uv, rep(0,length(uv)))
    for (j in seq(1,length(uv))){
      ap2 <- anno_paste2[anno_paste2[,1]==uv[j],2]
      ap3 <- c()
      for (k in seq(1,length(ap2))){
        ap3 <- c(ap3,unlist(str_split(ap2[k],",")))
      }
      anno_paste3[j,2] <- paste(ap3[!duplicated(ap3)],collapse=",")
    }
    
    anno.snps <- gds.df.umaf[gds_paste %in% anno_paste3[,1],]
    anno.snps$annotation <- anno_paste3[anno_paste3[,1] %in% gds_paste,2]
  }
  print(paste("Number of exonic snps: ",length(anno.snps[,1]),sep=""))
  
  ## get variants within the chromatin states
  states.subset_o <- states.gr[end(states.gr)>start(cur_gene.gr),]
  states.subset <- states.subset_o[start(states.subset_o)<end(cur_gene.gr),]
  states.hits <- findOverlaps(gds.gr,states.subset)
  states.snps_o <- gds.gr[queryHits(states.hits),]
  states.snps <- gds.df.umaf[queryHits(states.hits),]
  states.snps$annotation <- states.subset$state[subjectHits(states.hits)]
  print(paste("Number of tss snps: ",length(states.snps[,1]),sep=""))
  
  # findOverlaps(tfbs.gr,cur_gene.gr)
  tfbs.subset_o <- tfbs.gr[end(tfbs.gr)>start(cur_gene.gr),]
  tfbs.subset <- tfbs.subset_o[start(tfbs.subset_o)<end(cur_gene.gr),]
  tfbs.hits <- findOverlaps(gds.gr,tfbs.subset)
  tfbs.snps_o <- gds.gr[queryHits(tfbs.hits),]
  tfbs.snps <- gds.df.umaf[queryHits(tfbs.hits),]
  tfbs.snps$annotation <- tfbs.subset$tf[subjectHits(tfbs.hits)]
  tfbs.snps <- tfbs.snps[!duplicated(tfbs.snps$variant.id),]
  print(paste("Number of tfbs snps: ",length(tfbs.snps[,1]),sep=""))
  # print(length(tfbs.snps[,1]))
  
  ## now put all of the snps together for our gene
  if (length(anno.subset[,1])>0){
    gds.subset.df <- rbind(states.snps,anno.snps,tfbs.snps)
    print(paste("Total snps in group: ",length(gds.subset.df[,1])))
  } else {
    gds.subset.df <- states.snps
    print(paste("Total snps in group: ",length(gds.subset.df[,1])))
  }
  groups[[cur_gene.gr$trans_id]] <- rbind(gds.subset.df,tfbs.snps)
  
}

## close the geno file
seqClose(gds.data)

## remove any groups that dont have loci
groups <- groups[!sapply(groups, is.null)]

## save groups
save(groups, "groups", file=paste(label,"groups.RData",sep="_"))
