# groupByGene_v5 <- function( gds.file , genes.all.file , genes.pan.file , anno.file , states.file , genes.pad , anno.value , states.names , minmaf ){

args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
genes.all.file <- args[2]
genes.pan.file <- args[3]
anno.file <- args[4]
states.file <- args[5]
chain.file <- args[6]
genes.pad <- 5000
anno.value <- c("splice_acceptor_variant","splice_donor_variant","splice_region_variant","stop_gained","stop_lost", "start_gained", "start_lost", "frameshift_variant")
states.names <- c("active_enhancer_1","active_enhancer_2","active_tss")
minmaf <- 0.01



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
  library(SeqArray)
  gds.data <- seqOpen(gds.file)
  chr_num <- seqGetData(gds.data,"chromosome")[1]
  chr <- paste("chr",chr_num,sep="")
  
  
  # genes.all.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/ensembl_genes.csv"
  # genes.pan.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/gtex/gtex.v6p.pancreas.expression.min.rpkm.0.1.txt"
  # gds.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/freeze4.chr10.pass.gtonly.minDP10.genotypes.gds"
  # genes.pad <- 5000
  # anno.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/annotations.subset.csv"
  # anno.value <- c("splice_acceptor_variant","splice_donor_variant","splice_region_variant","stop_gained","stop_lost", "start_gained", "start_lost", "frameshift_variant")
  # states.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/varshney_islets_chromatin_state.aggregation.chr10.csv"
  # state.names <- c("active_enhancer_1","active_enhancer_2","active_tss")
  # minmaf <- 0.01
  
  ## load genes first
  library(data.table)
  genes.all <- fread(genes.all.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  genes.all <- genes.all[order(genes.all$hgnc_symbol),]
  genes.pan.raw <- fread(genes.pan.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  genes.pan <- genes.all[genes.all$ensembl_gene_id %in%  genes.pan.raw[genes.pan.raw$RPKM>=2 & genes.pan.raw$Gene_Type=="protein_coding",1],]
  genes.pan <- genes.pan[!duplicated(genes.pan$ensembl_gene_id),]
  
  ## remap to hg19
  library(rtracklayer)
  hg38tohg19 <- import.chain(chain.file) # (http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/)
  genes.hg38 <- GRanges(seqnames=rep(chr,length(genes.pan$chromosome_name)),ranges=IRanges(start=genes.pan$start_position, end=genes.pan$end_position),symbol=genes.pan$hgnc_symbol,id=genes.pan$ensembl_gene_id)
  genes.hg19 <- unlist(liftOver(genes.hg38, hg38tohg19))
  
  library(GenomicRanges)
  genes.pan.gr <- GRanges(seqnames=rep(chr_num,length(genes.hg19$id)),ranges=IRanges(start=start(genes.hg19)-genes.pad, end=end(genes.hg19)+genes.pad),symbol=genes.hg19$symbol,id=genes.hg19$id)
  genes.pan.gr_nopad <- GRanges(seqnames=rep(chr_num,length(genes.hg19$id)),ranges=IRanges(start=start(genes.hg19), end=end(genes.hg19)),symbol=genes.hg19$symbol,id=genes.hg19$id)
  
  ## load annotations
  anno.raw <- fread(anno.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  
  ## load chromatin states
  states.raw <- fread(states.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  states.subset <- states.raw[which(states.raw[,7] %in% states.names),] # 14363
  states.subset[,4] <- sub("chr","",states.subset[,4]) # 1329
  states.gr <- GRanges(seqnames=states.subset$chr,ranges=IRanges(start=states.subset$start,end=states.subset$stop),state=states.subset$type)
  
  ## collecting the groups
  groups <- list()
  
  ## loop through genes
  for (gind in seq(1,length(genes.pan.gr[,1]))){
    print(paste(gind,"/",length(genes.pan.gr[,1]),sep=""))
    
    
    ## set cur gene
    cur_gene.gr <- genes.pan.gr[gind,]
    cur_gene.gr_nopad <- genes.pan.gr_nopad[gind,]
    
    ## subset gds
    seqSetFilter(gds.data,cur_gene.gr)
    
    ## get the variants to a better format
    library(SeqVarTools)
    library(dplyr)
    library(tidyr)
    gds.df <- .expandAlleles(gds.data)
    gds.df$maf <- alleleFrequency(gds.data,n=1)
    
    ## subset the variants to those with maf<1% and maf > 0
    gds.df.umaf <- gds.df[gds.df$maf<minmaf,]
    gds.df.umaf <- gds.df.umaf[gds.df.umaf$maf>0,]
    gds.gr <- GRanges(seqnames=gds.df.umaf$chromosome, ranges=IRanges(start=gds.df.umaf$position,end=gds.df.umaf$position))
    names(gds.gr) <- gds.df.umaf$variant.id
    
    ## subset annotations to our range
    anno.subset <- anno.raw[anno.raw$VEP_ensembl_Gene_ID == cur_gene.gr$id,]
    
    ## get only those annotations that have our consequences
    library(stringr)
    anno.subvals <- anno.subset[apply(as.matrix(anno.subset$VEP_ensembl_Consequence), 1, function(x) any(sapply(anno.value, function(y) str_detect(unlist(str_split(x, pattern=",")),y)))),]
    anno.df <- data.frame(pos=anno.subvals$pos,alt=as.vector(anno.subvals$alt),annotation=as.vector(anno.subvals$VEP_ensembl_Consequence))
    gds_match <- data.frame(pos=gds.df.umaf$position,alt=as.vector(gds.df.umaf$allele))
    gds_paste = do.call("paste", gds_match[,1:2])
    anno_paste = do.call("paste", anno.df[,1:2])
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
    
    ## get variants within the chromatin states
    states.subset <- states.gr[end(states.gr)>start(cur_gene.gr),]
    states.subset <- states.subset[start(states.subset)<end(cur_gene.gr),]
    states.hits <- findOverlaps(gds.gr,states.subset)
    states.snps <- gds.gr[queryHits(states.hits),]
    states.snps <- gds.df.umaf[queryHits(states.hits),]
    states.snps$annotation <- states.subset$state[subjectHits(states.hits)]
    
    ## now put all of the snps together for our gene
    gds.subset.df <- rbind(states.snps,anno.snps)
    
    groups[[cur_gene.gr$symbol]] <- gds.subset.df
    
  }
  
  ## close the geno file
  seqClose(gds.data)
  
  ## remove any groups that dont have loci
  groups <- groups[!sapply(groups, is.null)]
  
  ## save groups
  save(groups, "groups", file="groups.RData")
  
  ## return the groups
# }

