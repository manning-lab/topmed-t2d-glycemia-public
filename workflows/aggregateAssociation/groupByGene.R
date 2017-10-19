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
# states.names <- c("active_enhancer_1","active_enhancer_2","active_tss")
states.names <- c("9_Active_enhancer_1","10_Active_enhancer_2","1_Active_TSS")
minmaf <- 0.01
library(SeqVarTools)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(SeqArray)
library(biomaRt)

# source("https://bioconductor.org/biocLite.R")
# biocLite("SeqVarTools")
# genes.all.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/ensembl_genes.csv"
# genes.pan.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/gtex/gtex.v6p.pancreas.expression.min.rpkm.0.1.txt"
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/freeze4.chr20.pass.gtonly.minDP10.genotypes.gds"
# # # gds.file <- "/Users/tmajaria/Documents/projects/topmed/results/varshney/genomewide_v1/group_check/freeze4.chr10.pass.gtonly.minDP10.genotypes.gds"
# # # genes.pad <- 5000
# anno.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/annotations.subset.csv"
# # anno.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/freezes_2a_3a_4.snp_indel.annotated.general20170422.subset.gz.chr20.csv"
# # # anno.value <- c("splice_acceptor_variant","splice_donor_variant","splice_region_variant","stop_gained","stop_lost", "start_gained", "start_lost", "frameshift_variant")
# states.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/local_cs_states/Islets.chromatinStates.bed"
# # # states.names <- c("active_enhancer_1","active_enhancer_2","active_tss")
# # minmaf <- 0.01
# chain.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/hg38ToHg19.over.chain"
# # 


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
  
  ## load genes first
  mart <- useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genes.all <- select(mart,keys = chr_num,keytype = "chromosome_name",
                    columns =  c( "ensembl_transcript_id","hgnc_symbol", "chromosome_name","transcript_start", "transcript_end", "ensembl_gene_id", "start_position","end_position"))
  # genes.all <- fread(genes.all.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  # genes.all <- genes.all[order(genes.all$hgnc_symbol),]
  genes.pan.raw <- fread(genes.pan.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  genes.pan <- genes.all[genes.all$ensembl_gene_id %in%  genes.pan.raw[genes.pan.raw$RPKM>=2 & genes.pan.raw$Gene_Type=="protein_coding",1],]
  genes.pan <- genes.pan[!duplicated(genes.pan$ensembl_gene_id),]

  rm(genes.all)
  rm(genes.pan.raw)
  
  ## remap to hg19
  
  hg38tohg19 <- import.chain(chain.file) # (http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/)
  genes.hg38 <- GRanges(seqnames=sub("^","chr",genes.pan$chromosome_name),ranges=IRanges(start=genes.pan$start_position, end=genes.pan$end_position),symbol=genes.pan$hgnc_symbol,id=genes.pan$ensembl_gene_id)
  genes.hg19 <- unlist(liftOver(genes.hg38, hg38tohg19))
  
  rm(hg38tohg19)
  rm(genes.hg38)
  
  genes.pan.gr <- GRanges(seqnames=rep(chr_num,length(genes.hg19$id)),ranges=IRanges(start=start(genes.hg19)-genes.pad, end=end(genes.hg19)+genes.pad),symbol=genes.hg19$symbol,id=genes.hg19$id)
  genes.pan.gr_nopad <- GRanges(seqnames=rep(chr_num,length(genes.hg19$id)),ranges=IRanges(start=start(genes.hg19), end=end(genes.hg19)),symbol=genes.hg19$symbol,id=genes.hg19$id)

  rm(genes.hg19)
  
  ## load annotations
  anno.raw <- fread(anno.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  # anno.raw <- anno.raw[,c(1:3,5,7)]
  print(head(anno.raw$VEP_ensembl_Gene_ID))
  # anno.raw <- anno.raw[order(anno.raw$VEP_ensembl_Gene_ID),]
  a <- anno.raw$VEP_ensembl_Gene_ID != "."
  anno.sub_one <- anno.raw[a,]

  rm(anno.raw)
  
  ## load chromatin states
  # states.raw <- fread(states.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  states.unparsed <- fread(states.file,sep="\t",header=F,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  print(head(states.unparsed))
  states.raw <- data.frame(name=states.unparsed[,4],start=states.unparsed[,2],stop=states.unparsed[,3], chr=states.unparsed[,1],tissue=rep('islets',length(states.unparsed[,1])),category=rep('chromatin_state',length(states.unparsed[,1])),type=states.unparsed[,4])
  states.subset <- states.raw[which(states.raw[,7] %in% states.names),] # 14363
  states.subset[,4] <- sub("chr","",states.subset[,4]) # 1329
  states.gr <- GRanges(seqnames=states.subset$chr,ranges=IRanges(start=states.subset$start,end=states.subset$stop),state=states.subset$type)
  
  rm(states.unparsed)
  rm(states.raw)
  rm(states.subset)

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
    
    num_var <- length(seqGetData(gds.data,"variant.id"))
    
    if(num_var==0){
      next
    }
    
    ## get the variants to a better format
    gds.df <- .expandAlleles(gds.data)
    # library(SeqVarTools)
    # gds.df$maf <- alleleFrequency(gds.data,n=1)
    ref.freq <- seqAlleleFreq(gds.data)
    gds.df$maf <- pmin(ref.freq, 1-ref.freq)

    ## subset the variants to those with maf<1% and maf > 0
    gds.df.umaf <- gds.df[gds.df$maf<minmaf,]
    gds.df.umaf <- gds.df.umaf[gds.df.umaf$maf>0,]
    gds.gr <- GRanges(seqnames=gds.df.umaf$chromosome, ranges=IRanges(start=gds.df.umaf$position,end=gds.df.umaf$position))
    names(gds.gr) <- gds.df.umaf$variant.id
    
    print(length(gds.df.umaf$variant.id))
    ## subset annotations to our range
    anno.subset <- anno.sub_one[anno.sub_one$VEP_ensembl_Gene_ID == cur_gene.gr$id,]
    head(anno.subset)
    print(length(anno.subset))
    
    ## get only those annotations that have our consequences
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
    states.subset_o <- states.gr[end(states.gr)>start(cur_gene.gr),]
    states.subset <- states.subset_o[start(states.subset_o)<end(cur_gene.gr),]
    states.hits <- findOverlaps(gds.gr,states.subset)
    states.snps_o <- gds.gr[queryHits(states.hits),]
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

