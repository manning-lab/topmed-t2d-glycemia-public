
library(data.table)

args <- commandArgs(trailingOnly=T)
results.file <- args[1]
groups.file <- args[2]
state.file <- args[3]
out.file_pref <- args[4]
state.names <- unlist(strsplit(args[5],','))

load(results.file)
results.raw <- assoc
results.res <- results.raw$results[!is.na(results.raw$results$n.sample.alt),]
results.top <- results.res[results.res$pval_0 < .005,]

load(groups.file)
groups.top <- groups[names(groups) %in% results.top$V1]

if (length(groups.top) == 0){
  pdf(paste(out.file_pref,"_top_hits",".pdf",sep=""),width=11)
  dev.off()
} else {

library(GenomicRanges)
source("https://bioconductor.org/biocLite.R")
biocLite("STAN")
library(STAN)
library(Gviz)
library(biomaRt)


grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl <- select(ensembl_75,keys = groups.top[[1]]$chr[1],keytype = "chromosome_name", columns =  c( "ensembl_transcript_id","hgnc_symbol", "chromosome_name","transcript_start","transcript_end", "start_position", "end_position", "ensembl_gene_id"))

gene.top <- ensembl[ensembl$hgnc_symbol %in% names(groups.top),]

state.raw <- fread(state.file,sep="\t",header=F,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
state.subset <- state.raw[which(state.raw[,4] %in% state.names),] # 14363
state.subset[,1] <- sub("chr","",state.subset[,1]) # 1329
state.gr <- GRanges(seqnames=state.subset[,1],ranges=IRanges(start=state.subset[,2],end=state.subset[,3]),state=state.subset[,4])

pdf(paste(out.file_pref,"_top_hits",".pdf",sep=""),width=11)
layout(matrix(seq(length(groups.top)),nrow=length(groups.top),ncol=1,byrow=T))

for (j in seq(1,length(groups.top))){
  
  # cur_group <- groups[[names(groups.top)[j]]]
  cur_group <- groups.top[[j]]
  all_group_info <- groups[[names(groups.top)[j]]]
  cur_gene <- gene.top[gene.top$hgnc_symbol==names(groups.top)[j],]
  
  var.gr <- GRanges(seqnames=cur_group$chr, ranges=IRanges(start=cur_group$pos ,end=cur_group$pos), maf=-log10(cur_group$freq))
  y_bounds <- c(floor(min(var.gr$maf)),ceiling(max(var.gr$maf)))
  
  cur_genes.gr <- GRanges(seqnames=cur_gene$chromosome_name,ranges=IRanges(start=cur_gene$start_position-5000,end=cur_gene$end_position+5000),gene_id=cur_gene$ensembl_gene_id)
  
  ovp <- findOverlaps(state.gr,cur_genes.gr)
  ovp
  cur_state.gr <- state.gr[queryHits(ovp),]
  
  from <- min(min(start(cur_genes.gr)),min(start(var.gr))) - 15000
  to <- max(max(end(cur_genes.gr)),max(start(var.gr))) + 15000
  
  chr <- groups.top[[1]]$chr[1]
  gen <- "hg19"
  
  dtrack <- DataTrack(data= var.gr$maf, start=start(var.gr),end=end(var.gr)+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF",background.title="orangered4")
  
  snpregtrack <- HighlightTrack(trackList = list(ovt), start = start(cur_state.gr), end=end(cur_state.gr),chromosome = groups.top[[1]]$chr[1], col="lightgrey", fill="lightgrey", alpha=0.7, inBackground=TRUE, background.title="orangered4")
  
  axTrack <- GenomeAxisTrack(genome = "hg19" , chromosome = groups.top[[1]]$chr[1])
  
  idxTrack <- IdeogramTrack(genome = "hg19" , chromosome = groups.top[[1]]$chr[1])
  
  cs.raw <- state.raw
  cs.raw[,1] <- sub("chr","",state.raw[,1])
  cs.raw$V1 <- as.numeric(cs.raw$V1)
  cs.raw <- cs.raw[which(cs.raw[,1] == groups.top[[1]]$chr[1]),]
  cs.gr <- GRanges(seqnames=cs.raw[,1],ranges=IRanges(start=cs.raw[,2],end=cs.raw[,3]),state=cs.raw[,4],color=cs.raw[,5])
  cs.sub <- cs.gr[which(seqnames(cs.gr)==groups.top[[1]]$chr[1]),]
  cs.sub <- cs.sub[which(end(cs.sub)>from),]
  cs.sub <- cs.sub[which(start(cs.sub)<to),]
  
  clrs <- unique(cs.sub$color)
  clrs1 <- unique(cs.sub$state)
  cs.cols <- c()
  cs.colnames <- c()
  for (i in seq(1,length(clrs))){
    cc <- as.integer(unlist(strsplit(clrs[i],",")))
    cs.cols <- c(cs.cols,rgb(cc[1],cc[2],cc[3],255,maxColorValue = 255))
    cs.colnames <- c(cs.colnames,clrs1[i])
  }
  names(cs.cols) <- cs.colnames
  mcols(cs.sub) <- DataFrame(name=cs.sub$state)
  options(ucscChromosomeNames=TRUE)
  cs.track = viterbi2Gviz(cs.sub, groups.top[[1]]$chr[1], "hg19", from,to, cs.cols)
  
  gtrack <- BiomartGeneRegionTrack(start = from, end = to, biomart=ensembl_75, strand="+-", genome = "hg19", chromosome=groups.top[[1]]$chr[1], name = "Islet", showId=TRUE, geneSymbol=FALSE, filters = list("ensembl_gene_id"=cur_genes.gr$gene_id), background.title="cadetblue4")
  snpregtrack <- HighlightTrack(trackList = list(gtrack), start = start(cur_state.gr), end=end(cur_state.gr),chromosome = groups.top[[1]]$chr[1], col="lightgrey", fill="lightgrey", alpha=0.7, inBackground=TRUE, background.title="orangered4")
  
  plotTracks(c(list(axTrack,idxTrack,snpregtrack2,snpregtrack),cs.track), showTitle = TRUE, sizes=c(2,1,7,2,1), from=from,to=to)
}

dev.off()
}
