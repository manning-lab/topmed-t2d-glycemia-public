
library(data.table)
install.packages("qqman",repos='http://cran.us.r-project.org')  #_0.14.tar.gz")
library(qqman)

args <- commandArgs(trailingOnly=T)
results.file <- args[1]
groups.file <- args[2]
state.file <- args[3]
out.file_pref <- args[4]
state.names <- c("9_Active_enhancer_1","10_Active_enhancer_2","1_Active_TSS")


# results.file <- "/Users/tmajaria/Documents/projects/topmed/results/varshney/genomewide_v2/freeze4_t2d.assoc.RData"
# results.file <- "/Users/tmajaria/Documents/projects/tfbs_chr1_test.assoc.RData"
# groups.file <- "/Users/tmajaria/Documents/projects/tfbs_chr1_groups.RData"
# state.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/local_cs_states/Islets.chromatinStates.bed"
# tfbs.file <- "/Users/tmajaria/Documents/projects/topmed/data/varshney/all_tfbs_abcu/all_tfbs_chr1.csv"
# out.file_pref <- "tfbs_v2"

# results.raw <- fread(results.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
load(results.file)
results.raw <- assoc
results.res <- results.raw$results[!is.na(results.raw$results$n.sample.alt),]
results.top <- results.res[results.res$pval_0 < .005,]

load(groups.file)
# groups.top <- groups[names(groups) %in% results.top$V1]
groups.top <- results.raw$variantInfo[names(results.raw$variantInfo) %in% row.names(results.top)]
# names(groups.top)

if (length(groups.top) == 0){
  pdf(paste(out.file_pref,"_top_hits",".pdf",sep=""),width=11)
  dev.off()
} else {

library(GenomicRanges)
# source("https://bioconductor.org/biocLite.R")
# biocLite("STAN")
library(STAN)
library(Gviz)
library(biomaRt)


grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl <- select(ensembl_75,keys = groups.top[[1]]$chr[1],keytype = "chromosome_name", columns =  c( "ensembl_transcript_id","hgnc_symbol", "chromosome_name","transcript_start","transcript_end", "start_position", "end_position", "ensembl_gene_id"))


# gene.raw <- fread(gene.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
# gene.prep <- ensembl[ensembl$ensembl_gene_id %in%  gene.raw[gene.raw$RPKM>=2 & gene.raw$Gene_Type=="protein_coding",1],] 
# gene.top <- gene.prep[gene.prep$hgnc_symbol %in% names(groups.top),]
gene.top <- ensembl[ensembl$ensembl_transcript_id %in% names(groups.top),]

state.raw <- fread(state.file,sep="\t",header=F,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
state.subset <- state.raw[which(state.raw[,4] %in% state.names),] # 14363
state.subset[,1] <- sub("chr","",state.subset[,1]) # 1329
state.gr <- GRanges(seqnames=state.subset[,1],ranges=IRanges(start=state.subset[,2],end=state.subset[,3]),state=state.subset[,4])
head(state.gr)

tfbs.unparsed <- fread(tfbs.file, sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
head(tfbs.unparsed)
# tfbs.gr_o <- GRanges(seqnames=sub("^","chr",tfbs.unparsed$V1),ranges=IRanges(start=tfbs.unparsed$V2-100,end=tfbs.unparsed$V3+100),tf=tfbs.unparsed$V4)
tfbs.gr_o <- GRanges(seqnames=sub("chr","",tfbs.unparsed$V1),ranges=IRanges(start=tfbs.unparsed$V2,end=tfbs.unparsed$V3),tf=tfbs.unparsed$V4)
rm(tfbs.unparsed)
tfbs.phits <- findOverlaps(tfbs.gr_o,state.gr)
tfbs.gr <- tfbs.gr_o[queryHits(tfbs.phits),]

pdf(paste(out.file_pref,"_top_hits",".pdf",sep=""),width=11)
layout(matrix(seq(length(groups.top)),nrow=length(groups.top),ncol=1,byrow=T))

for (j in seq(1,length(groups.top))){
  
  # cur_group <- groups[[names(groups.top)[j]]]
  cur_group <- groups.top[[j]]
  all_group_info <- groups[[names(groups.top)[j]]]
  cur_gene <- gene.top[gene.top$ensembl_transcript_id==names(groups.top)[j],]
  
  var.gr <- GRanges(seqnames=cur_group$chr, ranges=IRanges(start=cur_group$pos ,end=cur_group$pos), maf=-log10(cur_group$freq))
  y_bounds <- c(floor(min(var.gr$maf)),ceiling(max(var.gr$maf)))
  
  cur_genes.gr <- GRanges(seqnames=cur_gene$chromosome_name,ranges=IRanges(start=cur_gene$start_position-5000,end=cur_gene$end_position+5000),gene_id=cur_gene$ensembl_transcript_id)
  
  ovp <- findOverlaps(state.gr,cur_genes.gr)
  ovp
  cur_state.gr <- state.gr[queryHits(ovp),]
  
  ovp <- findOverlaps(tfbs.gr,cur_genes.gr)
  ovp
  cur_tfbs.gr <- tfbs.gr[queryHits(ovp),]
  
  from <- min(min(start(cur_genes.gr)),min(start(var.gr))) - 15000
  to <- max(max(end(cur_genes.gr)),max(start(var.gr))) + 15000
  
  
  # annos <- c("splice_acceptor_variant","splice_donor_variant","splice_region_variant","stop_gained","stop_lost", "start_gained", "start_lost", "frameshift_variant")
  # sav <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "splice_acceptor_variant" %in% unlist(strsplit(as.character(x),",")))),]
  # sdv <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "splice_donor_variant" %in% unlist(strsplit(as.character(x),",")))),]
  # srv <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "splice_region_variant" %in% unlist(strsplit(as.character(x),",")))),]
  # sg <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "stop_gained" %in% unlist(strsplit(as.character(x),",")))),]
  # sl <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "stop_lost" %in% unlist(strsplit(as.character(x),",")))),]
  # srg <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "start_gained" %in% unlist(strsplit(as.character(x),",")))),]
  # srl <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "start_lost" %in% unlist(strsplit(as.character(x),",")))),]
  # fv <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) "frameshift_variant" %in% unlist(strsplit(as.character(x),",")))),]
  chr <- groups.top[[1]]$chr[1]
  gen <- "hg19"
  
  # sav_track <- DataTrack(data=-log10(sav$maf), start=sav$position,end=sav$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # sdv_track <- DataTrack(data=-log10(sdv$maf), start=sdv$position,end=sdv$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # srv_track <- DataTrack(data=-log10(srv$maf), start=srv$position,end=srv$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # sg_track <- DataTrack(data=-log10(sg$maf), start=sg$position,end=sg$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # # sl_track <- DataTrack(data=-log10(sl$maf), start=sl$position,end=sl$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # # srg_track <- DataTrack(data=-log10(srg$maf), start=srg$position,end=srg$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # srl_track <- DataTrack(data=-log10(srl$maf), start=srl$position,end=srl$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # fv_track <- DataTrack(data=-log10(fv$maf), start=fv$position,end=fv$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  # 
  # displayPars(sav_track) <- list(groups=factor("sav", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="red",legend=T,box.legend=F,cex.legend=0.8)
  # displayPars(sdv_track) <- list(groups=factor("sdv", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="yellow",legend=T,box.legend=F,cex.legend=0.8)
  # displayPars(srv_track) <- list(groups=factor("srv", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="green",legend=T,box.legend=F,cex.legend=0.8)
  # displayPars(sg_track) <- list(groups=factor("sg", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="blue",legend=T,box.legend=F,cex.legend=0.8)
  # # displayPars(sl_track) <- list(groups=factor("sl", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="orange",legend=TRUE,box.legend=F,cex.legend=0.8)
  # # displayPars(srg_track) <- list(groups=factor("srg", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="pink",legend=TRUE,box.legend=F,cex.legend=0.8)
  # displayPars(srl_track) <- list(groups=factor("srl", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="purple",legend=T,box.legend=F,cex.legend=0.8)
  # displayPars(fv_track) <- list(groups=factor("fv", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="cyan",legend=T,box.legend=F,cex.legend=0.8)
  # 
  # reg <- all_group_info[unlist(lapply(all_group_info$annotation,function(x) x %in% c("9_Active_enhancer_1","10_Active_enhancer_2","1_Active_TSS"))),]
  
  # reg_track <- DataTrack(data=-log10(all_group_info$maf), start=all_group_info$position,end=all_group_info$position+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF")
  
  # displayPars(reg_track) <- list(groups=factor("reg", levels=c("sav","sdv","srv","sg","sl","srg","srl","fv","reg")), col=c("red","yellow","green","blue","orange","pink","purple","cyan","black"),fill="black",legend=T,box.legend=FALSE,cex.legend=0.8)
  
  # ovt <- OverlayTrack(trackList = list(reg_track,sav_track,sdv_track,srv_track,sg_track,srl_track,fv_track),legend=T,box.legend=FALSE, ylim=y_bounds, background.title="darkred", name="-log MAF",cex.legend=0.8)
  dtrack <- DataTrack(data= var.gr$maf, start=start(var.gr),end=end(var.gr)+1,chromosome = chr, genome = gen,ylim=y_bounds, name = "-logMAF",background.title="orangered4")
  # ovt <- OverlayTrack(trackList = list(reg_track,sav_track,sdv_track,srv_track,sg_track,srl_track,fv_track),legend=T,box.legend=FALSE, ylim=y_bounds, background.title="darkred", name="-log MAF",cex.legend=0.8)
  
  # snpregtrack <- HighlightTrack(trackList = list(ovt), start = start(cur_state.gr), end=end(cur_state.gr),chromosome = groups.top[[1]]$chr[1], col="lightgrey", fill="lightgrey", alpha=0.7, inBackground=TRUE, background.title="orangered4")
  
  # snpregtrack2 <- HighlightTrack(trackList = list(snpregtrack), start = start(cur_tfbs.gr), end=end(cur_tfbs.gr),chromosome = groups.top[[1]]$chr[1], col="red", fill="red", alpha=0.7, inBackground=TRUE, background.title="orangered4")
  snpregtrack2 <- HighlightTrack(trackList = list(dtrack), start = start(cur_tfbs.gr), end=end(cur_tfbs.gr),chromosome = groups.top[[1]]$chr[1], col="orangered4", fill="orangered4", alpha=0.3, inBackground=TRUE, background.title="orangered4")
  
  # plotTracks(snpregtrack)
  
  axTrack <- GenomeAxisTrack(genome = "hg19" , chromosome = groups.top[[1]]$chr[1])
  
  idxTrack <- IdeogramTrack(genome = "hg19" , chromosome = groups.top[[1]]$chr[1])
  
  # plotTracks(list(axTrack,idxTrack,snpregtrack))
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
  
  gtrack <- BiomartGeneRegionTrack(start = from, end = to, biomart=ensembl_75, strand="+-", genome = "hg19", chromosome=groups.top[[1]]$chr[1], name = "Islet", showId=TRUE, geneSymbol=FALSE, filters = list("ensembl_transcript_id"=cur_genes.gr$gene_id), background.title="cadetblue4")
  snpregtrack <- HighlightTrack(trackList = list(gtrack), start = start(cur_state.gr), end=end(cur_state.gr),chromosome = groups.top[[1]]$chr[1], col="lightgrey", fill="lightgrey", alpha=0.7, inBackground=TRUE, background.title="orangered4")
  # plotTracks(c(list(axTrack,idxTrack,snpregtrack,gtrack),cs.track), showTitle = TRUE, sizes=c(2,1,7,8,1), from=from,to=to)
  # plotTracks(c(list(axTrack,idxTrack,snpregtrack,snpregtrack2,gtrack),cs.track), showTitle = TRUE, sizes=c(2,1,7,7,8,1), from=from,to=to)
  plotTracks(c(list(axTrack,idxTrack,snpregtrack2,snpregtrack),cs.track), showTitle = TRUE, sizes=c(2,1,7,2,1), from=from,to=to)
}

dev.off()
}
