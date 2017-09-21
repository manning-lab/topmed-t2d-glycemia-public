
###########################################################
#         R script for visualization using Gviz           #
#         5/11/2017   Qing Liu & Xiaochen Lin             #
########################################################### 


### Step 2:  Read in tissue specific gene list (from GTEx, RPKM>0.1)
# adipose_sub <- read.table("gtex.v6p.adipose_subcutaneous.expression.min.rpkm.0.1.txt", h=T)
# adipose_vis <- read.table("gtex.v6p.adipose_visceral.expression.min.rpkm.0.1.txt", h=T)
# pancreas <- read.table("gtex.v6p.pancreas.expression.min.rpkm.0.1.txt", h=T)
# muscle <- read.table("gtex.v6p.muscle_skeletal.expression.min.rpkm.0.1.txt", h=T)
# liver <- read.table("gtex.v6p.liver.expression.min.rpkm.0.1.txt", h=T)

adipose_sub <- read.table(gtex.adipose.sub.file, h=T)
adipose_vis <- read.table(gtex.adipose.vis.file, h=T)
pancreas <- read.table(gtex.pancreas.file, h=T)
muscle <- read.table(gtex.muscle.file, h=T)
liver <- read.table(gtex.liver.file, h=T)



do.Qings.Gviz <- function(in.d,chr,chr_num,snp,P.value.column,MarkerID.column,plot=F) {
  #setwd("")
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("Gviz")
  #install.packages("Gviz")
  library(Gviz)
  library(biomaRt)
  #library(data.table)
  
  
  ### Step 1: Read in meta results (by chr) & Specify variant position
  #in.meta <- fread("meta_FG_common_chr2.txt",h=T)  
  #for_pos = sapply(as.character(in.meta$MarkerName),function(x) unlist(strsplit(x,":"))[2])
  #in.meta$pos = as.numeric(gsub("([0-9]+).*$", "\\1", for_pos))

    in.d.ordered <- in.d[order(in.d[,"chr"], in.d$pos),]
  min.p <- min(in.d[,P.value.column],na.rm=T)
  center <- in.d[which(in.d[,MarkerID.column]==snp),"pos"]
  
  
  region=500000  # Adjust region when necessary
  from <- center - region/2
  to <- center + region/2
  
  gen="hg19"
  #chr="chr11"    # Change when necessary
  #chr_num = 11   # Change when necessary
  
  # snp <- "11-92708710_C-G"
  plot.name = paste0(chr,"_",chartr(":/","--",snp),"_1MB.pdf")    # Change when necessary
  
  
  
  
  ### Step 3: Generate tracks
  ## axTrack
  axTrack <- GenomeAxisTrack(genome = gen , chromosome = chr)
  
  ## idxTrack
  idxTrack <- IdeogramTrack(genome = gen, chromosome = chr)
  
  ## cpgIslands track
  cpgIslands <- UcscTrack(genome = gen, chromosome = chr,
                          track = "cpgIslandExt", from = from, to = to,
                          trackType = "AnnotationTrack", start = "chromStart",
                          end = "chromEnd", id = "name", shape = "box",
                          fill = "#006400", name = "CpG")
  
  ## dtrack
  subset <- in.d.ordered[in.d.ordered$pos %in% c(from: to), ]
  coords <- subset$pos 
  neg_logP <- -log10(subset[,P.value.column])
  dtrack <- DataTrack(data = neg_logP, start = coords,
                      end = coords+1, chromosome = chr, genome = gen,
                      ylim=c(0,ceiling(-log10(min.p))), name = "WGS -logP")
  
  ## ENSEMBL gene track 
  mart <- useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  
  ensembl <- select(mart,keys = chr_num,keytype = "chromosome_name",
                    columns =  c( "ensembl_transcript_id","hgnc_symbol", "chromosome_name","transcript_start", "transcript_end", "ensembl_gene_id"))
  ensembl <- ensembl[order(ensembl$hgnc_symbol),]
  ensembl_transcript <- ensembl[!duplicated(ensembl$hgnc_symbol),][-1,"ensembl_transcript_id"]
  
  biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, start = from, end = to, 
                                      name = "ENSEMBL", filter = list(ensembl_transcript_id = ensembl_transcript))
  
  
  ## Tissue specific gene track
  # Limit to protein-coding genes & Map gene_IDs to transcript_IDs
  gene_adipose_sub <- ensembl[ensembl$ensembl_gene_id %in% adipose_sub[adipose_sub$RPKM>=2 & adipose_sub$Gene_Type=="protein_coding",1],]
  transcript_adipose_sub_chr <- gene_adipose_sub[!duplicated(gene_adipose_sub$ensembl_gene_id),"ensembl_transcript_id"] 
  
  gene_adipose_vis <- ensembl[ensembl$ensembl_gene_id %in% adipose_vis[adipose_vis$RPKM>=2 & adipose_vis$Gene_Type=="protein_coding",1],]
  transcript_adipose_vis_chr <- gene_adipose_vis[!duplicated(gene_adipose_vis$ensembl_gene_id),"ensembl_transcript_id"] # 976 genes
  
  gene_pancreas <- ensembl[ensembl$ensembl_gene_id %in%  pancreas[pancreas$RPKM>=2 & pancreas$Gene_Type=="protein_coding",1],]
  transcript_pancreas_chr <- gene_pancreas[!duplicated(gene_pancreas$ensembl_gene_id),"ensembl_transcript_id"] # 960 genes
  
  gene_muscle <- ensembl[ensembl$ensembl_gene_id %in% muscle[muscle$RPKM>=2 & muscle$Gene_Type=="protein_coding",1],]
  transcript_muscle_chr <- gene_muscle[!duplicated(gene_muscle$ensembl_gene_id),"ensembl_transcript_id"] # 872 genes
  
  gene_liver <- ensembl[ensembl$ensembl_gene_id %in% liver[liver$RPKM>=2 & liver$Gene_Type=="protein_coding",1],]
  transcript_liver_chr <- gene_liver[!duplicated(gene_liver$ensembl_gene_id),"ensembl_transcript_id"]  # 915 genes
  
  
  # Generate tracks 
  biomTrack.adipose_sub <- BiomartGeneRegionTrack(genome = gen, name = "Adi_sub",showId=FALSE, geneSymbol=FALSE, filters = list(ensembl_transcript_id=transcript_adipose_sub_chr))
  biomTrack.adipose_vis <- BiomartGeneRegionTrack(genome = gen, name = "Adi_vis",showId=FALSE, geneSymbol=FALSE, filters = list(ensembl_transcript_id=transcript_adipose_vis_chr))
  biomTrack.pancreas <- BiomartGeneRegionTrack(genome = gen, name = "Pancreas", showId=FALSE, geneSymbol=FALSE, filters = list(ensembl_transcript_id=transcript_pancreas_chr))
  biomTrack.muscle <- BiomartGeneRegionTrack(genome = gen, name = "Muscle", showId=FALSE, geneSymbol=FALSE, filters = list(ensembl_transcript_id=transcript_muscle_chr))
  biomTrack.liver <- BiomartGeneRegionTrack(genome = gen, name = "Liver",showId=FALSE, geneSymbol=FALSE, filters = list(ensembl_transcript_id=transcript_liver_chr))
  
  
  
  ### Step 4: Generate plots
  
  if(plot) {
    pdf(plot.name,height=12) 
  }
  plotTracks(list(axTrack,idxTrack,dtrack,cpgIslands,biomTrack,
                  biomTrack.adipose_sub,biomTrack.adipose_vis,biomTrack.pancreas,biomTrack.muscle,biomTrack.liver),
             from = from, to = to,transcriptAnnotation = "symbol",  showTitle = TRUE)
  if(plot) {
    dev.off()
  }
}
