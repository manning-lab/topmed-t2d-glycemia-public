## make aggregation units

getAggUnits <- function( state.file , peak.file , anno.file , gene.file , geno.gds , state.names , anno.type , anno.value ){
  
  ##### these functions are from the UWGAC topmed analysis pipeline #####
  #credit: https://github.com/UW-GAC/analysis_pipeline/blob/master/TopmedPipeline/R/filterVariants.R
  # use fruntion from wash u for subetting
  aggregateListByPosition <- function(gds, groups, indexOnly=FALSE) {
    stopifnot(all(c("group_id", "chromosome", "start", "end") %in% names(groups)))
    
    ## select only variants in requested regions
    filtOrig <- seqGetFilter(gds)
    gr <- GRanges(seqnames=groups$chromosome,
                  ranges=IRanges(groups$start, groups$end, names=groups$group_id))
    seqSetFilter(gds, gr, verbose=FALSE)
    
    variants <- .expandAlleles(gds)
    
    seqSetFilter(gds, sample.sel=filtOrig$sample.sel,
                 variant.sel=filtOrig$variant.sel, verbose=FALSE)
    
    ## find group_id for each variant
    vr <- GRanges(seqnames=variants$chromosome,
                  ranges=IRanges(variants$position, variants$position, names=variants$variant.id))
    ol <- findOverlaps(vr, gr)
    map <- data.frame(group_id=names(gr)[subjectHits(ol)],
                      variant.id=as.integer(names(vr))[queryHits(ol)],
                      stringsAsFactors=FALSE)
    
    variants <- distinct_(map) %>%
      left_join(variants, by="variant.id")
    .groupVariants(variants, indexOnly)
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
  .groupVariants <- function(variants, indexOnly) {
    
    ## columns to return
    extraCols <- if (indexOnly) character(0) else c("chromosome", "position", "ref", "nAlleles", "allele")
    
    groups <- unique(variants$group_id)
    lapply(setNames(groups, groups), function(g) {
      filter_(variants, ~(group_id == g)) %>%
        select_(~(one_of(c("variant.id", extraCols, "allele.index"))))
    })
  }
  
  state.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/data/varshney_islets_chromatin_state.aggregation.chr10.csv"
  peak.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/data/varshney_islets_abcu196_50_broad_peaks.aggregation.chr10.csv"
  anno.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/data/freezes_2a_3a_4.chr10.snp.annotated.general20170422.gz"
  gene.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/data/gtex.v6p.pancreas.expression.min.rpkm.0.1.txt"
  geno.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/data/freeze4.chr10.pass.gtonly.minDP10.genotypes.gds"
  state.names <- c("active_tss","active_enhancer_1")
  anno.value <- c("splice_acceptor_variant","splice_donor_variant","stop_gained","splice_region_variant")
  minmaf <- 0.01
  # anno.value <- c("stopgain", "splicing")
  
  # data.table to read fast
  library(data.table)
  library(SeqArray)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(SeqVarTools)
  

  geno.gds <- seqOpen(geno.file)
  
  ref.freq <- seqAlleleFreq(geno.gds, .progress=TRUE, parallel = TRUE)
  maf <- pmin(ref.freq, 1-ref.freq)
  vid.id <- seqGetData(geno.gds,"variant.id")
  vid.filt <- vid.id[which(maf <= minmaf)]
  seqSetFilter(geno.gds, variant.id=vid.filt, action="intersect", verbose=TRUE)
  
  
  chr_num <- seqGetData(geno.gds,"chromosome")[1]
  # subset to just the desired states # load chromostate file
  state.raw <- fread(state.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  state.raw <- state.raw[which(state.raw[,7] %in% state.names),]
  state.raw[,4] <- sub("chr","",state.raw[,4])
  
  # load peaks
  peak.raw <- fread(peak.file,sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  peak.raw[,4] <- sub("chr","",peak.raw[,4])
  
  # load annotation file
  # library(devtools)
  # devtools::install_github("UW-GAC/wgsaparsr@1.0.0.9003")
  # library(wgsaparsr)
  # desired_columns <- c("`#chr`","pos","VEP_ensembl_Gene_ID","VEP_ensembl_Consequence")
  # to_split <- c("VEP_ensembl_Gene_ID","VEP_ensembl_Consequence")
  # parse_to_file(anno.file, "anno_subset.tsv", desired_columns, to_split, verbose = TRUE)
  anno.raw <- fread("anno_subset.tsv",sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  
  # load gene file
  library(biomaRt)
  mart <- useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  # ensembl <- select(mart,keys = chr_num,keytype = "chromosome_name",
  #                   columns =  c( "ensembl_transcript_id","hgnc_symbol", "chromosome_name","transcript_start", "transcript_end", "ensembl_gene_id"))
  ensembl <- biomaRt::select(mart,keys = chr_num,keytype = "chromosome_name",
                    columns =  c("ensembl_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position"))
  ensembl <- ensembl[order(ensembl$hgnc_symbol),]
  ensembl <- ensembl[!duplicated(ensembl$hgnc_symbol),][-1,]
  gene.raw <- fread(gene.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  gene.prep <- ensembl[ensembl$ensembl_gene_id %in%  gene.raw[gene.raw$RPKM>=2 & gene.raw$Gene_Type=="protein_coding",1],]
  gene.prep <- gene.prep[!duplicated(gene.prep$ensembl_gene_id),]
  
  # get TADs
  ensembl <- select(mart,keys = chr_num,keytype = "chromosome_name",
                    columns =  c("ensembl_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position"))
  
  
  #granges
  peak.gr <- GRanges(seqnames=peak.raw$chr,ranges=IRanges(start=peak.raw$start,end=peak.raw$stop))
  states.gr <- GRanges(seqnames=state.raw$chr,ranges=IRanges(start=state.raw$start,end=state.raw$stop))
  gene.gr <- GRanges(seqnames=gene.prep$chromosome_name, ranges=IRanges(start=gene.prep$start_position-5000,end=gene.prep$end_position+5000),gene_id=gene.prep$ensembl_gene_id)
  gene.gr <- reduce(gene.gr)
  
  seqSetFilter(geno.gds,gene.gr)
  
  # filter snps for gene regions
  #seqSetFilter(geno.gds, gene.gr, verbose=FALSE)
  
  # define our groups by genes first
  # groups <- data.frame(group_id=seq(length(gene.prep[,1])), chromosome=gene.prep[,3], start=gene.prep[,4], end=gene.prep[,5])
  groups <- data.frame(group_id=seq(length(gene.gr)), chromosome=seqnames(gene.gr), start=start(gene.gr), end=end(gene.gr))
  
  # get list of variants with at least one of our tags
  
  anno.sub <- anno.raw[which(anno.raw$VEP_ensembl_Gene_ID != "."),]
  
  ininds <- c()
  for (i in seq(1,length(anno.sub[,1]))){
    if (any(str_detect(anno.value,anno.sub$VEP_ensembl_Consequence[i]))) {
      ininds <- c(ininds,i)
    }
  }
  
  # anno.sub <- anno.sub[apply(sapply(anno.raw$VEP_ensembl_Consequence, function(x) str_detect(x,anno.value)), 2, any),]
  anno.sub <- anno.sub[ininds,]
  
  # anno.sub <- anno.raw$pos[apply(sapply(anno.raw$VEP_ensembl_Consequence, function(x) str_detect(x,anno.value)), 2, any)]
  snp.prot <- data.frame(names=seqGetData(geno.gds,"variant.id"), pos=seqGetData(geno.gds,"position"))
  snp.prot <- snp.prot[which(snp.prot$pos %in% anno.sub$pos),]$names
  filtOrig <- seqGetFilter(geno.gds)
  filt_pregene <- filtOrig
  seqSetFilter(geno.gds,gene.gr) # 2261861 variants
  
  
  # get snps in open peaks
  filtOrig <- seqGetFilter(geno.gds)
  # seqSetFilter(geno.gds, variant.sel=peak.gr,action="intersect")
  seqSetFilter(geno.gds, peak.gr)
  
  snp.peaks <- seqGetData(geno.gds,"variant.id")
  seqResetFilter(geno.gds)
  seqSetFilter(geno.gds, sample.sel=filtOrig$sample.sel,variant.sel=filtOrig$variant.sel)
  
  # get snps in islet regions (chromatin states)
  seqSetFilter(geno.gds, states.gr)
  snp.states <- seqGetData(geno.gds,"variant.id")
  seqResetFilter(geno.gds)
  seqSetFilter(geno.gds, sample.sel=filtOrig$sample.sel,variant.sel=filtOrig$variant.sel)
  
  # get our variants lists for each mask
  
  # regulatory or stopgain/splice
  seqSetFilter(geno.gds, variant.id=c(snp.prot,snp.states),action="intersect")
  variants.reg_ptv <- seqGetData(geno.gds,"variant.id")
  print(length(variants.reg_ptv)) # 115762
  mask.reg_ptv <- aggregateListByPosition(geno.gds, groups)
  seqResetFilter(geno.gds)
  seqSetFilter(geno.gds, sample.sel=filtOrig$sample.sel,variant.sel=filtOrig$variant.sel)
  
  # regulatory or broad peak or stopgain/splice
  seqSetFilter(geno.gds, variant.id=c(snp.prot,snp.states, snp.peaks),action="intersect")
  variants.pks_reg_ptv <- seqGetData(geno.gds,"variant.id")
  print(length(variants.pks_reg_ptv)) # 261110
  mask.pks_reg_ptv <- aggregateListByPosition(geno.gds, groups)
  seqResetFilter(geno.gds)
  seqSetFilter(geno.gds, sample.sel=filtOrig$sample.sel,variant.sel=filtOrig$variant.sel)
  
  # write mask files as r data objects
  mask.file.split <- unlist(strsplit(geno.file,"\\."))
  mask.file.base <- paste(paste(mask.file.split,collapse = "."),"minmaf",minmaf,sep=".")
  save(mask.reg_ptv, "mask.reg_ptv", file=paste(mask.file.base,paste(state.names,collapse="."),paste(anno.value,collapse="."),"groups","RData",sep="."))
  save(mask.pks_reg_ptv, "mask.pks_reg_ptv", file=paste(mask.file.base,paste(state.names,collapse="."),paste(anno.value,collapse="."),"peaks","groups","RData",sep="."))
  write(variants.reg_ptv, file=paste(mask.file.base,paste(state.names,collapse="."),paste(anno.value,collapse="."),"variantlist","txt",sep="."),ncolumns=1)
  write(variants.pks_reg_ptv, file=paste(mask.file.base,paste(state.names,collapse="."),paste(anno.value,collapse="."),"peaks","variantlist","txt",sep="."),ncolumns=1)
  
  # close gds file
  #seqClose(geno.gds)
  
  return(mask)
}

plotAggUnit <- function( geno.gds , groups , pancrease ){
  
  library(Gviz)
  library(biomaRt)
  ##
  groups <- mask.reg_ptv
  gene_pancreas <- gene.prep
  pancreas <- gene.raw
  ##
  
  cur_group <- as.matrix(groups[1])[1][[1]]
  buff <- 20000
  from <- min(cur_group$position) - buff
  to <- max(cur_group$position) + buff
  chr <- cur_group$chromosome[1]
  gen <- "hg19"
  
  plot.name = "test.pdf"
  axTrack <- GenomeAxisTrack(genome = gen , chromosome = chr)
  
  ## idxTrack
  idxTrack <- IdeogramTrack(genome = gen, chromosome = chr)
  
  ## cpgIslands track
  # cpgIslands <- UcscTrack(genome = gen, chromosome = chr,
  #                         track = "cpgIslandExt", from = from, to = to,
  #                         trackType = "AnnotationTrack", start = "chromStart",
  #                         end = "chromEnd", id = "name", shape = "box",
  #                         fill = "#006400", name = "CpG")
  
  ## dtrack
  coords <- cur_group$position
  filtOrig <- seqGetFilter(geno.gds)
  # seqResetFilter(geno.gds)
  v.gr <- GRanges(seqnames=c(chr),ranges=IRanges(start=min(cur_group$position),en=max(cur_group$position)))
  seqSetFilter(geno.gds,v.gr)
  seqSetFilter(geno.gds,variant.id=cur_group$variant.id, action="intersect")
  ref.freq <- seqAlleleFreq(geno.gds)
  maf <- pmin(ref.freq, 1-ref.freq)
  vid.id <- seqGetData(geno.gds,"variant.id")
  wm <- which(maf <= minmaf)
  wm <- unique(wm)
  maf <- maf[wm]
  vid.filt <- vid.id[wm]
  seqSetFilter(geno.gds, variant.id=vid.filt, action="intersect", verbose=TRUE)
  # ref.freq <- seqAlleleFreq(geno.gds)
  # maf <- pmin(ref.freq, 1-ref.freq)
  coords <- seqGetData(geno.gds,"position")
  lm <- -log10(maf)
  keep <- !is.infinite(lm)
  lm <- lm[keep]
  coords <- coords[keep]
  # maf <- seqAlleleFreq(geno.gds)
  seqResetFilter(geno.gds)
  seqSetFilter(geno.gds, sample.sel=filtOrig$sample.sel,variant.sel=filtOrig$variant.sel)
  
  dtrack <- DataTrack(data = lm, start = coords,
                      end = coords+1, chromosome = chr, genome = gen,
                      ylim=c(floor(min(lm)),ceiling(max(lm))), name = "MAF")
  
  
  ######################## here ############################
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
  
  # gene_pancreas_chr <- gene_pancreas[!duplicated(gene_pancreas$ensembl_gene_id),"ensembl_gene_id"] # 960 genes
  # biomTrack.pancreas <- BiomartGeneRegionTrack(genome = gen, name = "Pancreas", showId=TRUE, geneSymbol=FALSE, filters = list(ensembl_gene_id=gene_pancreas_chr))
  # 
  gene_pancreas <- ensembl[ensembl$ensembl_gene_id %in%  pancreas[pancreas$RPKM>=2 & pancreas$Gene_Type=="protein_coding",1],]
  transcript_pancreas_chr <- gene_pancreas[!duplicated(gene_pancreas$ensembl_gene_id),"ensembl_transcript_id"] # 960 genes
  biomTrack.pancreas <- BiomartGeneRegionTrack(genome = gen, name = "Pancreas", showId=FALSE, geneSymbol=FALSE, filters = list(ensembl_transcript_id=transcript_pancreas_chr))
  
  pdf(plot.name) 
  
  # plotTracks(list(axTrack,idxTrack,dtrack,cpgIslands,biomTrack, biomTrack.pancreas),from = from, to = to,transcriptAnnotation = "symbol",  showTitle = TRUE)
  # plotTracks(list(axTrack,idxTrack,dtrack,biomTrack, biomTrack.pancreas),from = from, to = to,transcriptAnnotation = "symbol",  showTitle = TRUE)
  plotTracks(list(axTrack,idxTrack,dtrack),from = from, to = to,transcriptAnnotation = "symbol",  showTitle = TRUE)
  dev.off()
}
  