## make aggregation units

getAggUnits <- function( state.file , anno.file , geno.gds , state.names , anno.type , anno.value ){

  state.file <- "/Users/tmajaria/Documents/projects/topmed/code/varshney/data/local_cs_states/Islets.chromatinStates.bed"
  anno.file <- ""
  geno.file <- "/Users/tmajaria/Documents/projects/topmed/code/testing_inputs/singleVariantFull/freeze4.chr21.pass.gtonly.minDP10.genotypes.gds"
  state.names <- c("1_Active_TSS")
  
  # data.table to read fast
  library(data.table)

  # load chromostate file
  state.raw <- fread(state.file,sep="\t",header=F,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  
  # subset to just the desired states
  state.raw <- state.raw[which(state.raw[,4] %in% state.names),]
  
  # correct chromosome col
  state.raw[,1] <- sub("chr","",state.raw[,1])
  
  # seqArray for gds
  library(SeqArray)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(SeqVarTools)
  

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
  groups <- data.frame(group_id=seq(length(state.raw[,1])), chromosome=state.raw[,1], start=state.raw[,2], end=state.raw[,3])
  variants <- aggregateListByPosition(geno.gds, groups)
  
  # close gds file
  #seqClose(geno.gds)
  
  return(variants)
}
  