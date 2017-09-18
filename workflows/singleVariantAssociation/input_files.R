

#P.value.column <- get.pvalue.column(study)
#MarkerID.column <- get.markerID.column(study)
get.pvalue.column <- function(study) {
  if(study=="meta") return("P.value")
  if(study=="GENESIS_T2D") return("Score.pval") 
  if(study == "GENESIS-FIRECLOUD") return("P")
  if(study == "GENESIS-QT-WALD") return("Wald.pval")
}

get.markerID.column <- function(study) {
  if(study=="meta") return("MarkerName")
  if(study=="GENESIS_T2D") return("snpID")
  if(study == "GENESIS-FIRECLOUD" | study == "GENESIS-QT-WALD") return("snpID")
}

get.input <- function(in.file.name, study) {

    if(study=="meta") {
    in.file <- read.table(in.file.name,as.is=T,header=T,sep="\t")
    tmp <- sapply(in.file$MarkerName,function(x){strsplit(x,split=":|_")[[1]]})
    in.file$chr <- as.numeric(ifelse(tmp[1,]!="X",tmp[1,],23))
    
    in.file$g.chr <- paste("chr",tmp[1,],sep="")
    
    in.file$pos <- as.numeric(tmp[2,])
    in.file$log10P <- -log10(in.file$P.value)
    in.file$MAF <- ifelse(in.file$Freq1<0.5,in.file$Freq1,1-in.file$Freq1)
    return(in.file)
    
  } else if (study == "GENESIS_T2D") {
    in.file <- c()
    for(in.file.name in list.files(pattern="^T2D_Adj_chr")) {
      tmp <- read.csv(in.file.name,header=T,as.is=T)
      in.file <- rbind(in.file,tmp)
      print(dim(in.file))
    }    
    
    in.file$g.chr <- paste("chr",in.file$chr,sep="")
    return(in.file)
  } else if (study == "GENESIS-FIRECLOUD" | study == "GENESIS-QT-WALD") {
    in.file <- read.csv(in.file.name,header=T,as.is=T)
    in.file$g.chr <- paste("chr",in.file$chr,sep="")
    
    return(in.file)
  }

  
  
}


