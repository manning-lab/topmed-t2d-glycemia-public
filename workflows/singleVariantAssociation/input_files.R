

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

get.input <- function(in.file.name, study,in.file.path=NA,min.maf=NA) {

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
    library(readr)
    in.file <- c()
    print(getwd())
    print(list.files(path=in.file.path,pattern=in.file.name))
    for(in.file.name in list.files(path=in.file.path,pattern=in.file.name)) {
      print(in.file.name)
      tmp <- read_csv(paste(in.file.path,"/",in.file.name,sep=""),col_names=T,
                      col_types=list(
        snpID = col_character(),
        chr = col_integer(),
        n = col_integer(),
        MAF = col_double(),
        minor.allele = col_character(),
        Score = col_double(),
        Var = col_double(),
        Score.Stat = col_double(),
        Score.pval = col_double(),
        pos = col_integer(),
        ref = col_character(),
        alt = col_character()
      ))
      in.file <- rbind(in.file,tmp[which(tmp[,"MAF"]>min.maf),])
      tail(in.file)
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


