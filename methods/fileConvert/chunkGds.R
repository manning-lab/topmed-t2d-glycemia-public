# break full genome gds files into chromosome->chromosome chunk files

library(SeqArray)
input_args <- commandArgs(trailingOnly=T)
gds.file <- inputs_args[1]
outbase <- basename(gds.file)
outbase <- substr(outbase, 0, nchar(outbase)-4)
outbase <- paste(getwd(),"/",outbase,sep="")

# https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
chunk <- function(x,n) {
  f <- sort(rep(1:(trunc(length(x)/n)+1),n))[1:length(x)]
  return(split(x,f))
}

gds <- seqOpen(gds.file)
variant.id <- seqGetData(gds, "variant.id")
variant.id.chunks <- chunk(variant.id,100)

# loop through the chunks
for(j in 1:length(variant.id.chunks)) {
  print(paste("Chunk ",j," / ",length(variant.id.chunks),sep=""))
  seqSetFilter(gds, variant.id=variant.id.chunks[[j]])
  seqExport(gds, paste(outbase,".chunk",j,".gds",sep=""))
  # reset filter
  seqResetFilter(gds)
}
seqClose(gds)