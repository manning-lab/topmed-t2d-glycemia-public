input_args <- commandArgs(trailingOnly=T)
fname <- input_args[1] 
d <- readChar(fname, file.info(fname)$size)
# fileConn<-file("output.txt")
# writeLines(c("Hello","World"), fileConn)
# close(fileConn)
lapply(d, write, "output.txt", append=TRUE)
