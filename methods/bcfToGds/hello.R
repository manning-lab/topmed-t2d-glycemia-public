input_args <- commandArgs(trailingOnly=T)
fname <- input_args[1] 
d = read.table(fname)
# fileConn<-file("output.txt")
# writeLines(c("Hello","World"), fileConn)
# close(fileConn)
write(d, file = "output.txt")
