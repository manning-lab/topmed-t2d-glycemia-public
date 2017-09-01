fileConn<-file("output.txt")
writeLines(c("Hello","World"), fileConn)
close(fileConn)
