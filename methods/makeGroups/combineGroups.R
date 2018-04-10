args <- commandArgs(trailingOnly=T)
groups <- unlist(strsplit(args[1],","))
outpref <- args[2]

library(data.table)

all_groups <- list()
for (i in seq(1,length(groups))){
	all_groups[[i]] <- fread(groups[i], data.table = F)
}

all_groups <- do.call(rbind,all_groups)

fwrite(all_groups, paste(outpref,"all.groups.mask3.csv", sep = "."), sep = ",")