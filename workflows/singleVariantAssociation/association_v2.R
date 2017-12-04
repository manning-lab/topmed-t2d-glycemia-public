# GENESIS-test-v3.R .. accepts pre-configured null-model as input

# INPUT:
# (1) a genotype (gds) file
# (2) a label describing the workflow (for naming the assoc results),
# (3) the outcome type (dichot/cont),  
# (4) the null model (.RData file)

# OUTPUT:
# An assoc.RData file

# command line example
# R --vanilla --args \
#  GoT2D.chr22.biallelic.gds \
#  dichotomous \
#  GoT2D.nullModel.RData

library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
label <- input_args[6]
outcomeType <- input_args[8]
nullModel <- input_args[4]
minMAC <- 30 # hard coded

print(paste("GDS file:",gds))
print(paste("Ped file:",ped))
print(paste("Kinship file:",kinship_file))
print(paste("Common IDs file:",commonID_file))
print(paste("ID column name:",id.column.name))
print(paste("label:",label))
print(paste("outcome:",outcome))
print(paste("outcomeType:",outcomeType))


if(length(input_args)>8) {
	
	covariates <- strsplit(input_args[9],split=",")[[1]]
	print(paste("covariates:",paste(covariates,collapse=" ")))

}



# load ped_file and re-order columns
# ped_file <- read.table(ped, header = TRUE, as.is = FALSE)
ped_file <- fread(ped,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
ped_file <- ped_file[!duplicated(ped_file[,id.column.name]),]
rownames(ped_file) <- ped_file[,id.column.name]
head(ped_file)

# load common ID and kinship files
load(commonID_file)

library(gdsfmt)
geno.pcrelate <- openfn.gds(kinship_file)
# load(paste(kinship_file)) #ped_file.trimmed, PCRgrm, filteredSNPs, both
  
## trim phenotype file 
ped_file.trimmed <- ped_file[ped_file[,id.column.name] %in% both, ]
ped_file.trimmed <- ped_file.trimmed[both,]


## load gds file using GWASTools/GENESIS
#geno <- GdsGenotypeReader(filename=gds)
#genoData <- GenotypeData(geno)

# convert phenotype dataframe to ScanAnnotationDataFrame object
if(!is.null(covariates)) {
	df <- cbind(setNames(data.frame(both,ped_file.trimmed[,outcome]), c("scanID",outcome)), covs = ped_file.trimmed[,covariates])

	# change covariates = [age, sex, BMI] to covariates = [covs.age, covs.sex, covs.BMI]
	covariates <- paste("covs.", covariates, sep="")

	# double check sex column;
	if("covs.sex" %in% covariates) {
		print("Sex is one of the covariates")
		if(sum(c(1,2) %in% names(table(df[,"covs.sex"])))==2) {
			df$sex <- ifelse(df[,"covs.sex"]==2,"F","M")
			print("Converted sex from 1/2 to M/F")
			print(table(sex=df$sex,covs.sex=df$covs.sex,useNA="always"))
		} else {
			df$sex <- df$covs.sex
			print(table(sex=df$sex,covs.sex=df$covs.sex,useNA="always"))
		}
		covariates <- c(setdiff(covariates,"covs.sex"),"sex")
		print(paste("covariates:",paste(covariates,collapse=" ")))	
	}
} else {
	df <- setNames(data.frame(both,ped_file.trimmed[,outcome]), c("scanID",outcome))

}

scanAnnot <- ScanAnnotationDataFrame(df)
str(scanAnnot)

#### Null model - Future Versions of this script can be extended to input null models

# get GRM in correct format
PCRgrm <- pcrelateMakeGRM(geno.pcrelate, scan.include = both, scaleKin = 2)

if(outcomeType=="dichotomous" & is.null(covariates)) {
	        nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = binomial, covMatList = PCRgrm)
}
if(outcomeType=="dichotomous" & !is.null(covariates)) {
	        nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = binomial, covMatList = PCRgrm, covars = covariates)
}
if(outcomeType=="continuous" & !is.null(covariates)) {
	        nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = gaussian, covMatList = PCRgrm, covars = covariates)
}

if(outcomeType=="continuous" & is.null(covariates)) {
	        nullmod <- fitNullMM(scanData = scanAnnot, outcome = outcome, family = gaussian, covMatList = PCRgrm)
}


#### run association test
geno <- seqOpen(gds)

# filter by MAF
seqSetFilter(geno,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)

ref.freq <- seqAlleleFreq(geno, .progress=TRUE)
maf <- pmin(ref.freq, 1-ref.freq)
maf.filt <- 2 * maf * (1-maf) * length(nullmod$scanID) >= minMAC
print(table(maf.filt))

if(sum(maf.filt)==0) {
print("No SNPs pass MAC filter. Finished Association Step")
assoc <- NA

} else {
seqSetFilter(geno, variant.sel=maf.filt, action="intersect", verbose=TRUE)

## add position and rsID
pos <- seqGetData(geno,"position")
#chr.pos <- seqGetData(geno, "$chrom_pos") This line didn't seem to work after filtering geno
allele <- seqGetData(geno, "allele")
snps.pos <- cbind(pos,allele)
print("Filtered SNPs")
print(dim(snps.pos))

genoData <- SeqVarData(geno)

## dichotomous 
if(outcomeType=="dichotomous" ) {
	assoc <- assocTestMM(genoData = genoData, nullMMobj = nullmod, test = "Score")
}

## continuous 
if(outcomeType=="continuous" ) {
	assoc <- assocTestMM(genoData = genoData, nullMMobj = nullmod, test = "Wald")
}
print("Finished Association Step")
print(dim(assoc))
assoc <- cbind(snps.pos, assoc)

}
## save assoc object
save(assoc, file=paste(label, ".assoc.RData", sep=""))
