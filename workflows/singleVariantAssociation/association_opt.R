
# GENESIS-test-v3.R

# INPUT:
# (1) a genotype (gds) file
# (2) a phenotype (ped) file
# (3) a GRM kinship.objects.RData file (from kinshipMatrix.R file)
# (4) a commonIDs .RData file
# (5) the name of the column containing the IDs in the ped file
# (6) a label describing the workflow (for naming the assoc results),
# (7) the outcome under study
# (8) the outcome type (dichot/cont), 
# (9) covariates to use, separated by "," and no space. Example: "age,sexTrait,BMI"   

# OUTPUT:
# An assoc.RData file

# command line example
# R --vanilla --args \
#  GoT2D.chr22.biallelic.gds \
#  GoT2D.phenotype.ped \
#  GoT2D_T2D.kinship.objects.RData \ 
#  GoT2D_T2D.commonIDs.RData \
#  GoT2DsampleID \
#  GoT2D_T2D \
#  T2D \ 
#  dichotomous \
#  age sex BMI 

library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)

input_args <- commandArgs(trailingOnly=T)

gds <- input_args[1] #"GoT2D.chr22.biallelic.gds"
ped <- input_args[2] #"GoT2D.phenotype.ped"
kinship_file <- input_args[3] # GoT2DSampleID
commonID_file <- input_args[4]
id.column.name <- input_args[5]
label <- input_args[6]
outcome <- input_args[7]
outcomeType <- input_args[8]
minMAC <- 30 # hard coded
covariates <- NULL

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
ped_file <- read.table(ped, header = TRUE, as.is = FALSE)
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



 

