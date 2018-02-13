# LOADING SOURCE FILES FOR GLY PHENO HARMONIZATION
# NOTE : you need to organize all of these files into a single directory. The input is then the absolutely path (ex: /User/username/Documents/.../gly_phenotype_files/) to that folder
get_pheno_data <- function(f.dir){
	
	# define file paths
	linker.file <- paste(f.dir, "freeze5b_sample_annot_2017-10-16.txt", sep="/")
	fhs.file <- paste(f.dir, "FHS_glycemicTraits_20170214_CS.ped", sep="/")
	fhs.pedigree <- paste(f.dir, "share_unrel_comb_01012017.csv", sep="/")
	jhs.file <- paste(f.dir, "JHS_glycemictraits_29Mar2017.ped", sep="/")
	sas.file <- paste(f.dir, "SAS_LIULIN_20170212_glycemic_traits.ped", sep="/")
	cfs.file <- paste(f.dir, "CFS_12May17_glycemic_traits_Freeze4.ped", sep="/")
	amish.file <- paste(f.dir, "Amish_HuichunXu_May15th2017_glycemic_traits.ped", sep="/")
	gensalt.file <- paste(f.dir, "GenSalt_EA_XuenanMi_20170614_glycemic_traits.ped", sep="/")
	aric.EU.file <- paste(f.dir, "ARIC_EU_NDAYA_20161220_glycemic_traits_ped.csv", sep="/")
	aric.AF.file <- paste(f.dir, "ARIC_AA_NDAYA_20161220_glycemic_traits_ped.csv", sep="/")
	whi.HS.file <- paste(f.dir, "WHI_HA_MP_20170726_glycemic_traits.ped", sep="/")
	whi.EU.file <- paste(f.dir, "WHI_EU_MP_20170726_glycemic_traits.ped", sep="/")
	whi.AS.file <- paste(f.dir, "WHI_AS_MP_20170726_glycemic_traits.ped", sep="/")
	whi.AF.file <- paste(f.dir, "WHI_AA_MP_20170726_glycemic_traits.ped", sep="/")
	genestar.AF.file <- paste(f.dir, "GeneSTAR_AA_YANEK_20170720_glycemic_traits.ped", sep="/")
	genestar.EU.file <- paste(f.dir, "GeneSTAR_EU_YANEK_20170720_glycemic_traits.ped", sep="/")
	chs.file <- paste(f.dir, "CHS_FLOYD_201708018_glycemic_traits.PED", sep="/")
	mesa.HS.file <- paste(f.dir, "MESA_HA_ABIGAILBALDRIDGE_04JAN17_glycemic_traits_sidno.csv", sep="/")
	mesa.EU.file <- paste(f.dir, "freeze5b_sample_annot_2017-10-16.txt", sep="/")
	mesa.AS.file <- paste(f.dir, "MESA_SA_ABIGAILBALDRIDGE_04JAN17_glycemic_traits_sidno.csv", sep="/")
	mesa.AF.file <- paste(f.dir, "MESA_AA_ABIGAILBALDRIDGE_04JAN17_glycemic_traits_sidno.csv", sep="/")
	genoa.file <- paste(f.dir, "GENOA_AA_BIELAK_20171124_glycemic_traits.ped", sep="/")

	# load all the data
	linker <- read.table(linker.file)

	fhs <- read.table(fhs.file,header=T)
	share<-read.csv(fhs.pedigree,header=T)
	jhs <- read.table(jhs.file,na.strings="x",header=T)
	sas<-read.table(sas.file,header=T)
	cfs<-read.table(cfs.file,header=T)
	amish<-read.table(amish.file,header=T,na.strings=" ",sep="\t")
	gensalt<-read.table(gensalt.file,head=T)
	aric_ea<-read.table(aric.EU.file,header=T,sep=',')
	aric_aa<-read.csv(aric.AF.file,header=T
	whi_ha<-read.table(whi.HS.file,header=F,fill=T,na.strings=" ",sep="\t")
	whi_ea<-read.table(whi.EU.file,header=F,fill=T,na.strings=" ",sep="\t")
	whi_as<-read.table(whi.AS.file,header=F,fill=T,na.strings=" ",sep="\t")
	whi_aa<-read.table(whi.AF.file,header=F,fill=T,na.strings=" ",sep="\t")
	gs_aa<-read.table(genestar.AF.file,header=T)
	gs_ea<-read.table(genestar.EU.file,header=T)
	chs<-read.table(chs.file,header=T,fill=T,sep='\t')
	mesa_ha<-read.csv(mesa.HS.file,header=T)
	mesa_ea<-read.csv(mesa.EU.file,header=T)
	mesa_sa<-read.csv(mesa.AS.file,header=T)
	mesa_aa<-read.csv(mesa.AF.file,header=T)
	genoa<-read.table(genoa.file,header=T,fill=T,sep='\t',na.strings=" ")
	sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T)
}