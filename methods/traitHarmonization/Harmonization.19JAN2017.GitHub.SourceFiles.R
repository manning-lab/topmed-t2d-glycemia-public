# LOADING SOURCE FILES FOR T2D PHENO HARMONIZATION
# NOTE : you need to organize all of these files into a single directory. The input is then the absolutely path (ex: /User/username/Documents/.../t2d_phenotype_files/) to that folder

args <- commandArgs(trailingOnly=T)
f.dir <- args[1]

## T2D Harmonization  Source Files
# list of files:
sample.file <- paste(f.dir, "freeze5b_sample_annot_2017-12-01.txt", sep="/")
MGH.file <- paste(f.dir, "MGH_Parent_dbGaP_SubjectPhenotypesDS_v3.txt", sep="/")
MGH.pedigree <- paste(f.dir, "MGH_Parent_PedigreeDS_v1.txt", sep="/")
VU.file <- paste(f.dir, "phs001032.v3.pht005675.v3.p2.c1.TOPMed_WGS_VUH_AF_Subject_Phenotypes.GRU-IRB.txt", sep="/")
VU.pedigree <- paste(f.dir, "phs001032.v3.pht007135.v1.p2.TOPMed_WGS_VUH_AF_Pedigree.MULTI.txt", sep="/")
afp.file <- paste(f.dir, "phs001024.v2.pht005693.v1.p1.c1.TOPMed_WGS_Partners_AFGen_Subject_Phenotypes.HMB.txt", sep="/")
afvub.file <- paste(f.dir, "VUH_Ablation_dbGaP_SubjectPhenotypesDS_v2.txt'", sep="/")
afccaf.file <- paste(f.dir, "CCAF_dbGaP_SubjectPhenotypesDS_v5.txt", sep="/")
HVH.file <- paste(f.dir, "HVH_FLOYD_20170818_T2DPhenotype.PED", sep="/")
ARIC.AF.file <- paste(f.dir, "ARIC_AA_NMAYA_20161220_T2D.ped", sep="/")
ARIC.EU.file <- paste(f.dir, "ARIC_EU_NMAYA_20161220_T2D.ped", sep="/")
copd.1.file <- paste(f.dir, "phs000179.v5.pht002239.v4.p2.c1.COPDGene_Subject_Phenotypes.HMB.txt", sep="/")
copd.2.file <- paste(f.dir, "phs000179.v5.pht002239.v4.p2.c2.COPDGene_Subject_Phenotypes.DS", sep="/")
amish.file <- paste(f.dir, "Amish_HuichunXu_May152017_T2D.ped", sep="/")
cfs.file <- paste(f.dir, "CFS_12May17_T2D_Freeze4.ped", sep="/")
sas.file <- paste(f.dir, "SAS_LIULIN_20170212_T2D.ped", sep="/")
jhs.file <- paste(f.dir, "JHS_T2Dincidence_16Aug2017.ped", sep="/")
fhs.file <- paste(f.dir, "FHS_T2D_20170214_CS.ped", sep="/")
whi.EU.file <- paste(f.dir, "WHI_EU_MP_20170726_T2D.ped", sep="/")
whi.AF.file <- paste(f.dir, "WHI_AA_MP_20170726_T2D.ped", sep="/")
whi.AS.file <- paste(f.dir, "WHI_AS_MP_20170726_T2D.ped", sep="/")
whi.HS.file <- paste(f.dir, "WHI_HA_MP_20170726_T2D.ped", sep="/")
gensalt.EU.file <- paste(f.dir, "GenSalt_EA_XuenanMi_20170614_T2D.ped", sep="/")
genestar.AF.file <- paste(f.dir, "GeneSTAR_AA_YANEK_20170720_T2D.ped", sep="/")
genestar.EU.file <- paste(f.dir, "GeneSTAR_EU_YANEK_20170720_T2D.ped", sep="/")
mesa.AF.file <- paste(f.dir, "MESA_AA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv", sep="/")
mesa.EU.file <- paste(f.dir, "MESA_EU_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv", sep="/")
mesa.HS.file <- paste(f.dir, "MESA_HA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv", sep="/")
mesa.AS.file <- paste(f.dir, "MESA_SA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv", sep="/")
mesa.fam.file <- paste(f.dir, "MESA_Family_AA_YAO_20171031_T2D_sidno.txt", sep="/")
chs.file <- paste(f.dir, "CHS_FLOYD_201708018_T2DPhenotype.PED", sep="/")
genoa.file <- paste(f.dir, "GENOA_AA_BIELAK_20171030_T2D.ped", sep="/")
dhs.file <- paste(f.dir, "2a_dbGaP_SubjectPhenotypesDS_AACAC.txt", sep="/")
dhs.pedigree <- paste(f.dir, "6a_dbGaP_PedigreeDS_AACAC_revised072817.txt", sep="/")

map <- read.table(sample.file,
                  header=TRUE, as.is=T, sep="\t") #n=54499 & 16 variables:


raw.MGH =read.table(MGH.file,
                    header=T,sep='\t',as.is=T,fill = TRUE) #n=1025 & 21 variables
raw.MGHPed =read.table(MGH.pedigree,
                       header=T,sep='\t',as.is=T,fill = TRUE) #n=8

raw.VUdaw = read.table(VU.file,
                       header=T,sep='\t',as.is=T) #n=1134
raw.VUdawPed =read.table(VU.pedigree,
                         header=T,sep='\t',as.is=T,fill = TRUE) #n=2 & 6 variables

afp <- read.table(afp.file,
                  header=T,sep='\t',as.is=T) #n=128

afvub = read.table(afvub.file,
                   header=T,sep='\t',as.is=T) #n=171 & 17 variables

afccaf = read.table(afccaf.file, header=T,sep='\t',as.is=T) #n=363

raw.HVH = read.table(HVH.file,header=T,sep='\t',as.is=T) #n=708

aric_aa = read.table(ARIC.AF.file,header=T,sep='\t',as.is=T) ## n=353
aric_ea = read.table(ARIC.EU.file,header=T,sep='\t',as.is=T) ## n=3877
aric_ea$ancestry = "EU"
aric_aa$ancestry = "AF"
aric <- rbind(aric_ea,aric_aa) #n=4230

copd.c1 <- read.table(copd.1.file,
                      skip=10, header=TRUE, sep="\t", fill=TRUE) #n=10099
copd.c2 <- read.table(copd.2.file,
                      skip=10, header=TRUE, sep="\t", fill=TRUE) #n=272
copd <- rbind(copd.c1, copd.c2)  #n=10371

amish = read.table(amish.file, header=T,sep='\t',as.is=T) #n=1013 & 15 variables

cfs = read.table(cfs.file, header=T,sep='\t',as.is=T) #n=2531

sas = read.table(sas.file, header=T,sep='',as.is=T) #n=3470

jhs = read.table(jhs.file,
                 header=T,sep='\t',as.is=T,na.string="x") #n=3406 & 15 variables

fhs = read.table(fhs.file,
                 header=T,sep='\t',as.is=T)  #n=14329 & 14 variables

whi.eu = read.table(whi.EU.file, header=F,sep='\t',as.is=T)#n=8983
whi.aa = read.table(whi.AF.file, header=F,sep='\t',as.is=T)#n=1442
whi.as = read.table(whi.AS.file, header=F,sep='\t',as.is=T)#n=203
whi.ha = read.table(whi.HS.file, header=F,sep='\t',as.is=T)#n=310
whi.eu$ancestry = "EU"
whi.aa$ancestry = "AF"
whi.as$ancestry = "AS"
whi.ha$ancestry = "HS"
whi <- rbind(whi.eu,whi.ha,whi.as,whi.aa) #n=10938
colnames(whi) <- c("FamilyID","individual_id","PaternalID","MaternalID","sex","t2d",
                   "last_exam_age","last_exam_bmi","last_exam_fg","last_exam_hba1c",
                   "last_exam_t2d_treatment","last_exam_visit","t2d_age","t2d_bmi",
                   "sequenced","ascertainment_criteria","ancestry")

gensalt = read.table(gensalt.EU.file, header=T,sep='',as.is=T) #n=1906

genestarAA = read.table(genestar.AF.file, header=T,sep='\t',as.is=T) #n=1433
genestarAA$ancestry = 'AF'
genestarEU = read.table(genestar.EU.file, header=T,sep='\t',as.is=T) #n=1774
genestarEU$ancestry = 'EU'
genestar <- rbind(genestarAA, genestarEU)  #n=3207
rm(genestarAA)
rm(genestarEU)

mesaAA = read.csv(mesa.AF.file, header=T,sep=',',na.strings = '.') #N=1676
mesaAA$ancestry = 'AF'
mesaEU = read.csv(mesa.EU.file, header=T,sep=',',na.strings = '.') #=2528
mesaEU$ancestry = 'EU'
mesaHS = read.csv(mesa.HS.file, header=T,sep=',',na.strings = '.') #N=1450
mesaHS$ancestry = 'HS'
mesaSA = read.csv(mesa.AS.file, header=T,sep=',',na.strings = '.') #N=775
mesaSA$ancestry = 'AS'
mesa <- rbind(mesaAA, mesaEU,mesaHS,mesaSA)  #n=6429
rm(mesaAA)
rm(mesaEU)
rm(mesaHS)
rm(mesaSA)

mesafam = read.table(mesa.fam.file, header=T,sep='\t',as.is=T) #n=1134

chs = read.table(chs.file, header=T,sep='\t',as.is=T) #n=3929

genoa = read.table(genoa.file, header=T,sep='\t',as.is=T) #n=1854

dhs = read.table(dhs.file, header=T,sep='\t',as.is=T) #n=405

dhsPed =read.table(dhs.pedigree,
                   header=T,sep='\t',as.is=T,fill = TRUE) #n=972
