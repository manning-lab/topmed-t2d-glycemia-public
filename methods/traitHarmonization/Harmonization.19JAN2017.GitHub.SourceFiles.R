## T2D Harmonization  Source Files
# list of files:
# sample_annotation.file <- "freeze5b_sample_annot_2017-12-01.txt"
# MGH.file <- "MGH_Parent_dbGaP_SubjectPhenotypesDS_v3.txt"
# VU.file <- "phs001032.v3.pht005675.v3.p2.c1.TOPMed_WGS_VUH_AF_Subject_Phenotypes.GRU-IRB.txt"
# afp.file <- "phs001024.v2.pht005693.v1.p1.c1.TOPMed_WGS_Partners_AFGen_Subject_Phenotypes.HMB.txt"
# afvub.file <- "VUH_Ablation_dbGaP_SubjectPhenotypesDS_v2.txt'"
# afccaf.file <- "CCAF_dbGaP_SubjectPhenotypesDS_v5.txt"
# HVH.file <- "HVH_FLOYD_20170818_T2DPhenotype.PED"
# ARIC.AF.file <- "ARIC_AA_NMAYA_20161220_T2D.ped"
# ARIC.EU.file <- "ARIC_EU_NMAYA_20161220_T2D.ped"
# copd.1.file <- "phs000179.v5.pht002239.v4.p2.c1.COPDGene_Subject_Phenotypes.HMB.txt"
# copd.2.file <- "phs000179.v5.pht002239.v4.p2.c2.COPDGene_Subject_Phenotypes.DS"
# amish.file <- "Amish_HuichunXu_May152017_T2D.ped"
# cfs.file <- "CFS_12May17_T2D_Freeze4.ped"
# sas.file <- "SAS_LIULIN_20170212_T2D.ped"
# jhs.file <- "JHS_T2Dincidence_16Aug2017.ped"
# fhs.file <- "FHS_T2D_20170214_CS.ped"
# whi.EU.file <- "WHI_EU_MP_20170726_T2D.ped"
# whi.AF.file <- "WHI_AA_MP_20170726_T2D.ped"
# whi.AS.file <- "WHI_AS_MP_20170726_T2D.ped"
# whi.HS.file <- "WHI_HA_MP_20170726_T2D.ped"
# gensalt.EU.file <- "GenSalt_EA_XuenanMi_20170614_T2D.ped"
# genestar.AF.file <- "GeneSTAR_AA_YANEK_20170720_T2D.ped"
# genestar.EU.file <- "GeneSTAR_EU_YANEK_20170720_T2D.ped"
# mesa.AF.file <- "MESA_AA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv"
# mesa.EU.file <- "MESA_EU_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv"
# mesa.HS.file <- "MESA_HA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv"
# mesa.AS.file <- "MESA_SA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv"
# mesa.fam.file <- "MESA_Family_AA_YAO_20171031_T2D_sidno.txt"
# chs.file <- "CHS_FLOYD_201708018_T2DPhenotype.PED"
# genoa.file <- "GENOA_AA_BIELAK_20171030_T2D.ped"
# dhs.file <- "2a_dbGaP_SubjectPhenotypesDS_AACAC.txt"

##### Questions for Wessel
# 2 amish files, only using one (not using Amish_HuichunXu_May152017_T2D.DD.ncbi_enc)
# what AACAC file to use? need 2a and 6a
# MESA files dont match, does the new file have all of the MESA ancestry groups?
# gensalt files differ?
# difference between SA files? what does DD mean?
# no thrv files included
# WHI files differ


map <- read.table("freeze5b_sample_annot_2017-12-01.txt",
                  header=TRUE, as.is=T, sep="\t") #n=54499 & 16 variables:


raw.MGH =read.table('MGH_Parent_dbGaP_SubjectPhenotypesDS_v3.txt',
                    header=T,sep='\t',as.is=T,fill = TRUE) #n=1025 & 21 variables
raw.MGHPed =read.table('MGH_Parent_PedigreeDS_v1.txt',
                       header=T,sep='\t',as.is=T,fill = TRUE) #n=8

raw.VUdaw = read.table('phs001032.v3.pht005675.v3.p2.c1.TOPMed_WGS_VUH_AF_Subject_Phenotypes.GRU-IRB.txt',
                       header=T,sep='\t',as.is=T) #n=1134
raw.VUdawPed =read.table('phs001032.v3.pht007135.v1.p2.TOPMed_WGS_VUH_AF_Pedigree.MULTI.txt',
                         header=T,sep='\t',as.is=T,fill = TRUE) #n=2 & 6 variables

afp <- read.table('phs001024.v2.pht005693.v1.p1.c1.TOPMed_WGS_Partners_AFGen_Subject_Phenotypes.HMB.txt',
                  header=T,sep='\t',as.is=T) #n=128

afvub = read.table('VUH_Ablation_dbGaP_SubjectPhenotypesDS_v2.txt',
                   header=T,sep='\t',as.is=T) #n=171 & 17 variables

afccaf = read.table('CCAF_dbGaP_SubjectPhenotypesDS_v5.txt', header=T,sep='\t',as.is=T) #n=363

raw.HVH = read.table('HVH_FLOYD_20170818_T2DPhenotype.PED',header=T,sep='\t',as.is=T) #n=708

aric_aa = read.table('ARIC_AA_NMAYA_20161220_T2D.ped',header=T,sep='\t',as.is=T) ## n=353
aric_ea = read.table('ARIC_EU_NMAYA_20161220_T2D.ped',header=T,sep='\t',as.is=T) ## n=3877
aric_ea$ancestry = "EU"
aric_aa$ancestry = "AF"
aric <- rbind(aric_ea,aric_aa) #n=4230

copd.c1 <- read.table("phs000179.v5.pht002239.v4.p2.c1.COPDGene_Subject_Phenotypes.HMB.txt",
                      skip=10, header=TRUE, sep="\t", fill=TRUE) #n=10099
copd.c2 <- read.table("phs000179.v5.pht002239.v4.p2.c2.COPDGene_Subject_Phenotypes.DS-CS-RD.txt",
                      skip=10, header=TRUE, sep="\t", fill=TRUE) #n=272
copd <- rbind(copd.c1, copd.c2)  #n=10371

amish = read.table('Amish_HuichunXu_May152017_T2D.ped', header=T,sep='\t',as.is=T) #n=1013 & 15 variables

cfs = read.table('CFS_12May17_T2D_Freeze4.ped', header=T,sep='\t',as.is=T) #n=2531

sas = read.table('SAS_LIULIN_20170212_T2D.ped', header=T,sep='',as.is=T) #n=3470

jhs = read.table('JHS_T2Dincidence_16Aug2017.ped',
                 header=T,sep='\t',as.is=T,na.string="x") #n=3406 & 15 variables

fhs = read.table('FHS_T2D_20170214_CS.ped',
                 header=T,sep='\t',as.is=T)  #n=14329 & 14 variables

whi.eu = read.table('WHI_EU_MP_20170726_T2D.ped', header=F,sep='\t',as.is=T)#n=8983
whi.aa = read.table('WHI_AA_MP_20170726_T2D.ped', header=F,sep='\t',as.is=T)#n=1442
whi.as = read.table('WHI_AS_MP_20170726_T2D.ped', header=F,sep='\t',as.is=T)#n=203
whi.ha = read.table('WHI_HA_MP_20170726_T2D.ped', header=F,sep='\t',as.is=T)#n=310
whi.eu$ancestry = "EU"
whi.aa$ancestry = "AF"
whi.as$ancestry = "AS"
whi.ha$ancestry = "HS"
whi <- rbind(whi.eu,whi.ha,whi.as,whi.aa) #n=10938
colnames(whi) <- c("FamilyID","individual_id","PaternalID","MaternalID","sex","t2d",
                   "last_exam_age","last_exam_bmi","last_exam_fg","last_exam_hba1c",
                   "last_exam_t2d_treatment","last_exam_visit","t2d_age","t2d_bmi",
                   "sequenced","ascertainment_criteria","ancestry")

gensalt = read.table('GenSalt_EA_XuenanMi_20170614_T2D.ped', header=T,sep='',as.is=T) #n=1906

genestarAA = read.table('GeneSTAR_AA_YANEK_20170720_T2D.ped', header=T,sep='\t',as.is=T) #n=1433
genestarAA$ancestry = 'AF'
genestarEU = read.table('GeneSTAR_EU_YANEK_20170720_T2D.ped', header=T,sep='\t',as.is=T) #n=1774
genestarEU$ancestry = 'EU'
genestar <- rbind(genestarAA, genestarEU)  #n=3207
rm(genestarAA)
rm(genestarEU)

mesaAA = read.csv('MESA_AA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv', header=T,sep=',',na.strings = '.') #N=1676
mesaAA$ancestry = 'AF'
mesaEU = read.csv('MESA_EU_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv', header=T,sep=',',na.strings = '.') #=2528
mesaEU$ancestry = 'EU'
mesaHS = read.csv('MESA_HA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv', header=T,sep=',',na.strings = '.') #N=1450
mesaHS$ancestry = 'HS'
mesaSA = read.csv('MESA_SA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv', header=T,sep=',',na.strings = '.') #N=775
mesaSA$ancestry = 'AS'
mesa <- rbind(mesaAA, mesaEU,mesaHS,mesaSA)  #n=6429
rm(mesaAA)
rm(mesaEU)
rm(mesaHS)
rm(mesaSA)

mesafam = read.table('MESA_Family_AA_YAO_20171031_T2D_sidno.txt', header=T,sep='\t',as.is=T) #n=1134

chs = read.table('CHS_FLOYD_201708018_T2DPhenotype.PED', header=T,sep='\t',as.is=T) #n=3929

genoa = read.table('GENOA_AA_BIELAK_20171030_T2D.ped', header=T,sep='\t',as.is=T) #n=1854

dhs = read.table('2a_dbGaP_SubjectPhenotypesDS_AACAC.txt', header=T,sep='\t',as.is=T) #n=405

dhsPed =read.table('6a_dbGaP_PedigreeDS_AACAC_revised072817.txt',
                   header=T,sep='\t',as.is=T,fill = TRUE) #n=972
