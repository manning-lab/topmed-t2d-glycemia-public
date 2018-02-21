## GLY Harmonization freeze5b
# Studies:
# fhs
# jhs
# sas
# cfs
# amish
# gensalt
# aric
# whi
# genestar
# chs
# mesa
# genoa

args <- commandArgs(trailingOnly=T)
f.dir <- args[1]
source.file <- args[2]


#### testing inputs ####
# f.dir <- "/Users/tmajaria/Documents/projects/topmed/data/freeze5b_phenotypes"
# source.file <- "/Users/tmajaria/Documents/projects/topmed/code/topmed-t2d-glycemia-public/methods/traitHarmonization/gly_traits_harmonization_021318_filepaths.R"
########################

# load the data
source(source.file)
dat <- get_pheno_data(f.dir)

# get data to the right variable names to match script
linker <- dat$linker 
fhs <- dat$fhs 
share <- dat$share
jhs <- dat$jhs 
sas <- dat$sas
cfs <- dat$cfs
amish <- dat$amish
gensalt <- dat$gensalt
aric_ea <- dat$aric_ea
aric_aa <- dat$aric_aa
whi_ha <- dat$whi_ha
whi_ea <- dat$whi_ea
whi_as <- dat$whi_as
whi_aa <- dat$whi_aa
gs_aa <- dat$gs_aa
gs_ea <- dat$gs_ea
chs <- dat$chs
mesa_ha <- dat$mesa_ha
mesa_ea <- dat$mesa_ea
mesa_sa <- dat$mesa_sa
mesa_aa <- dat$mesa_aa
genoa <- dat$genoa

rm(dat)

############################## FHS ##############################
############################## FHS ##############################

fhs<-merge(fhs,linker[which(linker$study=='FHS'),],by.x='shareid',by.y="study_subject_id")

colnames(fhs)[which(colnames(fhs)=="shareid")]<-"Individual_ID"

colnames(fhs)[which(colnames(fhs)=="sample.id")]<-"TOPMEDID"

fhs.ped<-share[,c("famid", "shareid","mshare", "fshare")]

colnames(fhs.ped)<-c("Family_ID","Individual_ID","Mother_ID","Father_ID")
fhs<-merge(fhs,fhs.ped,by="Individual_ID",all.x=T)

fhs$STUDY_TOPMEDID<-paste("FHS",fhs$TOPMEDID,sep="_")
fhs$STUDY_ANCESTRY<-ifelse(fhs$idtype==0,'FHS_1_EU',ifelse(fhs$idtype==1,'FHS_2_EU','FHS_3_EU'))


####check T2D coding
table(fhs$T2D,useNA="always")
table(fhs$T2D_FG,useNA="always")
table(fhs$T2D_FI,useNA="always")
table(fhs$T2D_HbA1c,useNA="always")

###recode T2D=1 if T2D is NA 
for (j in 1:nrow(fhs)){
	fhs$T2D_FG[j]<-ifelse(is.na(fhs$T2D_FG[j])==T,1,fhs$T2D_FG[j])
	fhs$T2D_FI[j]<-ifelse(is.na(fhs$T2D_FI[j])==T,1,fhs$T2D_FI[j])
	fhs$T2D_HbA1c[j]<-ifelse(is.na(fhs$T2D_HbA1c[j])==T,1,fhs$T2D_HbA1c[j])
	} 


### set FG to NA if FastingGlucose>=7&T2D_FG==1
fhs$FastingGlucose<-ifelse(fhs$FastingGlucose>=7&fhs$T2D_FG==1,NA,fhs$FastingGlucose)
####set A1C>7 to NA
fhs$HbA1c<-ifelse(fhs$HbA1c>=6.5&fhs$T2D_HbA1c==1,NA,fhs$HbA1c)
fhs<-fhs[,c("TOPMEDID" ,              "Individual_ID",          "Family_ID" ,            
"Mother_ID" ,             "Father_ID"     ,         "T2D"        ,           
 "sex"       ,             "FastingGlucose",         "T2D_FG"     ,           
 "age_FG"     ,            "BMI_FG"         ,        "FastingInsulin",        
"T2D_FI"       ,          "age_FI"           ,      "BMI_FI"          ,      
 "HbA1c"        ,          "T2D_HbA1c"        ,      "age_HbA1c"       ,      
"BMI_HbA1c"      ,        "FastingGlucose_HbA1c",   "TwoHourGlucose_HbA1c",  
"Hb_HbA1c"        ,       "MCV_HbA1c"            ,  "MCH_HbA1c"            , 
"MCHC_HbA1c"       ,      "Fe_HbA1c"              , "Ferritin_HbA1c"        ,
"TSAT_HbA1c"        ,     "sequenced"              ,"ascertainment_criteria",
"consent"            ,    "topmed_phs"             ,"study"                 ,
"topmed_project"      ,   "CENTER"                 ,"geno.cntl"             ,
"TRIO.dups"            ,  "no.post"                ,"STUDY_TOPMEDID"        ,
"STUDY_ANCESTRY")]

############################## FHS ##############################
############################## FHS ##############################

############################## JHS ##############################
############################## JHS ##############################


#####JHS

jhs<-merge(jhs,linker[which(linker$study=='JHS'),],by.x='Individual_ID',by.y="study_subject_id")

colnames(jhs)[which(colnames(jhs)=="Gluc2Hr_HbA1c")]<-"TwoHourGlucose_HbA1c"
#colnames(jhs)[which(colnames(fhs)=="sample.id")]<-"NWDid"
jhs$FastingInsulin<-(jhs$FastingInsulin)/6


colnames(jhs)[which(colnames(jhs)=="Individual_ID")]<-"Individual_ID"

colnames(jhs)[which(colnames(jhs)=="sample.id")]<-"TOPMEDID"

jhs$STUDY_TOPMEDID<-paste("JHS",jhs$TOPMEDID,sep="_")
jhs$STUDY_ANCESTRY<-"JHS_AA"

jhs<-jhs[,names(fhs)]


####check T2D coding
table(jhs$T2D,useNA="always")
table(jhs$T2D_FG,useNA="always")
table(jhs$T2D_FI,useNA="always")
table(jhs$T2D_HbA1c,useNA="always")

###recode T2D=1 if T2D is NA 
for (j in 1:nrow(jhs)){
	jhs$T2D_FG[j]<-ifelse(is.na(jhs$T2D_FG[j])==T,1,jhs$T2D_FG[j])
	jhs$T2D_FI[j]<-ifelse(is.na(jhs$T2D_FI[j])==T,1,jhs$T2D_FI[j])
	jhs$T2D_HbA1c[j]<-ifelse(is.na(jhs$T2D_HbA1c[j])==T,1,jhs$T2D_HbA1c[j])
	} 


### set FG to NA if FastingGlucose>=7&T2D_FG==1
jhs$FastingGlucose<-ifelse(jhs$FastingGlucose>=7&jhs$T2D_FG==1,NA,jhs$FastingGlucose)
####set A1C>7 to NA
jhs$HbA1c<-ifelse(jhs$HbA1c>=6.5&jhs$T2D_HbA1c==1,NA,jhs$HbA1c)

############################## JHS ##############################
############################## JHS ##############################

############################## SAS ##############################
############################## SAS ##############################


###SAS

sas<-merge(sas,linker[which(linker$study=='SAS'),],by.x='SUBJECT_ID',by.y="study_subject_id")
sas<-sas[,-c(2)]
colnames(sas)[which(colnames(sas)=="SUBJECT_ID")]<-"Individual_ID"


#colnames(jhs)[which(colnames(jhs)=="Individual_ID")]<-"Individual_ID"

colnames(sas)[which(colnames(sas)=="sample.id")]<-"TOPMEDID"


sas$STUDY_TOPMEDID<-paste("SAS",sas$TOPMEDID,sep="_")
sas$STUDY_ANCESTRY<-"SAS_SA"

sas<-sas[,-c(30)]

sas<-sas[,names(fhs)]

####check T2D coding
table(sas$T2D,useNA="always")
table(sas$T2D_FG,useNA="always")
table(sas$T2D_FI,useNA="always")
table(sas$T2D_HbA1c,useNA="always")

###recode T2D
for (j in 1:nrow(sas)){
	sas$T2D_FG[j]<-ifelse(sas$T2D_FG[j]==-9|sas$T2D_FG[j]==0,1,sas$T2D_FG[j])
	sas$T2D_FI[j]<-ifelse(sas$T2D_FI[j]==-9|sas$T2D_FI[j]==0,1,sas$T2D_FI[j])
	sas$T2D_HbA1c[j]<-ifelse(is.na(sas$T2D_HbA1c[j])==T,1,sas$T2D_HbA1c[j])}


### set FG to NA if FastingGlucose>=7&T2D_FG==1
sas$FastingGlucose<-ifelse(sas$FastingGlucose>=7&sas$T2D_FG==1,NA,sas$FastingGlucose)
####set A1C>7 to NA
sas$HbA1c<-ifelse(sas$HbA1c>=6.5&sas$T2D_HbA1c==1,NA,sas$HbA1c)

for (j in 1:nrow(sas)){
	sas$FastingGlucose[j]<-ifelse(sas$T2D_FG[j]==2,NA,sas$FastingGlucose[j])
	}
for (j in 1:nrow(sas)){
	sas$FastingInsulin[j]<-ifelse(sas$T2D_FI[j]==2,NA,sas$FastingInsulin[j])
	}	

############################## SAS ##############################
############################## SAS ##############################

############################## CFS ##############################
############################## CFS ##############################

####cfs

cfs<-merge(cfs,linker[which(linker$study=='CFS'),],by.x='IID',by.y="study_subject_id")


colnames(cfs)[which(colnames(cfs)=="IID")]<-"Individual_ID"
colnames(cfs)[which(colnames(cfs)=="NWDid")]<-"TOPMEDID"
colnames(cfs)[which(colnames(cfs)=="FamID")]<-"Family_ID"
colnames(cfs)[which(colnames(cfs)=="FID")]<-"Father_ID"
colnames(cfs)[which(colnames(cfs)=="MID")]<-"Mother_ID"
colnames(cfs)[which(colnames(cfs)=="Sequenced")]<-"sequenced"
cfs<-cfs[,!(names(cfs)%in%c("sample.id","T2D_all_ages"))]

cfs$HbA1c<-NA
cfs$T2D_HbA1c<-1
cfs$age_HbA1c<-NA
cfs$BMI_HbA1c<-NA
cfs$FastingGlucose_HbA1c<-NA
cfs$TwoHourGlucose_HbA1c<-NA
cfs$Hb_HbA1c<-NA
cfs$MCV_HbA1c<-NA
cfs$MCH_HbA1c<-NA
cfs$MCHC_HbA1c<-NA
cfs$Fe_HbA1c<-NA
cfs$Ferritin_HbA1c<-NA
cfs$TSAT_HbA1c<-NA
cfs$ascertainment_criteria<-NA
#cfs$FHS_cohorts<-NA
#cfs$topmed_project<-"CFS"
#cfs$study<-"CFS"
cfs$STUDY_ANCESTRY<-ifelse(cfs$Population=="CFS-whites",paste("CFS","EU",sep="_"),paste("CFS","AA",sep="_"))
cfs$STUDY_TOPMEDID<-paste("CFS",cfs$TOPMEDID,sep="_")
cfs<-cfs[,!(names(cfs)%in%c("Population"))]

cfs<-cfs[,names(jhs)]

####check T2D coding
table(cfs$T2D,useNA="always")
table(cfs$T2D_FG,useNA="always")
table(cfs$T2D_FI,useNA="always")
table(cfs$T2D_HbA1c,useNA="always")

for (j in 1:nrow(cfs)){
	cfs$T2D_FG[j]<-ifelse(is.na(cfs$T2D_FG[j])==T|cfs$T2D_FG[j]==0,1,cfs$T2D_FG[j])
	cfs$T2D_FI[j]<-ifelse(is.na(cfs$T2D_FI[j])==T|cfs$T2D_FI[j]==0,1,cfs$T2D_FI[j])
	cfs$T2D[j]<-ifelse(is.na(cfs$T2D[j])==T,0,cfs$T2D[j])

	#cfs$T2D_HbA1c[j]<-ifelse(is.na(cfs$T2D_HbA1c[j])==T,1,cfs$T2D_HbA1c[j])
	}

for (j in 1:nrow(cfs)){
	cfs$FastingGlucose[j]<-ifelse(cfs$T2D_FG[j]==2,NA,cfs$FastingGlucose[j])
	}
for (j in 1:nrow(cfs)){
	cfs$FastingInsulin[j]<-ifelse(cfs$T2D_FI[j]==2,NA,cfs$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
cfs$FastingGlucose<-ifelse(cfs$FastingGlucose>=7&cfs$T2D_FG==1,NA,cfs$FastingGlucose)
####set A1C>7 to NA
cfs$HbA1c<-ifelse(cfs$HbA1c>=6.5&cfs$T2D_HbA1c==1,NA,cfs$HbA1c)

############################## CFS ##############################
############################## CFS ##############################

############################## Amish ##############################
############################## Amish ##############################



####amish

amish<-merge(amish,linker[which(linker$study=='Amish'),],by.x='dbgap_id',by.y="study_subject_id")



names(amish)<-c("Individual_ID","Family_ID","Father_ID","Mother_ID","sex","T2D","FastingGlucose","T2D_FG","age_FG","BMI_FG","amishstudy_fg","FastingInsulin","T2D_FI","age_FI","BMI_FI","amishstudy_fi","HbA1c","T2D_HbA1c","age_HbA1c", "BMI_HbA1c","FastingGlucose_HbA1c","TwoHourGlucose_HbA1c","sequenced","TOPMEDID","consent"   ,           "topmed_phs"  ,         "study" ,"topmed_project",       "CENTER"   ,            "geno.cntl","TRIO.dups" ,           "no.post" )
amish<-amish[,!(names(amish)%in%c("amishstudy_fi",'amishstudy_fg'))]
amish$Hb_HbA1c<-NA
amish$MCV_HbA1c<-NA
amish$MCH_HbA1c<-NA
amish$MCHC_HbA1c<-NA
amish$Fe_HbA1c<-NA
amish$Ferritin_HbA1c<-NA
amish$TSAT_HbA1c<-NA


amish$STUDY_TOPMEDID<-paste("Amish",amish$TOPMEDID,sep="_")
amish$STUDY_ANCESTRY<-"Amish_EU"
amish$ascertainment_criteria<-NA

amish<-amish[,names(jhs)]


####check T2D coding
table(amish$T2D,useNA="always")
table(amish$T2D_FG,useNA="always")
table(amish$T2D_FI,useNA="always")
table(amish$T2D_HbA1c,useNA="always")


amish$T2D<-ifelse(amish$T2D==1,2,0)

amish$T2D_FG<-ifelse(is.na(amish$T2D_FG)==T,1,1)

amish$T2D_FI<-ifelse(is.na(amish$T2D_FI)==T,1,1)

amish$T2D_HbA1c<-ifelse(is.na(amish$T2D_HbA1c)==T,1,1)


for (j in 1:nrow(amish)){
	amish$FastingGlucose[j]<-ifelse(amish$T2D_FG[j]==2,NA,amish$FastingGlucose[j])
	}
for (j in 1:nrow(amish)){
	amish$FastingInsulin[j]<-ifelse(amish$T2D_FI[j]==2,NA,amish$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
amish$FastingGlucose<-ifelse(amish$FastingGlucose>=7&amish$T2D_FG==1,NA,amish$FastingGlucose)
####set A1C>7 to NA
amish$HbA1c<-ifelse(amish$HbA1c>=6.5&amish$T2D_HbA1c==1,NA,amish$HbA1c)

############################## Amish ##############################
############################## Amish ##############################


############################## Gensalt ##############################
############################## Gensalt ##############################


####Gensalt

gensalt<-merge(gensalt,linker[which(linker$study=='GenSalt'),],by.x='Individual_ID',by.y="study_subject_id")

colnames(gensalt)[which(colnames(gensalt)=="SEX")]<-"sex"

colnames(gensalt)[which(colnames(gensalt)=="Age_FG")]<-"age_FG"
colnames(gensalt)[which(colnames(gensalt)=="Age_FI")]<-"age_FI"
colnames(gensalt)[which(colnames(gensalt)=="Age_HbA1c")]<-"age_HbA1c"

colnames(gensalt)[which(colnames(gensalt)=="sample.id")]<-"TOPMEDID"


gensalt<-gensalt[,!(names(gensalt)%in%c('Baseline_SBP', 'Baseline_DBP' ,'Proband'))]

gensalt$STUDY_TOPMEDID<-paste("GenSalt",gensalt$TOPMEDID,sep="_")
gensalt$STUDY_ANCESTRY<-"GenSalt_AS"
gensalt$ascertainment_criteria<-NA



gensalt<-gensalt[,names(jhs)]


table(gensalt$T2D,useNA="always")
table(gensalt$T2D_FG,useNA="always")
table(gensalt$T2D_FI,useNA="always")
table(gensalt$T2D_HbA1c,useNA="always")
gensalt$T2D<-ifelse(gensalt$T2D==9,-9,gensalt$T2D)
gensalt$T2D_FI<-1
gensalt$T2D_HbA1c<-1

gensalt0<-gensalt[which(gensalt$FastingGlucose=='.'),]


gensalt1<-gensalt[which(!(gensalt$TOPMEDID%in%gensalt0$TOPMEDID)),]

gensalt1$FastingGlucose<- as.numeric(as.character(gensalt1$FastingGlucose))


gensalt0$FastingGlucose<-NA


gensalt<-rbind(gensalt1,gensalt0)





for (j in 1:nrow(gensalt)){
	gensalt$FastingGlucose[j]<-ifelse(gensalt$T2D_FG[j]==2,NA,gensalt$FastingGlucose[j])
	}
for (j in 1:nrow(gensalt)){
	gensalt$FastingInsulin[j]<-ifelse(gensalt$T2D_FI[j]==2,NA,gensalt$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
gensalt$FastingGlucose<-ifelse(gensalt$FastingGlucose>=7&gensalt$T2D_FG==1,NA,gensalt$FastingGlucose)
####set A1C>7 to NA
gensalt$HbA1c<-ifelse(gensalt$HbA1c>=6.5&gensalt$T2D_HbA1c==1,NA,gensalt$HbA1c)

############################## Gensalt ##############################
############################## Gensalt ##############################

############################## ARIC ##############################
############################## ARIC ##############################



#####ARIC

###ARIC_EA



aric_ea<-merge(aric_ea,linker[which(linker$study=='ARIC'),],by.x='Individual_ID',by.y="study_subject_id")
colnames(aric_ea)[which(colnames(aric_ea)=="Diabetes_Status")]<-"T2D"

colnames(aric_ea)[which(colnames(aric_ea)=="sample.id")]<-"TOPMEDID"

aric_ea$TwoHourGlucose_HbA1c<-NA
aric_ea$MCH_HbA1c<-NA
aric_ea$MCHC_HbA1c<-NA
aric_ea$Fe_HbA1c<-NA
aric_ea$Ferritin_HbA1c<-NA
aric_ea$TSAT_HbA1c<-NA
aric_ea$sequenced<-NA
aric_ea$ascertainment_criteria<-NA
aric_ea$STUDY_TOPMEDID<-paste("ARIC",aric_ea$TOPMEDID,sep="_")
aric_ea$STUDY_ANCESTRY<-"ARIC_EU"


aric_ea<-aric_ea[,names(jhs)]

table(aric_ea$T2D,useNA="always")
table(aric_ea$T2D_FG,useNA="always")
table(aric_ea$T2D_FI,useNA="always")
table(aric_ea$T2D_HbA1c,useNA="always")
aric_ea$T2D_FG<-aric_ea$T2D_FG+1
aric_ea$T2D_FI<-aric_ea$T2D_FI+1
aric_ea$T2D_HbA1c<-aric_ea$T2D_HbA1c+1
for (j in 1:nrow(aric_ea)){
	aric_ea$T2D_FG[j]<-ifelse(is.na(aric_ea$T2D_FG[j])==T,1,aric_ea$T2D_FG[j])
	aric_ea$T2D_FI[j]<-ifelse(is.na(aric_ea$T2D_FI[j])==T,1,aric_ea$T2D_FI[j])
	aric_ea$T2D_HbA1c[j]<-ifelse(is.na(aric_ea$T2D_HbA1c[j])==T,1,aric_ea$T2D_HbA1c[j])
	} 

for (j in 1:nrow(aric_ea)){
	aric_ea$FastingGlucose[j]<-ifelse(aric_ea$T2D_FG[j]==2,NA,aric_ea$FastingGlucose[j])
	}
for (j in 1:nrow(aric_ea)){
	aric_ea$FastingInsulin[j]<-ifelse(aric_ea$T2D_FI[j]==2,NA,aric_ea$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
aric_ea$FastingGlucose<-ifelse(aric_ea$FastingGlucose>=7&aric_ea$T2D_FG==1,NA,aric_ea$FastingGlucose)
####set A1C>7 to NA
aric_ea$HbA1c<-ifelse(aric_ea$HbA1c>=6.5&aric_ea$T2D_HbA1c==1,NA,aric_ea$HbA1c)


###ARIC_AA



aric_aa<-merge(aric_aa,linker[which(linker$study=='ARIC'),],by.x='Individual_ID',by.y="study_subject_id")
colnames(aric_aa)[which(colnames(aric_aa)=="Diabetes_Status")]<-"T2D"

colnames(aric_aa)[which(colnames(aric_aa)=="sample.id")]<-"TOPMEDID"

aric_aa$TwoHourGlucose_HbA1c<-NA
aric_aa$MCH_HbA1c<-NA
aric_aa$MCHC_HbA1c<-NA
aric_aa$Fe_HbA1c<-NA
aric_aa$Ferritin_HbA1c<-NA
aric_aa$TSAT_HbA1c<-NA
aric_aa$sequenced<-NA
aric_aa$ascertainment_criteria<-NA
aric_aa$STUDY_TOPMEDID<-paste("ARIC",aric_aa$TOPMEDID,sep="_")
aric_aa$STUDY_ANCESTRY<-"ARIC_AA"


aric_aa<-aric_aa[,names(jhs)]

table(aric_aa$T2D,useNA="always")
table(aric_aa$T2D_FG,useNA="always")
table(aric_aa$T2D_FI,useNA="always")
table(aric_aa$T2D_HbA1c,useNA="always")
aric_aa$T2D_FG<-aric_aa$T2D_FG+1
aric_aa$T2D_FI<-aric_aa$T2D_FI+1
aric_aa$T2D_HbA1c<-aric_aa$T2D_HbA1c+1
for (j in 1:nrow(aric_aa)){
	aric_aa$T2D_FG[j]<-ifelse(is.na(aric_aa$T2D_FG[j])==T,1,aric_aa$T2D_FG[j])
	aric_aa$T2D_FI[j]<-ifelse(is.na(aric_aa$T2D_FI[j])==T,1,aric_aa$T2D_FI[j])
	aric_aa$T2D_HbA1c[j]<-ifelse(is.na(aric_aa$T2D_HbA1c[j])==T,1,aric_aa$T2D_HbA1c[j])
	} 

for (j in 1:nrow(aric_aa)){
	aric_aa$FastingGlucose[j]<-ifelse(aric_aa$T2D_FG[j]==2,NA,aric_aa$FastingGlucose[j])
	}
for (j in 1:nrow(aric_aa)){
	aric_aa$FastingInsulin[j]<-ifelse(aric_aa$T2D_FI[j]==2,NA,aric_aa$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
aric_aa$FastingGlucose<-ifelse(aric_aa$FastingGlucose>=7&aric_aa$T2D_FG==1,NA,aric_aa$FastingGlucose)
####set A1C>7 to NA
aric_aa$HbA1c<-ifelse(aric_aa$HbA1c>=6.5&aric_aa$T2D_HbA1c==1,NA,aric_aa$HbA1c)

aric<-rbind(aric_ea,aric_aa)

############################## ARIC ##############################
############################## ARIC ##############################


############################## WHI ##############################
############################## WHI ##############################


#####WHI




###HA

names(whi_ha)<-c("Family_ID","Individual_ID","Father_ID","Mother_ID","sex","T2D","FastingGlucose","T2D_FG","age_FG","BMI_FG"         ,        "FastingInsulin",        
"T2D_FI"       ,          "age_FI"           ,      "BMI_FI","HbA1c"        ,          "T2D_HbA1c"        ,      "age_HbA1c"       ,      
"BMI_HbA1c"      ,        "FastingGlucose_HbA1c",   "TwoHourGlucose_HbA1c",  
"Hb_HbA1c"        ,       "MCV_HbA1c"            ,  "MCH_HbA1c"            , 
"MCHC_HbA1c"       ,      "Fe_HbA1c"              , "Ferritin_HbA1c"        ,
"TSAT_HbA1c",'sequenced',             'ascertainment_criteria')



whi_ha<-merge(whi_ha,linker[which(linker$study=='WHI'),],by.x='Individual_ID',by.y="study_subject_id")


colnames(whi_ha)[which(colnames(whi_ha)=="sample.id")]<-"TOPMEDID"




whi_ha$STUDY_TOPMEDID<-paste("WHI",whi_ha$TOPMEDID,sep="_")
whi_ha$STUDY_ANCESTRY<-"WHI_HA"

whi_ha<-whi_ha[,names(jhs)]



table(whi_ha$T2D,useNA="always")
table(whi_ha$T2D_FG,useNA="always")
table(whi_ha$T2D_FI,useNA="always")
table(whi_ha$T2D_HbA1c,useNA="always")

for (j in 1:nrow(whi_ha)){
	whi_ha$T2D_FG[j]<-ifelse(is.na(whi_ha$T2D_FG[j])==T,1,whi_ha$T2D_FG[j])
	whi_ha$T2D_FI[j]<-ifelse(is.na(whi_ha$T2D_FI[j])==T,1,whi_ha$T2D_FI[j])
	whi_ha$T2D_HbA1c[j]<-ifelse(is.na(whi_ha$T2D_HbA1c[j])==T,1,whi_ha$T2D_HbA1c[j])
	whi_ha$T2D[j]<-ifelse(is.na(whi_ha$T2D[j])==T,0,whi_ha$T2D[j])
	
	} 

for (j in 1:nrow(whi_ha)){
	whi_ha$FastingGlucose[j]<-ifelse(whi_ha$T2D_FG[j]==2,NA,whi_ha$FastingGlucose[j])
	}
for (j in 1:nrow(whi_ha)){
	whi_ha$FastingInsulin[j]<-ifelse(whi_ha$T2D_FI[j]==2,NA,whi_ha$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
whi_ha$FastingGlucose<-ifelse(whi_ha$FastingGlucose>=7&whi_ha$T2D_FG==1,NA,whi_ha$FastingGlucose)
####set A1C>7 to NA
whi_ha$HbA1c<-ifelse(whi_ha$HbA1c>=6.5&whi_ha$T2D_HbA1c==1,NA,whi_ha$HbA1c)

####EU

names(whi_ea)<-c("Family_ID","Individual_ID","Father_ID","Mother_ID","sex","T2D","FastingGlucose","T2D_FG","age_FG","BMI_FG"         ,        "FastingInsulin",        
"T2D_FI"       ,          "age_FI"           ,      "BMI_FI","HbA1c"        ,          "T2D_HbA1c"        ,      "age_HbA1c"       ,      
"BMI_HbA1c"      ,        "FastingGlucose_HbA1c",   "TwoHourGlucose_HbA1c",  
"Hb_HbA1c"        ,       "MCV_HbA1c"            ,  "MCH_HbA1c"            , 
"MCHC_HbA1c"       ,      "Fe_HbA1c"              , "Ferritin_HbA1c"        ,
"TSAT_HbA1c",'sequenced',             'ascertainment_criteria')

whi_ea<-merge(whi_ea,linker[which(linker$study=='WHI'),],by.x='Individual_ID',by.y="study_subject_id")


colnames(whi_ea)[which(colnames(whi_ea)=="sample.id")]<-"TOPMEDID"




whi_ea$STUDY_TOPMEDID<-paste("WHI",whi_ea$TOPMEDID,sep="_")
whi_ea$STUDY_ANCESTRY<-"WHI_EU"

whi_ea<-whi_ea[,names(jhs)]



table(whi_ea$T2D,useNA="always")
table(whi_ea$T2D_FG,useNA="always")
table(whi_ea$T2D_FI,useNA="always")
table(whi_ea$T2D_HbA1c,useNA="always")

	whi_ea$T2D_FG<-ifelse(is.na(whi_ea$T2D_FG)==T,1,whi_ea$T2D_FG)
	whi_ea$T2D_FI<-ifelse(is.na(whi_ea$T2D_FI)==T,1,whi_ea$T2D_FI)
	whi_ea$T2D_HbA1c<-ifelse(is.na(whi_ea$T2D_HbA1c)==T,1,whi_ea$T2D_HbA1c)
	whi_ea$T2D<-ifelse(is.na(whi_ea$T2D)==T,0,whi_ea$T2D)
	


for (j in 1:nrow(whi_ea)){
	whi_ea$FastingGlucose[j]<-ifelse(whi_ea$T2D_FG[j]==2,NA,whi_ea$FastingGlucose[j])
	}
for (j in 1:nrow(whi_ea)){
	whi_ea$FastingInsulin[j]<-ifelse(whi_ea$T2D_FI[j]==2,NA,whi_ea$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
whi_ea$FastingGlucose<-ifelse(whi_ea$FastingGlucose>=7&whi_ea$T2D_FG==1,NA,whi_ea$FastingGlucose)
####set A1C>7 to NA
whi_ea$HbA1c<-ifelse(whi_ea$HbA1c>=6.5&whi_ea$T2D_HbA1c==1,NA,whi_ea$HbA1c)

###AS


names(whi_as)<-c("Family_ID","Individual_ID","Father_ID","Mother_ID","sex","T2D","FastingGlucose","T2D_FG","age_FG","BMI_FG"         ,        "FastingInsulin",        
"T2D_FI"       ,          "age_FI"           ,      "BMI_FI","HbA1c"        ,          "T2D_HbA1c"        ,      "age_HbA1c"       ,      
"BMI_HbA1c"      ,        "FastingGlucose_HbA1c",   "TwoHourGlucose_HbA1c",  
"Hb_HbA1c"        ,       "MCV_HbA1c"            ,  "MCH_HbA1c"            , 
"MCHC_HbA1c"       ,      "Fe_HbA1c"              , "Ferritin_HbA1c"        ,
"TSAT_HbA1c",'sequenced',             'ascertainment_criteria')

whi_as<-merge(whi_as,linker[which(linker$study=='WHI'),],by.x='Individual_ID',by.y="study_subject_id")


colnames(whi_as)[which(colnames(whi_as)=="sample.id")]<-"TOPMEDID"




whi_as$STUDY_TOPMEDID<-paste("WHI",whi_as$TOPMEDID,sep="_")
whi_as$STUDY_ANCESTRY<-"WHI_AS"

whi_as<-whi_as[,names(jhs)]



table(whi_as$T2D,useNA="always")
table(whi_as$T2D_FG,useNA="always")
table(whi_as$T2D_FI,useNA="always")
table(whi_as$T2D_HbA1c,useNA="always")

	whi_as$T2D_FG<-ifelse(is.na(whi_as$T2D_FG)==T,1,whi_as$T2D_FG)
	whi_as$T2D_FI<-ifelse(is.na(whi_as$T2D_FI)==T,1,whi_as$T2D_FI)
	whi_as$T2D_HbA1c<-ifelse(is.na(whi_as$T2D_HbA1c)==T,1,whi_as$T2D_HbA1c)
	whi_as$T2D<-ifelse(is.na(whi_as$T2D)==T,0,whi_as$T2D)
	


for (j in 1:nrow(whi_as)){
	whi_as$FastingGlucose[j]<-ifelse(whi_as$T2D_FG[j]==2,NA,whi_as$FastingGlucose[j])
	}
for (j in 1:nrow(whi_as)){
	whi_as$FastingInsulin[j]<-ifelse(whi_as$T2D_FI[j]==2,NA,whi_as$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
whi_as$FastingGlucose<-ifelse(whi_as$FastingGlucose>=7&whi_as$T2D_FG==1,NA,whi_as$FastingGlucose)
####set A1C>7 to NA
whi_as$HbA1c<-ifelse(whi_as$HbA1c>=6.5&whi_as$T2D_HbA1c==1,NA,whi_as$HbA1c)

###AA

names(whi_aa)<-c("Family_ID","Individual_ID","Father_ID","Mother_ID","sex","T2D","FastingGlucose","T2D_FG","age_FG","BMI_FG"         ,        "FastingInsulin",        
"T2D_FI"       ,          "age_FI"           ,      "BMI_FI","HbA1c"        ,          "T2D_HbA1c"        ,      "age_HbA1c"       ,      
"BMI_HbA1c"      ,        "FastingGlucose_HbA1c",   "TwoHourGlucose_HbA1c",  
"Hb_HbA1c"        ,       "MCV_HbA1c"            ,  "MCH_HbA1c"            , 
"MCHC_HbA1c"       ,      "Fe_HbA1c"              , "Ferritin_HbA1c"        ,
"TSAT_HbA1c",'sequenced',             'ascertainment_criteria')

whi_aa<-merge(whi_aa,linker[which(linker$study=='WHI'),],by.x='Individual_ID',by.y="study_subject_id")


colnames(whi_aa)[which(colnames(whi_aa)=="sample.id")]<-"TOPMEDID"




whi_aa$STUDY_TOPMEDID<-paste("WHI",whi_aa$TOPMEDID,sep="_")
whi_aa$STUDY_ANCESTRY<-"WHI_AA"

whi_aa<-whi_aa[,names(jhs)]



table(whi_aa$T2D,useNA="always")
table(whi_aa$T2D_FG,useNA="always")
table(whi_aa$T2D_FI,useNA="always")
table(whi_aa$T2D_HbA1c,useNA="always")

	whi_aa$T2D_FG<-ifelse(is.na(whi_aa$T2D_FG)==T,1,whi_aa$T2D_FG)
	whi_aa$T2D_FI<-ifelse(is.na(whi_aa$T2D_FI)==T,1,whi_aa$T2D_FI)
	whi_aa$T2D_HbA1c<-ifelse(is.na(whi_aa$T2D_HbA1c)==T,1,whi_aa$T2D_HbA1c)
	whi_aa$T2D<-ifelse(is.na(whi_aa$T2D)==T,0,whi_aa$T2D)
	


for (j in 1:nrow(whi_aa)){
	whi_aa$FastingGlucose[j]<-ifelse(whi_aa$T2D_FG[j]==2,NA,whi_aa$FastingGlucose[j])
	}
for (j in 1:nrow(whi_aa)){
	whi_aa$FastingInsulin[j]<-ifelse(whi_aa$T2D_FI[j]==2,NA,whi_aa$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
whi_aa$FastingGlucose<-ifelse(whi_aa$FastingGlucose>=7&whi_aa$T2D_FG==1,NA,whi_aa$FastingGlucose)
####set A1C>7 to NA
whi_aa$HbA1c<-ifelse(whi_aa$HbA1c>=6.5&whi_aa$T2D_HbA1c==1,NA,whi_aa$HbA1c)

whi<-rbind(whi_ha,whi_aa,whi_ea,whi_as)

############################## WHI ##############################
############################## WHI ##############################


############################## Genestar ##############################
############################## Genestar ##############################


####Genestar
##AA
names(gs_aa)[1:4]<-c("Family_ID","Individual_ID","Father_ID","Mother_ID")

gs_aa<-merge(gs_aa,linker[which(linker$study=='GeneSTAR'),],by.x='Individual_ID',by.y="study_subject_id")



colnames(gs_aa)[which(colnames(gs_aa)=="sample.id")]<-"TOPMEDID"

gs_aa$T2D_HbA1c<-1
gs_aa$age_HbA1c<-NA
gs_aa$BMI_HbA1c<-NA
gs_aa$FastingGlucose_HbA1c<-NA
gs_aa$TwoHourGlucose_HbA1c<-NA
gs_aa$Hb_HbA1c<-NA
gs_aa$MCV_HbA1c<-NA
gs_aa$MCH_HbA1c<-NA
gs_aa$MCHC_HbA1c<-NA
gs_aa$Fe_HbA1c<-NA
gs_aa$Ferritin_HbA1c<-NA
gs_aa$TSAT_HbA1c<-NA
gs_aa$ascertainment_criteria<-NA


gs_aa$STUDY_TOPMEDID<-paste("GeneSTAR",gs_aa$TOPMEDID,sep="_")
gs_aa$STUDY_ANCESTRY<-"GeneSTAR_AA"


gs_aa<-gs_aa[,names(jhs)]


table(gs_aa$T2D,useNA="always")
table(gs_aa$T2D_FG,useNA="always")
table(gs_aa$T2D_FI,useNA="always")
table(gs_aa$T2D_HbA1c,useNA="always")

	gs_aa$T2D_FG<-ifelse(is.na(gs_aa$T2D_FG)==T,1,gs_aa$T2D_FG)
	gs_aa$T2D_FI<-ifelse(is.na(gs_aa$T2D_FI)==T,1,gs_aa$T2D_FI)
	gs_aa$T2D_HbA1c<-ifelse(is.na(gs_aa$T2D_HbA1c)==T,1,gs_aa$T2D_HbA1c)
	gs_aa$T2D<-ifelse(is.na(gs_aa$T2D)==T,0,gs_aa$T2D)
	


for (j in 1:nrow(gs_aa)){
	gs_aa$FastingGlucose[j]<-ifelse(gs_aa$T2D_FG[j]==2,NA,gs_aa$FastingGlucose[j])
	}
for (j in 1:nrow(gs_aa)){
	gs_aa$FastingInsulin[j]<-ifelse(gs_aa$T2D_FI[j]==2,NA,gs_aa$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
gs_aa$FastingGlucose<-ifelse(gs_aa$FastingGlucose>=7&gs_aa$T2D_FG==1,NA,gs_aa$FastingGlucose)
####set A1C>7 to NA
gs_aa$HbA1c<-ifelse(gs_aa$HbA1c>=6.5&gs_aa$T2D_HbA1c==1,NA,gs_aa$HbA1c)

####EU

names(gs_ea)[1:4]<-c("Family_ID","Individual_ID","Father_ID","Mother_ID")

gs_ea<-merge(gs_ea,linker[which(linker$study=='GeneSTAR'),],by.x='Individual_ID',by.y="study_subject_id")



colnames(gs_ea)[which(colnames(gs_ea)=="sample.id")]<-"TOPMEDID"
colnames(gs_ea)[which(colnames(gs_ea)=="HbA1C")]<-"HbA1c"

gs_ea$T2D_HbA1c<-1
gs_ea$age_HbA1c<-NA
gs_ea$BMI_HbA1c<-NA
gs_ea$FastingGlucose_HbA1c<-NA
gs_ea$TwoHourGlucose_HbA1c<-NA
gs_ea$Hb_HbA1c<-NA
gs_ea$MCV_HbA1c<-NA
gs_ea$MCH_HbA1c<-NA
gs_ea$MCHC_HbA1c<-NA
gs_ea$Fe_HbA1c<-NA
gs_ea$Ferritin_HbA1c<-NA
gs_ea$TSAT_HbA1c<-NA
gs_ea$ascertainment_criteria<-NA


gs_ea$STUDY_TOPMEDID<-paste("GeneSTAR",gs_ea$TOPMEDID,sep="_")
gs_ea$STUDY_ANCESTRY<-"GeneSTAR_EU"


gs_ea<-gs_ea[,names(jhs)]


table(gs_ea$T2D,useNA="always")
table(gs_ea$T2D_FG,useNA="always")
table(gs_ea$T2D_FI,useNA="always")
table(gs_ea$T2D_HbA1c,useNA="always")

	gs_ea$T2D_FG<-ifelse(is.na(gs_ea$T2D_FG)==T,1,gs_ea$T2D_FG)
	gs_ea$T2D_FI<-ifelse(is.na(gs_ea$T2D_FI)==T,1,gs_ea$T2D_FI)
	gs_ea$T2D_HbA1c<-ifelse(is.na(gs_ea$T2D_HbA1c)==T,1,gs_ea$T2D_HbA1c)
	gs_ea$T2D<-ifelse(is.na(gs_ea$T2D)==T,0,gs_ea$T2D)
	


for (j in 1:nrow(gs_ea)){
	gs_ea$FastingGlucose[j]<-ifelse(gs_ea$T2D_FG[j]==2,NA,gs_ea$FastingGlucose[j])
	}
for (j in 1:nrow(gs_ea)){
	gs_ea$FastingInsulin[j]<-ifelse(gs_ea$T2D_FI[j]==2,NA,gs_ea$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
gs_ea$FastingGlucose<-ifelse(gs_ea$FastingGlucose>=7&gs_ea$T2D_FG==1,NA,gs_ea$FastingGlucose)
####set A1C>7 to NA
gs_ea$HbA1c<-ifelse(gs_ea$HbA1c>=6.5&gs_ea$T2D_HbA1c==1,NA,gs_ea$HbA1c)

gs<-rbind(gs_ea,gs_aa)


############################## Genestar ##############################
############################## Genestar ##############################


############################## CHS ##############################
############################## CHS ##############################

#####CHS

names(chs)[1:4]<-c("Family_ID","Individual_ID","Father_ID","Mother_ID")

chs<-merge(chs,linker[which(linker$study=='CHS'),],by.x='Individual_ID',by.y="study_subject_id")


chs$STUDY_ANCESTRY<-ifelse(chs$race01==1,paste("CHS","EU",sep="_"),ifelse(chs$race01==5,paste("CHS","Other",sep="_"),paste("CHS","AA",sep="_")))


chs<-chs[,!(names(chs)%in%c("race01","clinic"))]


colnames(chs)[which(colnames(chs)=="sample.id")]<-"TOPMEDID"

chs$HbA1c<-NA
chs$T2D_HbA1c<-1
chs$age_HbA1c<-NA
chs$BMI_HbA1c<-NA
chs$FastingGlucose_HbA1c<-NA
chs$TwoHourGlucose_HbA1c<-NA
chs$Hb_HbA1c<-NA
chs$MCV_HbA1c<-NA
chs$MCH_HbA1c<-NA
chs$MCHC_HbA1c<-NA
chs$Fe_HbA1c<-NA
chs$Ferritin_HbA1c<-NA
chs$TSAT_HbA1c<-NA
chs$ascertainment_criteria<-NA


chs$STUDY_TOPMEDID<-paste("CHS",chs$TOPMEDID,sep="_")


chs<-chs[,names(jhs)]


table(chs$T2D,useNA="always")
table(chs$T2D_FG,useNA="always")
table(chs$T2D_FI,useNA="always")
table(chs$T2D_HbA1c,useNA="always")

	chs$T2D_FG<-ifelse(is.na(chs$T2D_FG)==T,1,chs$T2D_FG)
	chs$T2D_FI<-ifelse(is.na(chs$T2D_FI)==T,1,chs$T2D_FI)
	chs$T2D_HbA1c<-ifelse(is.na(chs$T2D_HbA1c)==T,1,chs$T2D_HbA1c)
	chs$T2D<-ifelse(is.na(chs$T2D)==T,0,chs$T2D)
	


for (j in 1:nrow(chs)){
	chs$FastingGlucose[j]<-ifelse(chs$T2D_FG[j]==2,NA,chs$FastingGlucose[j])
	}
for (j in 1:nrow(chs)){
	chs$FastingInsulin[j]<-ifelse(chs$T2D_FI[j]==2,NA,chs$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
chs$FastingGlucose<-ifelse(chs$FastingGlucose>=7&chs$T2D_FG==1,NA,chs$FastingGlucose)
####set A1C>7 to NA
chs$HbA1c<-ifelse(chs$HbA1c>=6.5&chs$T2D_HbA1c==1,NA,chs$HbA1c)

############################## CHS ##############################
############################## CHS ##############################


############################## MESA ##############################
############################## MESA ##############################


####MESA
###ha
names(mesa_ha)[1]<-"Family_ID" 



mesa_ha<-merge(mesa_ha,linker[which(linker$study=='MESA'),],by.x='Individual_ID',by.y="study_subject_id")

colnames(mesa_ha)[which(colnames(mesa_ha)=="sample.id")]<-"TOPMEDID"
colnames(mesa_ha)[which(colnames(mesa_ha)=="Sex")]<-"sex"

mesa_ha$TwoHourGlucose_HbA1c<-NA
mesa_ha$Hb_HbA1c<-NA
mesa_ha$MCV_HbA1c<-NA
mesa_ha$MCH_HbA1c<-NA
mesa_ha$MCHC_HbA1c<-NA
mesa_ha$Fe_HbA1c<-NA
mesa_ha$Ferritin_HbA1c<-NA
mesa_ha$TSAT_HbA1c<-NA


mesa_ha$STUDY_TOPMEDID<-paste("MESA",mesa_ha$TOPMEDID,sep="_")
mesa_ha$STUDY_ANCESTRY<-"MESA_HA"


mesa_ha<-mesa_ha[,names(jhs)]



table(mesa_ha$T2D,useNA="always")
table(mesa_ha$T2D_FG,useNA="always")
table(mesa_ha$T2D_FI,useNA="always")
table(mesa_ha$T2D_HbA1c,useNA="always")

	mesa_ha$T2D_FG<-mesa_ha$T2D_FG+1
	mesa_ha$T2D_FI<-mesa_ha$T2D_FI+1
	mesa_ha$T2D_HbA1c<-mesa_ha$T2D_HbA1c+1

	


for (j in 1:nrow(mesa_ha)){
	mesa_ha$FastingGlucose[j]<-ifelse(mesa_ha$T2D_FG[j]==2,NA,mesa_ha$FastingGlucose[j])
	}
for (j in 1:nrow(mesa_ha)){
	mesa_ha$FastingInsulin[j]<-ifelse(mesa_ha$T2D_FI[j]==2,NA,mesa_ha$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
mesa_ha$FastingGlucose<-ifelse(mesa_ha$FastingGlucose>=7&mesa_ha$T2D_FG==1,NA,mesa_ha$FastingGlucose)
####set A1C>7 to NA
mesa_ha$HbA1c<-ifelse(mesa_ha$HbA1c>=6.5&mesa_ha$T2D_HbA1c==1,NA,mesa_ha$HbA1c)


###EA
names(mesa_ea)[1]<-"Family_ID" 

mesa_ea<-merge(mesa_ea,linker[which(linker$study=='MESA'),],by.x='Individual_ID',by.y="study_subject_id")

colnames(mesa_ea)[which(colnames(mesa_ea)=="sample.id")]<-"TOPMEDID"
colnames(mesa_ea)[which(colnames(mesa_ea)=="Sex")]<-"sex"

mesa_ea$TwoHourGlucose_HbA1c<-NA
mesa_ea$Hb_HbA1c<-NA
mesa_ea$MCV_HbA1c<-NA
mesa_ea$MCH_HbA1c<-NA
mesa_ea$MCHC_HbA1c<-NA
mesa_ea$Fe_HbA1c<-NA
mesa_ea$Ferritin_HbA1c<-NA
mesa_ea$TSAT_HbA1c<-NA


mesa_ea$STUDY_TOPMEDID<-paste("MESA",mesa_ea$TOPMEDID,sep="_")
mesa_ea$STUDY_ANCESTRY<-"MESA_EU"


mesa_ea<-mesa_ea[,names(jhs)]



table(mesa_ea$T2D,useNA="always")
table(mesa_ea$T2D_FG,useNA="always")
table(mesa_ea$T2D_FI,useNA="always")
table(mesa_ea$T2D_HbA1c,useNA="always")

	mesa_ea$T2D_FG<-mesa_ea$T2D_FG+1
	mesa_ea$T2D_FI<-mesa_ea$T2D_FI+1
	mesa_ea$T2D_HbA1c<-mesa_ea$T2D_HbA1c+1

	


for (j in 1:nrow(mesa_ea)){
	mesa_ea$FastingGlucose[j]<-ifelse(mesa_ea$T2D_FG[j]==2,NA,mesa_ea$FastingGlucose[j])
	}
for (j in 1:nrow(mesa_ea)){
	mesa_ea$FastingInsulin[j]<-ifelse(mesa_ea$T2D_FI[j]==2,NA,mesa_ea$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
mesa_ea$FastingGlucose<-ifelse(mesa_ea$FastingGlucose>=7&mesa_ea$T2D_FG==1,NA,mesa_ea$FastingGlucose)
####set A1C>7 to NA
mesa_ea$HbA1c<-ifelse(mesa_ea$HbA1c>=6.5&mesa_ea$T2D_HbA1c==1,NA,mesa_ea$HbA1c)


#####SA



names(mesa_sa)[1]<-"Family_ID" 




mesa_sa<-merge(mesa_sa,linker[which(linker$study=='MESA'),],by.x='Individual_ID',by.y="study_subject_id")

colnames(mesa_sa)[which(colnames(mesa_sa)=="sample.id")]<-"TOPMEDID"
colnames(mesa_sa)[which(colnames(mesa_sa)=="Sex")]<-"sex"

mesa_sa$TwoHourGlucose_HbA1c<-NA
mesa_sa$Hb_HbA1c<-NA
mesa_sa$MCV_HbA1c<-NA
mesa_sa$MCH_HbA1c<-NA
mesa_sa$MCHC_HbA1c<-NA
mesa_sa$Fe_HbA1c<-NA
mesa_sa$Ferritin_HbA1c<-NA
mesa_sa$TSAT_HbA1c<-NA


mesa_sa$STUDY_TOPMEDID<-paste("MESA",mesa_sa$TOPMEDID,sep="_")
mesa_sa$STUDY_ANCESTRY<-"MESA_AS"


mesa_sa<-mesa_sa[,names(jhs)]



table(mesa_sa$T2D,useNA="always")
table(mesa_sa$T2D_FG,useNA="always")
table(mesa_sa$T2D_FI,useNA="always")
table(mesa_sa$T2D_HbA1c,useNA="always")

	mesa_sa$T2D_FG<-mesa_sa$T2D_FG+1
	mesa_sa$T2D_FI<-mesa_sa$T2D_FI+1
	mesa_sa$T2D_HbA1c<-mesa_sa$T2D_HbA1c+1

	


for (j in 1:nrow(mesa_sa)){
	mesa_sa$FastingGlucose[j]<-ifelse(mesa_sa$T2D_FG[j]==2,NA,mesa_sa$FastingGlucose[j])
	}
for (j in 1:nrow(mesa_sa)){
	mesa_sa$FastingInsulin[j]<-ifelse(mesa_sa$T2D_FI[j]==2,NA,mesa_sa$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
mesa_sa$FastingGlucose<-ifelse(mesa_sa$FastingGlucose>=7&mesa_sa$T2D_FG==1,NA,mesa_sa$FastingGlucose)
####set A1C>7 to NA
mesa_sa$HbA1c<-ifelse(mesa_sa$HbA1c>=6.5&mesa_sa$T2D_HbA1c==1,NA,mesa_sa$HbA1c)

####AA
names(mesa_aa)[1]<-"Family_ID" 

mesa_aa<-merge(mesa_aa,linker[which(linker$study=='MESA'),],by.x='Individual_ID',by.y="study_subject_id")

colnames(mesa_aa)[which(colnames(mesa_aa)=="sample.id")]<-"TOPMEDID"
colnames(mesa_aa)[which(colnames(mesa_aa)=="Sex")]<-"sex"

mesa_aa$TwoHourGlucose_HbA1c<-NA
mesa_aa$Hb_HbA1c<-NA
mesa_aa$MCV_HbA1c<-NA
mesa_aa$MCH_HbA1c<-NA
mesa_aa$MCHC_HbA1c<-NA
mesa_aa$Fe_HbA1c<-NA
mesa_aa$Ferritin_HbA1c<-NA
mesa_aa$TSAT_HbA1c<-NA


mesa_aa$STUDY_TOPMEDID<-paste("MESA",mesa_aa$TOPMEDID,sep="_")
mesa_aa$STUDY_ANCESTRY<-"MESA_AA"


mesa_aa<-mesa_aa[,names(jhs)]



table(mesa_aa$T2D,useNA="always")
table(mesa_aa$T2D_FG,useNA="always")
table(mesa_aa$T2D_FI,useNA="always")
table(mesa_aa$T2D_HbA1c,useNA="always")

	mesa_aa$T2D_FG<-mesa_aa$T2D_FG+1
	mesa_aa$T2D_FI<-mesa_aa$T2D_FI+1
	mesa_aa$T2D_HbA1c<-mesa_aa$T2D_HbA1c+1

	


for (j in 1:nrow(mesa_aa)){
	mesa_aa$FastingGlucose[j]<-ifelse(mesa_aa$T2D_FG[j]==2,NA,mesa_aa$FastingGlucose[j])
	}
for (j in 1:nrow(mesa_aa)){
	mesa_aa$FastingInsulin[j]<-ifelse(mesa_aa$T2D_FI[j]==2,NA,mesa_aa$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
mesa_aa$FastingGlucose<-ifelse(mesa_aa$FastingGlucose>=7&mesa_aa$T2D_FG==1,NA,mesa_aa$FastingGlucose)
####set A1C>7 to NA
mesa_aa$HbA1c<-ifelse(mesa_aa$HbA1c>=6.5&mesa_aa$T2D_HbA1c==1,NA,mesa_aa$HbA1c)


mesa<-rbind(mesa_ha,mesa_aa,mesa_ea,mesa_sa)

############################## MESA ##############################
############################## MESA ##############################

############################## Genoa ##############################
############################## Genoa ##############################

p0<-rbind(fhs,jhs,sas,cfs,amish,gensalt,aric,whi,chs,mesa,gs)


p0$HbA1c<-ifelse(p0$T2D_HbA1c==2,NA,p0$HbA1c)


####genoa
names(genoa)[1:4]<-c("Family_ID","Individual_ID","Father_ID","Mother_ID")

genoa<-merge(genoa,linker[which(linker$study=='GENOA'),],by.x='Individual_ID',by.y="study_subject_id")

colnames(genoa)[which(colnames(genoa)=="sample.id")]<-"TOPMEDID"


genoa$HbA1c<-NA
genoa$T2D_HbA1c<-1
genoa$age_HbA1c<-NA
genoa$BMI_HbA1c<-NA
genoa$FastingGlucose_HbA1c<-NA
genoa$TwoHourGlucose_HbA1c<-NA
genoa$Hb_HbA1c<-NA
genoa$MCV_HbA1c<-NA
genoa$MCH_HbA1c<-NA
genoa$MCHC_HbA1c<-NA
genoa$Fe_HbA1c<-NA
genoa$Ferritin_HbA1c<-NA
genoa$TSAT_HbA1c<-NA

genoa$STUDY_TOPMEDID<-paste("GENOA",genoa$TOPMEDID,sep="_")
genoa$STUDY_ANCESTRY<-"GENOA_AA"

genoa<-genoa[,names(jhs)]

table(genoa$T2D,useNA="always")
table(genoa$T2D_FG,useNA="always")
table(genoa$T2D_FI,useNA="always")
table(genoa$T2D_HbA1c,useNA="always")

genoa$T2D_FG<-ifelse(genoa$T2D_FG==2,genoa$T2D_FG,1)
genoa$T2D_FI<-ifelse(genoa$T2D_FI==2,genoa$T2D_FI,1)
	
for (j in 1:nrow(genoa)){
	genoa$FastingGlucose[j]<-ifelse(genoa$T2D_FG[j]==2,NA,genoa$FastingGlucose[j])
	}
for (j in 1:nrow(genoa)){
	genoa$FastingInsulin[j]<-ifelse(genoa$T2D_FI[j]==2,NA,genoa$FastingInsulin[j])
	}	

### set FG to NA if FastingGlucose>=7&T2D_FG==1
genoa$FastingGlucose<-ifelse(genoa$FastingGlucose>=7&genoa$T2D_FG==1,NA,genoa$FastingGlucose)
####set A1C>7 to NA
genoa$HbA1c<-ifelse(genoa$HbA1c>=6.5&genoa$T2D_HbA1c==1,NA,genoa$HbA1c)

############################## Genoa ##############################
############################## Genoa ##############################


p0<-rbind(p0,genoa)

write.table(p0,"/restricted/projectnb/glycemic/peitao/phenotype_harmonization/pooled_analysis/Pooled_Glycemic_Traits_freeze5_duplicate_ID_20180116.ped",sep="\t",col=T,row=F,quote=FALSE)
