## T2D Harmonization  
## GENESIS needs sex as M,F
# FamilyID = SUBJECT_ID or individual_id depending on dataset
## PaternalID = 0 if not provided
## MaternalID = 0 if not provided
# create BMI variable
# BMI = 703·weight(lb)/height2(in2)
##

setwd("/N/dc2/scratch/wesselj/OriginalFiles/")

## MAP
map <- read.table("freeze5b_sample_annot_2017-12-01.txt",
                  header=TRUE, as.is=T, sep="\t") #n=54499 & 16 variables:
names(map)
# "sample.id"            "unique_subject_key"   "submitted_subject_id" "consent"             
# "sex"                  "sexchr.kary"          "topmed_phs"           "study"               
# "topmed_project"       "CENTER"               "geno.cntl"            "TRIO.dups"           
# "MZtwinID"             "keep"                 "unique.geno"          "unique.subj" 

# n=289 have study unique identifier=NA but have a NWD id assigned, dropping these for now
#map<-subset(map, (!is.na(map$unique_subject_key))) 

# per Josee Dupuis, keep trios DEC 2017 !! 

# n=54210
# Variables to filter on
#table(map$keep,useNA='always')
#FALSE  TRUE  <NA> 
#  375    54124 0 
#map <- map[map$keep == "TRUE",]
#n=54124 exclude n=375

table(map$study)
# Amish        ARIC        BAGS        CCAF         CFS         CHS    COPDGene         CRA         DHS 
# 1028        3602         962         329         920          69        8640        1041         337 
# EOCOPD         FHS      GALAII    GeneSTAR       GENOA     GenSalt       GOLDN         HVH    HyperGEN 
# 66        3749         913        1634        1141        1684         902          66        1774 
# JHS    Mayo_VTE        MESA      MGH_AF    Partners        SAFS        SAGE Sarcoidosis         SAS 
# 3128        1250        4814         918         111        1502         451         601        1208 
# VAFAR       VU_AF        WGHS         WHI 
# 154        1018          98       10024 

## Create variables
map$individual_id = map$submitted_subject_id
map$topmedid = map$sample.id
map$study_topmedid <- paste(map$study,map$sample.id, sep = "_")

map <- subset(map, select=c("study_topmedid","topmedid","individual_id","sample.id", "unique_subject_key",
                            "submitted_subject_id", "consent", "sexchr.kary",  "topmed_phs", "study", 
                            "topmed_project", "CENTER", "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", 
                            "unique.geno",  "unique.subj" ))

##
# AFib MGH
## new file downloaded from dbGaP 5AUG2017
# 7DEC2017 checked for new files - none

#names(raw.MGHPed)
# "FAMILY_ID"  "SUBJECT_ID" "FATHER"     "MOTHER"     "SEX" 
raw.MGHPed <- raw.MGHPed[,c("FAMILY_ID",  "SUBJECT_ID", "FATHER", "MOTHER")]
raw.MGH <- merge(raw.MGH,raw.MGHPed,by.x='SUBJECT_ID',by.y='SUBJECT_ID',all.x=TRUE) #n=1025 & 25 variables
rm(raw.MGHPed)
raw.MGH$JWsource = "dbGaP"
mapMGH <- map[which(map$study == "MGH_AF"),] #n=918
raw.MGH <- merge(mapMGH,raw.MGH,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=918
rm(mapMGH)

# recode & check variable names & distributions
table(raw.MGH$diabetes,useNA='always')
raw.MGH$t2d = ifelse(is.na(raw.MGH$diabetes),NA,raw.MGH$diabetes)
raw.MGH$t2d[raw.MGH$t2d == 'yes'] = 2
raw.MGH$t2d[raw.MGH$t2d == 'no'] = 0
with(raw.MGH,table(diabetes,t2d,useNA='always')) 
table(raw.MGH$race,useNA='always')
table(raw.MGH$ethnicity,useNA='always')
raw.MGH$ancestry[raw.MGH$race == 'american_indian' & raw.MGH$ethnicity == 'no'] = "AMR"
raw.MGH$ancestry[raw.MGH$race == 'american_indian_black' & raw.MGH$ethnicity == 'no'] = "MIXED"
raw.MGH$ancestry[raw.MGH$race == 'white' & raw.MGH$ethnicity == 'no'] = "EU"
raw.MGH$ancestry[raw.MGH$race == 'white' & raw.MGH$ethnicity == 'yes'] = "HS"
table(raw.MGH$ancestry,useNA='always')
table(raw.MGH$sex,useNA='always')
raw.MGH$origsex = raw.MGH$sex
raw.MGH$sex[raw.MGH$sex == 'male'] = 'M'
raw.MGH$sex[raw.MGH$sex == 'female'] = 'F'
with(raw.MGH,table(origsex,sex,useNA='always'))
summary(raw.MGH$age)
raw.MGH$last_exam_age = raw.MGH$age
#raw.MGH = subset(raw.MGH, age >= 25) #drop 11 individuals
summary(raw.MGH$height)
summary(raw.MGH$weight)
raw.MGH$last_exam_bmi = ((703*raw.MGH$weight)/(raw.MGH$height*raw.MGH$height))
summary(raw.MGH$last_exam_bmi) # low BMI - exclude?
raw.MGH$last_exam_fg = NA
raw.MGH$last_exam_hba1c = NA
raw.MGH$last_exam_t2d_treatment = NA
raw.MGH$last_exam_visit = NA
raw.MGH$t2d_age = NA
raw.MGH$t2d_bmi = NA
raw.MGH$FamilyID = raw.MGH$FAMILY_ID
raw.MGH$PaternalID = raw.MGH$FATHER
raw.MGH$MaternalID = raw.MGH$MOTHER
raw.MGH$study_ancestry <- paste(raw.MGH$study,raw.MGH$ancestry, sep = "_")

raw.MGH <- raw.MGH[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                      'study_topmedid','study_ancestry','JWsource','ancestry',
                      "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                      "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                      "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=918
#6 missing diabetes status

#
## AFib VU
# downloaded new dataset 7DEC2017

raw.VUdawPed <- raw.VUdawPed[,c("FAMILY_ID",  "SUBJECT_ID", "FATHER", "MOTHER")]
raw.VUdaw <- merge(raw.VUdaw,raw.VUdawPed,by.x='SUBJECT_ID',by.y='SUBJECT_ID',all.x=TRUE) #n=1134
rm(raw.VUdawPed)
raw.VUdaw$JWsource = "dbGaP"
mapV <- map[which(map$study == "VU_AF"),] #n=1018
raw.VUdaw <- merge(mapV,raw.VUdaw,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE)
rm(mapV)

# recode & check variable names & distributions
table(raw.VUdaw$diabetes,useNA='always')
raw.VUdaw$t2d = ifelse(is.na(raw.VUdaw$diabetes),NA,raw.VUdaw$diabetes)
raw.VUdaw$t2d[raw.VUdaw$t2d == 'yes'] = 2
raw.VUdaw$t2d[raw.VUdaw$t2d == 'no'] = 0
with(raw.VUdaw,table(diabetes,t2d,useNA='always'))
table(raw.VUdaw$race,useNA='always')
table(raw.VUdaw$ethnicity,useNA='always')
raw.VUdaw$ancestry[raw.VUdaw$race == 'native_american' & raw.VUdaw$ethnicity == 'no'] = "AMR"
raw.VUdaw$ancestry[raw.VUdaw$race == 'black' & raw.VUdaw$ethnicity == 'no'] = "AF"
raw.VUdaw$ancestry[raw.VUdaw$race == 'white' & raw.VUdaw$ethnicity == 'no'] = "EU"
raw.VUdaw$ancestry[raw.VUdaw$race == 'asian' & raw.VUdaw$ethnicity == 'no'] = "AS"
raw.VUdaw$ancestry[raw.VUdaw$race == 'white' & raw.VUdaw$ethnicity == 'yes'] = "HS"
table(raw.VUdaw$ancestry,useNA='always')
table(raw.VUdaw$sex,useNA='always')
raw.VUdaw$origsex = raw.VUdaw$sex
raw.VUdaw$sex[raw.VUdaw$sex == 'male'] = 'M'
raw.VUdaw$sex[raw.VUdaw$sex == 'female'] = 'F'
with(raw.VUdaw,table(origsex,sex,useNA='always'))
summary(raw.VUdaw$age)  ## Contains some < 18 yrs
raw.VUdaw$last_exam_age = raw.VUdaw$age
#raw.VUdaw = subset(raw.VUdaw, age >= 25) #n=1085, drop 24 individuals
summary(raw.VUdaw$height)
summary(raw.VUdaw$weight)
raw.VUdaw$last_exam_bmi = ((703*raw.VUdaw$weight)/(raw.VUdaw$height*raw.VUdaw$height))
summary(raw.VUdaw$last_exam_bmi) # low BMI - exclude?
raw.VUdaw$last_exam_fg = NA
raw.VUdaw$last_exam_hba1c = NA
raw.VUdaw$last_exam_t2d_treatment = NA
raw.VUdaw$last_exam_visit = NA
raw.VUdaw$t2d_age = NA
raw.VUdaw$t2d_bmi = NA
raw.VUdaw$FamilyID = raw.VUdaw$FAMILY_ID
raw.VUdaw$PaternalID = raw.VUdaw$FATHER
raw.VUdaw$MaternalID = raw.VUdaw$MOTHER
raw.VUdaw$study_ancestry <- paste(raw.VUdaw$study,raw.VUdaw$ancestry, sep = "_")

raw.VUdaw <- raw.VUdaw[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                          'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                          'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                          'study_topmedid','study_ancestry','JWsource','ancestry',
                          "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                          "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                          "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=1018

#
## AFib Partners
# 7DEC2017 using downloaded newly released data 6OCT2017

afp$JWsource = "dbGaP"
mapP <- map[which(map$study == "Partners"),] #n=111
afp <- merge(mapP,afp,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE)
rm(mapP)

# recode & check variable names & distributions
table(afp$diabetes,useNA='always')
afp$t2d = ifelse(is.na(afp$diabetes),NA,afp$diabetes)
afp$t2d[afp$diabetes == 'yes'] = 2
afp$t2d[afp$diabetes == 'no'] = 0
with(afp,table(diabetes,t2d,useNA='always'))
table(afp$race,useNA='always')
table(afp$ethnicity,useNA='always')
afp$ancestry[afp$race == 'black' & afp$ethnicity == 'no'] = "AF"
afp$ancestry[afp$race == 'white' & afp$ethnicity == 'no'] = "EU"
afp$ancestry[afp$race == 'white'] = "EU"
afp$ancestry[afp$race == 'asian' & afp$ethnicity == 'no'] = "AS"
afp$ancestry[afp$race == 'hispanic' & afp$ethnicity == 'yes'] = "HS"
afp$ancestry[afp$race == 'other' & afp$ethnicity == 'no'] = "OTHER"
table(afp$ancestry,useNA='always')
table(afp$sex,useNA='always')
afp$origsex = afp$sex
afp$sex[afp$sex == 'male'] = 'M'
afp$sex[afp$sex == 'female'] = 'F'
with(afp,table(origsex,sex,useNA='always'))
summary(afp$age)
afp$last_exam_age = afp$age
#afp = subset(afp, last_exam_age >= 25) #drop 2 individuals
afp = subset(afp, weight <= 401) #drop 1 individual
## exclude 1 individual with weight = 400+ ##
afp$weight <- as.numeric(afp$weight)
class(afp$weight)
summary(afp$height)
summary(afp$weight)
afp$last_exam_bmi = ((703*afp$weight)/(afp$height*afp$height))
summary(afp$last_exam_bmi)
afp$last_exam_fg = NA
afp$last_exam_hba1c = NA
afp$last_exam_t2d_treatment = NA
afp$last_exam_visit = NA
afp$t2d_age = NA
afp$t2d_bmi = NA
afp$FamilyID = afp$individual_id
afp$PaternalID = 0
afp$MaternalID = 0
afp$study_ancestry <- paste(afp$study,afp$ancestry, sep = "_")

afp <- afp[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=111


#
## AFib VU_Ben
# downlaoded new dataset 5AUG2017
# 7DEC2017 checked for new files - none

afvub$JWsource = "dbGaP"
mapB <- map[which(map$study == "VAFAR"),] #n=154
afvub <- merge(mapB,afvub,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=
rm(mapB)

# recode & check variable names & distributions
table(afvub$diabetes,useNA='always') 
afvub$t2d = ifelse(is.na(afvub$diabetes),NA,afvub$diabetes)
afvub$t2d[afvub$diabetes == 'yes'] = 2
afvub$t2d[afvub$diabetes == 'no'] = 0
with(afvub,table(diabetes,t2d,useNA='always'))
table(afvub$race,useNA='always')
table(afvub$ethnicity,useNA='always')
afvub$ancestry[afvub$race == 'white' & afvub$ethnicity == 'no'] = "EU"
table(afvub$ancestry,useNA='always')
table(afvub$sex,useNA='always')
afvub$origsex = afvub$sex
afvub$sex[afvub$sex == 'male'] = 'M'
afvub$sex[afvub$sex == 'female'] = 'F'
with(afvub,table(origsex,sex,useNA='always'))
afvub$last_exam_bmi = ((703*afvub$weight)/(afvub$height*afvub$height))
summary(afvub$last_exam_bmi)
summary(afvub$age)
afvub$last_exam_age = afvub$age
#afvub = subset(afvub, age >= 25) #n=117, no drops
summary(afvub$height)
summary(afvub$weight)
afvub$last_exam_bmi = ((703*afvub$weight)/(afvub$height*afvub$height))
summary(afvub$last_exam_bmi)
afvub$last_exam_fg = NA
afvub$last_exam_hba1c = NA
afvub$last_exam_t2d_treatment = NA
afvub$last_exam_visit = NA
afvub$t2d_age = NA
afvub$t2d_bmi = NA
afvub$FamilyID = afvub$individual_id
afvub$PaternalID = 0
afvub$MaternalID = 0
afvub$study_ancestry <- paste(afvub$study,afvub$ancestry, sep = "_")

afvub <- afvub[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                  'study_topmedid','study_ancestry','JWsource','ancestry',
                  "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                  "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                  "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=154

#
## AFib CCAF
# downlaoded new dataset 5AUG2017
# 7DEC2017 checked for new files - none

afccaf$JWsource = "dbGaP"
mapCCAF <- map[which(map$study == "CCAF"),] #n=329
afccaf <- merge(mapCCAF,afccaf,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE)
rm(mapCCAF) #n=329

# recode & check variable names & distributions
table(afccaf$diabetes,useNA='always')
afccaf$t2d = ifelse(is.na(afccaf$diabetes),NA,afccaf$diabetes)
afccaf$t2d[afccaf$diabetes == 'yes'] = 2
afccaf$t2d[afccaf$diabetes == 'no'] = 0
with(afccaf,table(diabetes,t2d,useNA='always')) #n=20 unknown t2d
table(afccaf$race,useNA='always')
table(afccaf$ethnicity,useNA='always')
afccaf$ancestry[afccaf$race == 'white' & afccaf$ethnicity == 'no'] = "EU"
afccaf$ancestry[afccaf$race == 'white' & afccaf$ethnicity == 'yes'] = "HS"
table(afccaf$ancestry,useNA='always')
table(afccaf$sex,useNA='always')
afccaf$origsex = afccaf$sex
afccaf$sex[afccaf$sex == 'male'] = 'M'
afccaf$sex[afccaf$sex == 'female'] = 'F'
with(afccaf,table(origsex,sex,useNA='always'))
summary(afccaf$age)
afccaf$last_exam_age = afccaf$age
#afccaf = subset(afccaf, age >= 25) #n=352, lose 3 individuals
summary(afccaf$height)
summary(afccaf$weight)
afccaf$last_exam_bmi = ((703*afccaf$weight)/(afccaf$height*afccaf$height))
summary(afccaf$last_exam_bmi) # ! low BMI
afccaf$last_exam_fg = NA
afccaf$last_exam_hba1c = NA
afccaf$last_exam_t2d_treatment = NA
afccaf$last_exam_visit = NA
afccaf$t2d_age = NA
afccaf$t2d_bmi = NA 
afccaf$FamilyID = afccaf$individual_id
afccaf$PaternalID = 0
afccaf$MaternalID = 0
afccaf$study_ancestry <- paste(afccaf$study,afccaf$ancestry, sep = "_")

afccaf <- afccaf[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                    'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                    'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                    'study_topmedid','study_ancestry','JWsource','ancestry',
                    "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                    "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                    "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=329

## HVH part of AFib and VTE projects
# HVH updated 18AUG2017

raw.HVH$JWsource = "dbGaP_Ex"
mapHVH <- map[which(map$study == "HVH"),] #n=66
raw.HVH <- merge(mapHVH,raw.HVH,by.x='individual_id',by.y='SampleID',all.x=TRUE)
rm(mapHVH)
#n=66

# recode & check variable names & distributions
table(raw.HVH$sequenced,useNA='always')
table(raw.HVH$T2D,useNA='always')
raw.HVH$t2d = raw.HVH$T2D
raw.HVH$t2d = ifelse(is.na(raw.HVH$T2D),NA,raw.HVH$T2D)
with(raw.HVH,table(T2D,t2d,useNA='always')) #
table(raw.HVH$RaceCalc,useNA='always')
raw.HVH$ancestry = raw.HVH$RaceCalc
raw.HVH$ancestry[raw.HVH$ancestry == 0] = 'AF'
raw.HVH$ancestry[raw.HVH$ancestry == 1] = 'EU'
raw.HVH$ancestry[raw.HVH$ancestry == 2] = 'OTHER'
table(raw.HVH$ancestry,useNA='always')
table(raw.HVH$sex,useNA='always')
raw.HVH$origsex = raw.HVH$sex
raw.HVH$sex[raw.HVH$sex == 1] = 'M'
raw.HVH$sex[raw.HVH$sex == 2] = 'F'
with(raw.HVH,table(origsex,sex,useNA='always'))
summary(raw.HVH$last_exam_age)
raw.HVH$last_exam_bmi = raw.HVH$last_exam_BMI
summary(raw.HVH$last_exam_bmi)
raw.HVH$last_exam_fg = NA
raw.HVH$last_exam_hba1c = NA
raw.HVH$last_exam_t2d_treatment = NA
raw.HVH$last_exam_visit = NA
raw.HVH$t2d_age = raw.HVH$T2D_age
summary(raw.HVH$T2D_age)
raw.HVH$t2d_bmi = NA
raw.HVH$study_ancestry <- paste(raw.HVH$study,raw.HVH$ancestry, sep = "_")

raw.HVH <- raw.HVH[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                      'study_topmedid','study_ancestry','JWsource','ancestry',
                      "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                      "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                      "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=66



## ARIC
# NEED TO FIX CODING OF last exam t2d treatment (Y/N/U) if used in analyses

#table(aric$T2D,useNA='always') 
aric$JWsource = "dbGaP_Ex"
mapAr <- map[which(map$study == "ARIC"),] #n=3602
aric <- merge(mapAr,aric,by.x='individual_id',by.y='Individual_ID',all.x=TRUE)
rm(mapAr)
rm(aric_ea)
rm(aric_aa)

# recode & check variable names & distributions
table(aric$T2D,useNA='always')
aric$t2d = ifelse(is.na(aric$T2D),-9,aric$T2D)
aric$t2d[aric$T2D == 2] = 2
aric$t2d[aric$T2D == 1] = 1
aric$t2d[aric$T2D == 0] = 0
with(aric,table(T2D,t2d,useNA='always')) #n=37 -9 & n=14 NA
table(aric$ancestry,useNA='always')
table(aric$sex,useNA='always')
aric$origsex = aric$sex
aric$sex[aric$sex == 1] = 'M'
aric$sex[aric$sex == 2] = 'F'
with(aric,table(origsex,sex,useNA='always'))
summary(aric$last_exam_age)
aric$last_exam_bmi = aric$last_exam_BMI
summary(aric$last_exam_bmi) #!low BMI
aric$last_exam_fg = NA
aric$last_exam_t2d_treatment = aric$last_exam_T2D_treatment
aric$last_exam_hba1c = aric$last_exam_HbA1c
aric$t2d_age = aric$T2D_age
summary(aric$T2D_age)
aric$t2d_bmi = aric$T2D_BMI
summary(aric$last_exam_hba1c)
summary(aric$t2d_age)
summary(aric$t2d_bmi)
table(aric$last_exam_visit)
table(aric$last_exam_t2d_treatment)
aric$FamilyID = aric$Family_ID
aric$PaternalID = aric$Father_ID
aric$MaternalID = aric$Mother_ID
aric$study_ancestry <- paste(aric$study,aric$ancestry, sep = "_")

aric <- aric[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                'study_topmedid','study_ancestry','JWsource','ancestry',
                "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=3602
# n=14 missing phenotype data

#####################################################

##COPD Gene Sample & C1 & C2

table(copd$Diabetes,useNA='always') #n=0 NA
rm(copd.c1)
rm(copd.c2)
copd$JWsource = "dbGaP"
mapS <- map[which(map$study == "COPDGene"),] ##n=8640
copd <- merge(mapS,copd,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE)
rm(mapS)

# recode & check variable names & distributions
table(copd$Diabetes,useNA='always') #n=102 NA
copd$t2d = ifelse(is.na(copd$Diabetes),-9,copd$Diabetes)
copd$t2d[copd$t2d == 1] = 2
copd$t2d[copd$t2d == 0] = 0
with(copd,table(Diabetes,t2d,useNA='always'))
table(copd$race,useNA='always')
table(copd$ethnic,useNA='always')
copd$ancestry[copd$race == 2 & copd$ethnic == 2] = "AF"
copd$ancestry[copd$race == 1 & copd$ethnic == 2] = "EU"
table(copd$ancestry,useNA='always')
table(copd$gender,useNA='always')
copd$sex = copd$gender
copd$sex[copd$sex == 1] = 'M'
copd$sex[copd$sex == 2] = 'F'
with(copd,table(gender,sex,useNA='always'))
summary(copd$Age_Enroll)
copd$last_exam_age = copd$Age_Enroll
copd$last_exam_bmi = copd$BMI
summary(copd$last_exam_bmi) # low BMI - exclude?
copd$last_exam_fg = NA
copd$last_exam_hba1c = NA
copd$last_exam_t2d_treatment = NA
copd$last_exam_visit = NA
copd$t2d_age = NA
copd$t2d_bmi = NA
copd$FamilyID = copd$individual_id
copd$PaternalID = 0
copd$MaternalID = 0
copd$study_ancestry <- paste(copd$study,copd$ancestry, sep = "_")

copd <- copd[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                'study_topmedid','study_ancestry','JWsource','ancestry',
                "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=8742
#n=102 with missing phenotype data

## Amish
## Amish has 2 pairs of twins. they left the co-twin out "were removed on our end as
#twins caused some issues with some of our phenotypes in our in house analysis (but not lipids)"
# 15AUG2017

#Diabetes_Status NA=0
amish$JWsource = "dbGaP_Ex"
mapAm <- map[which(map$study == "Amish"),] #n=1028
amish <- merge(mapAm,amish,by.x='individual_id',by.y='SUBJID',all.x=TRUE) #n=1028
rm(mapAm)

# recode & check variable names & distributions
table(amish$sequenced,useNA='always')
table(amish$Diabetes_Status,useNA='always') # NA=92
amish$t2d = ifelse(is.na(amish$Diabetes_Status),-9,amish$Diabetes_Status)
amish$t2d[amish$Diabetes_Status == 1] = 2
amish$t2d[amish$Diabetes_Status == 0] = 0
with(amish,table(Diabetes_Status,t2d,useNA='always'))
amish$ancestry = 'EU'
table(amish$ancestry,useNA='always')
table(amish$sex,useNA='always')
amish$origsex = amish$sex
amish$sex[amish$sex == 1] = 'M'
amish$sex[amish$sex == 2] = 'F'
with(amish,table(origsex,sex,useNA='always'))
summary(amish$last_exam_age)
#amish = subset(amish, last_exam_age >= 25) #n=962
summary(amish$last_exam_BMI)
summary(amish$T2D_BMI)
summary(amish$T2D_age)
amish$last_exam_bmi = amish$last_exam_BMI
amish$last_exam_fg = amish$last_exam_FG
amish$last_exam_hba1c = amish$last_exam_HbA1c
amish$last_exam_t2d_treatment = amish$last_exam_T2D_treatment
amish$t2d_age = amish$T2D_age
amish$t2d_bmi = amish$T2D_BMI
amish$last_exam_visit = NA
amish$FamilyID = amish$FAMID
amish$PaternalID = amish$FATHER
amish$MaternalID = amish$MOTHER
amish$study_ancestry <- paste(amish$study,amish$ancestry, sep = "_")
amish <- amish[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                  'study_topmedid','study_ancestry','JWsource','ancestry',
                  "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                  "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                  "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=1030

#####################################################
## CFS
# imported data 5AUG2017
# CFS has 1 MZ twin and both are included in pheno, although they end up being excluded because of age
# 15AUG2017

#T2D n=880 NA
cfs$JWsource = "dbGaP_Ex"
mapCFS <- map[which(map$study == "CFS"),] #n=920
cfs <- merge(mapCFS,cfs,by.x='individual_id',by.y='IID',all.x=TRUE) #n=920
rm(mapCFS)

# recode & check variable names & distributions
table(cfs$Sequenced,useNA='always')
table(cfs$T2D,useNA='always') # n=253 unknown
cfs$t2d = ifelse(is.na(cfs$T2D),-9,cfs$T2D)
cfs$t2d[cfs$T2D == 2] = 2
cfs$t2d[cfs$T2D == 1] = 1
cfs$t2d[cfs$T2D == 0] = 0
with(cfs,table(T2D,t2d,useNA='always')) #n=253 with unknown t2d
table(cfs$Population,useNA='always')
cfs$ancestry[cfs$Population == 'CFS-blacks'] = "AF"
cfs$ancestry[cfs$Population == 'CFS-whites'] = "EU"
with(cfs,table(Population,ancestry,useNA='always'))
table(cfs$sex,useNA='always')
cfs$origsex = cfs$sex
cfs$sex[cfs$sex == 1] = 'M'
cfs$sex[cfs$sex == 2] = 'F'
with(cfs,table(origsex,sex,useNA='always'))
summary(cfs$last_exam_age)
#cfs = subset(cfs, last_exam_age >= 25) #n=727 drop 267
summary(cfs$T2D_age)
#cfs = cfs[(cfs$T2D_age >= 25 | is.na(cfs$T2D_age)),] #n=723 drop 4
summary(cfs$last_exam_BMI)
cfs$last_exam_bmi = cfs$last_exam_BMI
cfs$last_exam_fg = cfs$FastingGlucose
cfs$last_exam_hba1c = NA
cfs$last_exam_t2d_treatment = cfs$last_exam_T2D_treatment
cfs$t2d_age = cfs$T2D_age
cfs$t2d_bmi = cfs$T2D_BMI
cfs$last_exam_visit = NA
cfs$FamilyID = cfs$FamID
cfs$PaternalID = cfs$FID
cfs$MaternalID = cfs$MID
cfs$study_ancestry <- paste(cfs$study,cfs$ancestry, sep = "_")
cfs <- cfs[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=920

#####################################################
## SAS

#T2D NA=0 & -9=497
sas$JWsource = "dbGaP_Ex"
mapSAS <- map[which(map$study == "SAS"),] #n=1208
sas <- merge(mapSAS,sas,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=1208
rm(mapSAS)


#sas2 <- sas[is.na(sas$T2D),]

#n=8 with WGS 7 no phenotype

# recode & check variable names & distributions
table(sas$sequenced,useNA='always')
table(sas$T2D,useNA='always')
sas$t2d = ifelse(is.na(sas$T2D),-9,sas$T2D)
sas$t2d[sas$T2D == 2] = 2
sas$t2d[sas$T2D == 1] = 1
sas$t2d[sas$T2D == 0] = 0
with(sas,table(T2D,t2d,useNA='always'))
sas$ancestry = 'AS'
table(sas$ancestry,useNA='always')
table(sas$sex,useNA='always')
sas$origsex = sas$sex
sas$sex[sas$sex == 1] = 'M'
sas$sex[sas$sex == 2] = 'F'
with(sas,table(origsex,sex,useNA='always'))
summary(sas$last_exam_age)
summary(sas$last_exam_BMI)
sas$last_exam_bmi = sas$last_exam_BMI
summary(sas$last_exam_BMI)
sas$last_exam_fg = sas$last_exam_FG
summary(sas$last_exam_fg)
sas$last_exam_hba1c = NA
sas$last_exam_t2d_treatment = sas$last_exam_T2D_treatment
sas$t2d_age = sas$T2D_age
summary(sas$T2D_age)
sas$t2d_bmi = NA
table(sas$last_exam_visit,useNA = 'always')
sas$FamilyID = sas$Family_ID
sas$PaternalID = sas$Father_ID
sas$MaternalID = sas$Mother_ID
sas$study_ancestry <- paste(sas$study,sas$ancestry, sep = "_")

sas <- sas[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=1208


###############################################
##JHS
# 15AUG2017 JHS has 2 MZ twins, currently co-twins are being excluded
# 22AUG2017 now all co-twins are being inlcuded

jhs$JWsource = "dbGaP_Ex"
mapjhs <- map[which(map$study == "JHS"),] #n=3128
jhs <- merge(mapjhs,jhs,by.x='individual_id',by.y='Individual_ID',all.x=TRUE) #n=3128
rm(mapjhs)

# recode & check variable names & distributions
table(jhs$T2D,useNA='always')
jhs$t2d = ifelse(is.na(jhs$T2D),-9,jhs$T2D)
with(jhs,table(T2D,t2d,useNA='always'))
jhs$ancestry = 'AF'
table(jhs$ancestry,useNA='always')
table(jhs$sex,useNA='always')
jhs$origsex = jhs$sex
jhs$sex[jhs$sex == 1] = 'M'
jhs$sex[jhs$sex == 2] = 'F'
with(jhs,table(origsex,sex,useNA='always'))
jhs$last_exam_age[jhs$last_exam_age == 'x'] = 'NA'
jhs$last_exam_age <- as.numeric(jhs$last_exam_age)
summary(jhs$last_exam_age)
jhs$last_exam_BMI[jhs$last_exam_BMI == 'x'] = 'NA'
jhs$last_exam_bmi <- as.numeric(jhs$last_exam_BMI)
summary(jhs$last_exam_bmi)
jhs$last_exam_FG[jhs$last_exam_FG == 'x'] = 'NA'
jhs$last_exam_fg <- as.numeric(jhs$last_exam_FG)
summary(jhs$last_exam_fg)
jhs$last_exam_HbA1c[jhs$last_exam_HbA1c == 'x'] = 'NA'
jhs$last_exam_hba1c <- as.numeric(jhs$last_exam_HbA1c)
summary(jhs$last_exam_hba1c)
jhs$last_exam_T2D_treatment[jhs$last_exam_T2D_treatment == 'x'] = 'NA'
jhs$last_exam_t2d_treatment <- as.numeric(jhs$last_exam_T2D_treatment)
table(jhs$last_exam_t2d_treatment)
jhs$T2D_age[jhs$T2D_age == 'x'] = 'NA'
jhs$t2d_age <- as.numeric(jhs$T2D_age)
summary(jhs$t2d_age)
jhs$T2D_BMI[jhs$T2D_BMI == 'x'] = 'NA'
jhs$t2d_bmi <- as.numeric(jhs$T2D_BMI)
summary(jhs$t2d_bmi)
jhs$last_exam_visit = NA
jhs$FamilyID = jhs$Family_ID
jhs$PaternalID = jhs$Father_ID
jhs$MaternalID = jhs$Mother_ID
jhs$study_ancestry <- paste(jhs$study,jhs$ancestry, sep = "_")

jhs <- jhs[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=3128

###############################################
##FHS
# #13AUG2017 currently n=4024 since 1 individual has not been WGSed
# 15AUG2017 FHS has 6 MZ twin pairs. 
#
###############################################
## harmonized dataset downloaded from dbGaP Ex

#NA = 0
fhs$JWsource = "dbGaP_Ex"
mapfhs <- map[which(map$study == "FHS"),] #n=3749
fhs <- merge(mapfhs,fhs,by.x='individual_id',by.y='shareid',all.x=TRUE) #n=3749
rm(mapfhs)

fhs2 <- fhs[is.na(fhs$T2D),]
#n=107 with WGS & no phenotype

# recode & check variable names & distributions
table(fhs$T2D,useNA='always') #NA=107 
#fhs$t2d = ifelse(is.na(fhs$T2D),-9,fhs$T2D)
fhs$t2d[fhs$T2D == 2] = 2
fhs$t2d[fhs$T2D == 1] = 1
fhs$t2d[fhs$T2D == 0] = 0
with(fhs,table(T2D,t2d,useNA='always'))
table(fhs$sex,useNA='always')
fhs$origsex = fhs$sex
fhs$sex[fhs$sex == 1] = 'M'
fhs$sex[fhs$sex == 2] = 'F'
with(fhs,table(sex,origsex,useNA='always'))
summary(fhs$last_exam_age)
fhs$last_exam_bmi = fhs$last_exam_BMI
summary(fhs$last_exam_bmi)
fhs$last_exam_fg = fhs$last_exam_FG
summary(fhs$last_exam_fg)
fhs$last_exam_hba1c = fhs$last_exam_HbA1c
summary(fhs$last_exam_hba1c)
fhs$last_exam_t2d_treatment = fhs$last_exam_T2D_treatment
table(fhs$last_exam_t2d_treatment,useNA='always')
fhs$t2d_age = fhs$T2D_age
summary(fhs$t2d_age)
fhs$t2d_bmi = fhs$T2D_BMI
summary(fhs$t2d_bmi)
table(fhs$last_exam_visit,useNA = 'always')
fhs$FamilyID = fhs$individual_id
fhs$PaternalID = 0
fhs$MaternalID = 0
fhs$ancestry = 'EU'
table(fhs$ancestry,useNA='always')
table(fhs$idtype,useNA = 'always')
fhs$studygen = ifelse(is.na(fhs$idtype),"fhsNA",fhs$idtype)
copd$ancestry[copd$race == 2 & copd$ethnic == 2] = "AF"
fhs$studygen[fhs$idtype == 0] = "fhs0"
fhs$studygen[fhs$idtype == 1] = "fhs1"
fhs$studygen[fhs$idtype == 3] = "fhs3"
table(fhs$studygen,useNA = 'always')
fhs$study_ancestry <- paste(fhs$studygen,fhs$ancestry, sep = "_")
table(fhs$study_ancestry,useNA = 'always')

fhs <- fhs[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]

#n=3758 but 117 are missing phenotype data????
table(fhs$t2d,useNA='always')


## WHI 

table(whi$t2d,useNA='always') #Na=0
table(whi$sex,useNA='always') #Na=0
table(whi$ancestry,useNA='always') #Na=0

rm(whi.eu)
rm(whi.aa)
rm(whi.as)
rm(whi.ha)
whi$JWsource = "dbGaP_Ex"
mapWHI <- map[which(map$study == "WHI"),] #n=10024
whi <- merge(mapWHI,whi,by.x='individual_id',by.y='individual_id',all.x=TRUE) #n=10024
rm(mapWHI)

whi2 <- whi[is.na(whi$t2d),] 
#n=133 with WGS & no phenotype

# recode & check variable names & distributions
table(whi$t2d,useNA='always') # ! LOTS of missing data !!!! n=133
table(whi$ancestry,useNA='always')
table(whi$sex,useNA='always')
whi$origsex = whi$sex
whi$sex[whi$sex == 2] = 'F'
with(whi,table(origsex,sex,useNA='always'))
summary(whi$last_exam_age)
summary(whi$last_exam_bmi) #!low BMI
whi$study_ancestry <- paste(whi$study,whi$ancestry, sep = "_")

whi <- whi[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=10024
table(whi$study_ancestry,useNA='always')
#
# GenSalt Phase 2
#

gensalt$JWsource = "dbGaP_Ex"
mapGS <- map[which(map$study == "GenSalt"),] #n=1684
gensalt <- merge(mapGS,gensalt,by.x='individual_id',by.y='Individual_ID',all.x=TRUE) #n=1684
rm(mapGS)

# recode & check variable names & distributions
table(gensalt$Sequenced,useNA='always')
#gensalt = subset(gensalt, Sequenced == 1) #n=1872 drop 34 individuals
table(gensalt$T2D,useNA='always')
#gensalt$t2d = ifelse(is.na(gensalt$T2D),-9,gensalt$T2D)
gensalt$t2d[gensalt$T2D == 9] = -9
with(gensalt,table(T2D,t2d,useNA='always'))
gensalt$ancestry = 'AS'
table(gensalt$ancestry,useNA='always')
table(gensalt$SEX,useNA='always')
gensalt$sex[gensalt$SEX == 1] = 'M'
gensalt$sex[gensalt$SEX == 2] = 'F'
with(gensalt,table(SEX,sex,useNA='always'))
#describeBy(gensalt$Last_exam_age, gensalt$Sequenced, mat = TRUE)
summary(gensalt$Last_exam_age,useNA='always') #! younger than 18
#gensalt = subset(gensalt, Last_exam_age >= 18) 
gensalt$Last_exam_BMI <- as.numeric(as.character(gensalt$Last_exam_BMI)) #convert class character to numeric
summary(gensalt$Last_exam_BMI)
summary(gensalt$T2D_age) #drop 0 individuals
gensalt$last_exam_age = gensalt$Last_exam_age
gensalt$last_exam_bmi = gensalt$Last_exam_BMI
gensalt$last_exam_fg = gensalt$Last_exam_FG
gensalt$last_exam_hba1c = gensalt$Last_exam_HbA1c
gensalt$last_exam_t2d_treatment = gensalt$Last_exam_T2D_treatment
gensalt$t2d_age = gensalt$T2D_age
gensalt$t2d_bmi = gensalt$T2D_BMI
table(gensalt$Last_exam_visit,useNA = 'always')
gensalt$FamilyID = gensalt$Family_ID
gensalt$PaternalID = gensalt$Father_ID
gensalt$MaternalID = gensalt$Mother_ID
gensalt$study_ancestry <- paste(gensalt$study,gensalt$ancestry, sep = "_")

gensalt <- gensalt[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                      'study_topmedid','study_ancestry','JWsource','ancestry',
                      "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                      "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                      "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=1684

#
# GeneSTAR
#

genestar$JWsource = "dbGaP_Ex"
mapgenestar <- map[which(map$study == "GeneSTAR"),] #n=1634
genestar <- merge(mapgenestar,genestar,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=1634
rm(mapgenestar)

# recode & check variable names & distributions
table(genestar$sequenced,useNA='always')
table(genestar$T2D,useNA='always') #NA=0
genestar$t2d = ifelse(is.na(genestar$T2D),-9,genestar$T2D)
with(genestar,table(T2D,t2d,useNA='always'))
table(genestar$ancestry,useNA='always')
table(genestar$sex,useNA='always')
genestar$origsex = genestar$sex
genestar$sex[genestar$origsex == 1] = 'M'
genestar$sex[genestar$origsex == 2] = 'F'
with(genestar,table(origsex,sex,useNA='always'))
summary(genestar$last_exam_age,useNA='always') # 66 MISSING
#genestar = subset(genestar, last_exam_age >= 25) #n=, drop XX individuals
summary(genestar$T2D_age)
#genestar = subset(genestar, T2D_age >= 25) #n=, drop XX individuals
summary(genestar$last_exam_BMI) # ! low BMI & missing
genestar$last_exam_bmi = genestar$last_exam_BMI
genestar$last_exam_fg = genestar$last_exam_FG
genestar$last_exam_hba1c = genestar$last_exam_HbA1C
genestar$last_exam_t2d_treatment = genestar$last_exam_T2D_treatment
genestar$t2d_age = genestar$T2D_age
genestar$t2d_bmi = genestar$T2D_BMI
table(genestar$last_exam_visit,useNA = 'always')
genestar$FamilyID = genestar$FAMILY_ID
genestar$PaternalID = genestar$FATHER
genestar$MaternalID = genestar$MOTHER
genestar$study_ancestry <- paste(genestar$study,genestar$ancestry, sep = "_")
genestar <- genestar[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                        'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                        'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                        'study_topmedid','study_ancestry','JWsource','ancestry',
                        "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                        "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                        "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]

#n=1634
table(genestar$sex,useNA='always')
#
# MESA 
# download 1NOV2017

mesa$JWsource = "dbGaP_Ex"
mapmesa <-map[which(map$topmed_project == 'MESA' & map$study == "MESA"),]  #n=4124
mesa <- merge(mapmesa,mesa,by.x='individual_id',by.y='Individual_ID',all.x=TRUE) #n=4254
rm(mapmesa)

# recode & check variable names & distributions
table(mesa$sequenced,useNA='always')
table(mesa$T2D,useNA='always')
mesa$t2d = ifelse(is.na(mesa$T2D),NA,mesa$T2D)
with(mesa,table(T2D,t2d,useNA='always'))
table(mesa$ancestry,useNA='always')
table(mesa$Sex,useNA='always')
mesa$sex[mesa$Sex == 1] = 'M'
mesa$sex[mesa$Sex == 2] = 'F'
with(mesa,table(Sex,sex,useNA='always'))
summary(mesa$last_exam_age,useNA='always')
summary(mesa$last_exam_bmi) # ! low BMI
mesa$last_exam_fg = mesa$last_exam_FG
summary(mesa$last_exam_fg)
summary(mesa$last_exam_HBa1c)
mesa$last_exam_hba1c = mesa$last_exam_HBa1c
mesa$last_exam_t2d_treatment = mesa$last_exam_T2D_treatment
mesa$t2d_age <- as.numeric(as.character(mesa$T2D_Age))
summary(mesa$t2d_age)
mesa$t2d_bmi <- as.numeric(as.character(mesa$T2D_BMI))
summary(mesa$t2d_bmi)
table(mesa$last_exam_visit,useNA = 'always')
mesa$FamilyID = mesa$Family_ID
mesa$PaternalID = mesa$Father_ID
mesa$MaternalID = mesa$Mother_ID
mesa$study_ancestry <- paste(mesa$study,mesa$ancestry, sep = "_")
mesa <- mesa[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                'study_topmedid','study_ancestry','JWsource','ancestry',
                "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]

#n=4124
table(mesa$sex,useNA='always')

#
# MESA Family 
## download 1NOV2017

mesafam$JWsource = "dbGaP_Ex"
mesafam$ancestry = 'AA'
mapmesafam <- map[which(map$topmed_project == 'AA_CAC' & map$study == "MESA"),] #n=690
mesafam <- merge(mapmesafam,mesafam,by.x='individual_id',by.y='INDID',all.x=TRUE) #n=690
rm(mapmesafam)

# recode & check variable names & distributions
table(mesafam$sequenced,useNA='always')
#mesafam = subset(mesafam, sequenced == 1) #n=694 drop 8 individuals
table(mesafam$T2D,useNA='always')
mesafam$t2d = ifelse(is.na(mesafam$T2D),-9,mesafam$T2D)
with(mesafam,table(T2D,t2d,useNA='always'))
table(mesafam$ancestry,useNA='always')
table(mesafam$SEX,useNA='always')
mesafam$sex[mesafam$SEX == 1] = 'M'
mesafam$sex[mesafam$SEX == 2] = 'F'
with(mesafam,table(SEX,sex,useNA='always'))
summary(mesafam$last_exam_age,useNA='always')
summary(mesafam$last_exam_BMI) # ! low BMI
mesafam$last_exam_bmi = mesafam$last_exam_BMI
mesafam$last_exam_fg = mesafam$last_exam_FG
mesafam$last_exam_hba1c = NA
mesafam$last_exam_t2d_treatment = mesafam$last_exam_T2D_treatment
table(mesafam$last_exam_visit,useNA = 'always')
mesafam$FamilyID = mesafam$FAMID
mesafam$PaternalID = mesafam$FAID
mesafam$MaternalID = mesafam$MOID
mesafam$study_ancestry <- paste(mesafam$study,mesafam$ancestry, sep = "_")

mesafam <- mesafam[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                      'study_topmedid','study_ancestry','JWsource','ancestry',
                      "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                      "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                      "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=690

#
# CHS
# 

chs$JWsource = "dbGaP_Ex"
mapchs <- map[which(map$study == "CHS"),] #n=69
chs <- merge(mapchs,chs,by.x='individual_id',by.y='SampleID',all.x=TRUE) #n=69
rm(mapchs)

# recode & check variable names & distributions
table(chs$sequenced,useNA='always')
table(chs$T2D,useNA='always')
chs$t2d = ifelse(is.na(chs$T2D),-9,chs$T2D)
with(chs,table(T2D,t2d,useNA='always'))
table(chs$race01,useNA='always')
chs$ancestry[chs$race01 == 2] = 'AF'
chs$ancestry[chs$race01 == 1] = 'EU'
chs$ancestry[chs$race01 == 5] = 'OTHER'
with(chs,table(race01,ancestry,useNA='always'))
table(chs$sex,useNA='always')
chs$sex[chs$sex == 1] = 'M'
chs$sex[chs$sex == 2] = 'F'
table(chs$sex,useNA='always')
summary(chs$last_exam_age)
summary(chs$last_exam_BMI)
summary(chs$last_exam_FG)
summary(chs$T2D_age)
chs$last_exam_bmi = chs$last_exam_BMI
chs$last_exam_fg = chs$last_exam_FG
chs$last_exam_hba1c = NA
chs$t2d_age = chs$T2D_age
chs$t2d_bmi = chs$T2D_BMI
chs$last_exam_t2d_treatment = chs$last_exam_T2D_treatment
chs$last_exam_visit = NA
chs$study_ancestry <- paste(chs$study,chs$ancestry, sep = "_")

chs <- chs[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=69

## GENOA download 31OCT2017

genoa$JWsource = "dbGaP_Ex"
mapgenoa <- map[which(map$study == "GENOA"),] #n=1141
genoa <- merge(mapgenoa,genoa,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=1141
rm(mapgenoa)

# recode & check variable names & distributions
table(genoa$sequenced,useNA='always')
table(genoa$T2D,useNA='always')
genoa$t2d = ifelse(is.na(genoa$T2D),-9,genoa$T2D)
with(genoa,table(T2D,t2d,useNA='always'))
genoa$ancestry = 'AF'
table(genoa$ancestry,useNA='always')
table(genoa$sex,useNA='always')
genoa$sex[genoa$sex == 1] = 'M'
genoa$sex[genoa$sex == 2] = 'F'
table(genoa$sex,useNA='always')
summary(genoa$last_exam_age,useNA='always') # used -9 for missing & NA
summary(genoa$last_exam_BMI)
summary(genoa$T2D_AGE)
summary(genoa$T2D_BMI)
genoa$t2d_age = genoa$T2D_AGE
genoa$t2d_bmi = genoa$T2D_BMI
genoa$last_exam_bmi = genoa$last_exam_BMI
genoa$last_exam_fg = genoa$last_exam_FG
genoa$last_exam_hba1c = genoa$last_exam_HbA1c
genoa$last_exam_t2d_treatment = genoa$last_exam_T2D_treatment
genoa$FamilyID = genoa$FAMILY_ID
genoa$PaternalID = genoa$FATHER
genoa$MaternalID = genoa$MOTHER
genoa$study_ancestry <- paste(genoa$study,genoa$ancestry, sep = "_")

genoa <- genoa[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                  'study_topmedid','study_ancestry','JWsource','ancestry',
                  "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                  "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                  "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=1141


## DHS download

dhs$JWsource = "dbGaP_Ex_General"

names(dhsPed)
# "FAMILY_ID"  "SUBJECT_ID" "FATHER"     "MOTHER"     "SEX" 
dhsPed <- dhsPed[,c("FAMILY_ID",  "SUBJECT_ID", "FATHER", "MOTHER")]
dhs <- merge(dhs,dhsPed,by.x='SUBJECT_ID',by.y='SUBJECT_ID',all.x=TRUE) #n=405 & 33 variables
rm(dhsPed)
mapdhs <- map[which(map$study == "DHS"),] #n=337
dhs <- merge(mapdhs,dhs,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=337
rm(mapdhs)

# recode & check variable names & distributions
table(dhs$DIABETES,useNA='always')
dhs$t2d = ifelse(is.na(dhs$DIABETES),-9,dhs$DIABETES)
dhs$t2d[dhs$DIABETES == 1] = 2
dhs$t2d[dhs$DIABETES == 0] = 0
with(dhs,table(DIABETES,t2d,useNA='always'))
table(dhs$RACE,useNA='always')
dhs$ancestry[dhs$RACE == 'AA'] = 'AF'
table(dhs$ancestry,useNA='always')
table(dhs$SEX,useNA='always')
dhs$sex = dhs$SEX
table(dhs$sex,useNA='always')
dhs$AGE<- as.numeric(dhs$AGE)
summary(dhs$AGE)
summary(dhs$BMI)
summary(dhs$HBA1C)
summary(dhs$GLUCOSE)
dhs$last_exam_fg = dhs$GLUCOSE*0.0555  # !! convert !!
summary(dhs$last_exam_fg)
dhs$t2d_age = NA
dhs$t2d_bmi = NA
dhs$last_exam_age = dhs$AGE
dhs$last_exam_bmi = dhs$BMI
dhs$last_exam_hba1c = dhs$HBA1C
dhs$last_exam_t2d_treatment = NA
dhs$FamilyID = dhs$FAMILY_ID
dhs$PaternalID = dhs$FATHER
dhs$MaternalID = dhs$MOTHER
dhs$study_ancestry <- paste(dhs$study,dhs$ancestry, sep = "_")

dhs <- dhs[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=337


###### ADDED THIS RECODE OF STUDIES  13FEB2018
## BEGIN

# AWAITING PHENOTYPE DATA
#
# SAFS CVD
## 9 SEP 2017 DOWNLOADED new dataset
# 7DEC2017 checked for new files - none
safs = read.csv('SAFSCVD_HA_MAHANEY_20170807_T2D.ped.csv', header=T,sep=',',as.is=TRUE) #n=2457 (n=2 Sequenced=0)
safs$Individual_ID_pheno = safs$Individual_ID
safs$JWsource = "dbGaP_Ex"

safsid = read.csv("SAFSCVD_PERALTA_09262017_nwd_mappingtable.csv", header=T,sep=',',as.is=TRUE) #n=2713
safs2 <- merge(safsid,safs,by.x='NWD.ID',by.y='Individual_ID',all.x=TRUE) #n=2715
write.csv(safs2,row.names=F,quote=F,file='SAFS.test2.17NOV2017.T2D.csv')

mapSAFS <- map[which(map$study == "SAFS"),] #n=1509
safs3 <- merge(mapSAFS,safs,by.x='sample.id',by.y='Individual_ID',all.x=TRUE,all.y=TRUE) #n=1509
rm(mapSAFS)
write.csv(safs3,row.names=F,quote=F,file='SAFS.test.25JAN2018.T2D.csv')

# recode & check variable names & distributions
table(safs$Sequenced,useNA='always')
table(safs$T2D,useNA='always')
safs$t2d = ifelse(is.na(safs$T2D),-9,safs$T2D)
with(safs,table(T2D,t2d,useNA='always'))
safs$ancestry = 'HS'
table(safs$ancestry,useNA='always')
table(safs$Sex,useNA='always')
safs$sex[safs$Sex == 1] = 'M'
safs$sex[safs$Sex == 2] = 'F'
with(safs,table(Sex,sex,useNA='always'))
summary(safs$Last_Exam_Age,useNA='always')
safs = subset(safs, Last_Exam_Age >= 25) #n=1903, drop 552 individuals
summary(safs$Last_Exam_BMI) # ! low BMI
safs$last_exam_age = safs$Last_Exam_Age
safs$last_exam_bmi = safs$Last_Exam_BMI
safs$last_exam_fg = safs$Last_Exam._G
safs$last_exam_hba1c = NA
## variables are missing contacted MM 26JUL2017
safs$last_exam_t2d_treatment = safs$last_exam_T2D_treatment
safs$t2d_age = safs$T2D_age
safs$t2d_bmi = safs$T2D_BMI
table(safs$Last_exam_visit,useNA = 'always')
safs$FamilyID = safs$Family_ID
safs$PaternalID = safs$Father_ID
safs$MaternalID = safs$Mother_ID
#safs$study_ancestry <- paste(safs$study,safs$ancestry, sep = "_")

safs <- safs[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                'study_topmedid','study_ancestry','JWsource','ancestry',
                "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=????

### AWAITING GENOTYPE DATA
##
# AFib Australia
## new file downloaded from dbGaP 30OCT2017
Aus =read.table('Fatkin_Australia_dbGaP_SubjectPhenotypesDS_v1.txt',
                header=T,sep='\t',as.is=T,fill = TRUE) #n=120 & 17 variables
Aus$JWsource = "dbGaP"
mapAus <- map[map$study == "",] #n=
Aus2 <- merge(mapAus,Aus,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=789
rm(mapAus)

# recode & check variable names & distributions
table(Aus$diabetes,useNA='always')
Aus$t2d = ifelse(is.na(Aus$diabetes),-9,Aus$diabetes)
Aus$t2d[Aus$t2d == 'yes'] = 2
Aus$t2d[Aus$t2d == 'no'] = 0
with(Aus,table(diabetes,t2d,useNA='always')) #n=5 unknown T2D
table(Aus$race,useNA='always')
table(Aus$ethnicity,useNA='always')
Aus$ancestry[Aus$race == 'white' & Aus$ethnicity == 'no'] = "EU"
table(Aus$ancestry,useNA='always')
table(Aus$sex,useNA='always')
Aus$origsex = Aus$sex
Aus$sex[Aus$sex == 'male'] = 'M'
Aus$sex[Aus$sex == 'female'] = 'F'
with(Aus,table(origsex,sex,useNA='always'))
summary(Aus$age)
Aus = subset(Aus, age >= 25) #drop 1 individuals
Aus$last_exam_age = Aus$age
summary(Aus$height)
summary(Aus$weight)
Aus$last_exam_bmi = ((703*Aus$weight)/(Aus$height*Aus$height))
summary(Aus$last_exam_bmi) 
Aus$last_exam_fg = NA
Aus$last_exam_hba1c = NA
Aus$last_exam_t2d_treatment = NA
Aus$last_exam_visit = NA
Aus$t2d_age = NA
Aus$t2d_bmi = NA
Aus$FamilyID = Aus$individual_id
Aus$PaternalID = 0
Aus$MaternalID = 0
Aus$study_ancestry <- paste(Aus$study,Aus$ancestry, sep = "_")

Aus <- Aus[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
              'study_topmedid','study_ancestry','JWsource','ancestry',
              "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
              "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
              "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=65
write.csv(Aus,row.names=F,quote=F,file='Harmonized.AFib.Aus.XXXX2017.T2D.csv')


##
# AFib UMass
## new file downloaded from dbGaP 30OCT2017
UMass =read.table('UMASS_dbGaP_SubjectPhenotypesDS_v1.txt',
                  header=T,sep='\t',as.is=T,fill = TRUE) #n=65 & 17 variables
UMass$JWsource = "dbGaP"
mapUMass <- map[map$study == "",] #n=
UMass2 <- merge(mapUMass,UMass,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=789
rm(mapUMass)

# recode & check variable names & distributions
table(UMass$diabetes,useNA='always')
UMass$t2d = ifelse(is.na(UMass$diabetes),-9,UMass$diabetes)
UMass$t2d[UMass$t2d == 'yes'] = 2
UMass$t2d[UMass$t2d == 'no'] = 0
with(UMass,table(diabetes,t2d,useNA='always')) #n=5 unknown T2D
table(UMass$race,useNA='always')
table(UMass$ethnicity,useNA='always')
UMass$ancestry[UMass$race == 'white' & UMass$ethnicity == 'no'] = "EU"
table(UMass$ancestry,useNA='always')
table(UMass$sex,useNA='always')
UMass$origsex = UMass$sex
UMass$sex[UMass$sex == 'male'] = 'M'
UMass$sex[UMass$sex == 'female'] = 'F'
with(UMass,table(origsex,sex,useNA='always'))
summary(UMass$age)
UMass$last_exam_age = UMass$age
summary(UMass$height)
summary(UMass$weight)
UMass$last_exam_bmi = ((703*UMass$weight)/(UMass$height*UMass$height))
summary(UMass$last_exam_bmi) 
UMass$last_exam_fg = NA
UMass$last_exam_hba1c = NA
UMass$last_exam_t2d_treatment = NA
UMass$last_exam_visit = NA
UMass$t2d_age = NA
UMass$t2d_bmi = NA
UMass$FamilyID = UMass$individual_id
UMass$PaternalID = 0
UMass$MaternalID = 0
UMass$study_ancestry <- paste(UMass$study,UMass$ancestry, sep = "_")


UMass <- UMass[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                  'study_topmedid','study_ancestry','JWsource','ancestry',
                  "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                  "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                  "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=65
write.csv(UMass,row.names=F,quote=F,file='Harmonized.AFib.UMass.XXXX2017.T2D.csv')


## THRV download
thrv = read.table('THRV_Chinese_KS_09142017_T2D.ped', header=T,sep='\t',as.is=T) #n=2346
thrv$JWsource = "dbGaP_Ex"
mapthrv <- map[map$study == "THRV",] #n=1144
thrv <- merge(mapthrv,thrv,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=1144
rm(mapthrv)

# recode & check variable names & distributions
table(thrv$sequenced,useNA='always')
table(thrv$T2D,useNA='always')
thrv$t2d = ifelse(is.na(thrv$T2D),-9,thrv$T2D)
with(thrv,table(T2D,t2d,useNA='always'))
thrv$ancestry = 'AS'
table(thrv$ancestry,useNA='always')
table(thrv$sex,useNA='always')
thrv$sex[thrv$sex == 1] = 'M'
thrv$sex[thrv$sex == 2] = 'F'
table(thrv$sex,useNA='always')
summary(thrv$last_exam_age,useNA='always')
summary(thrv$last_exam_BMI)
summary(thrv$T2D_age)
summary(thrv$T2D_BMI)
summary(thrv$last_exam_HbA1c)
thrv$t2d_age = thrv$T2D_AGE
thrv$t2d_bmi = thrv$T2D_BMI
thrv$last_exam_bmi = thrv$last_exam_BMI
thrv$last_exam_fg = thrv$last_exam_FG
thrv$last_exam_hba1c = thrv$last_exam_HbA1c
thrv$last_exam_t2d_treatment = thrv$last_exam_T2D_treatment
thrv$FamilyID = thrv$family_ID
thrv$PaternalID = thrv$father_id
thrv$MaternalID = thrv$mother_id
thrv$study_ancestry <- paste(thrv$study,thrv$ancestry, sep = "_")

thrv <- thrv[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
                'study_topmedid','study_ancestry','JWsource','ancestry',
                "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
                "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
                "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
#n=

##END OF RECODE



########### CREATE POOLED DATASET ###########################
# 
# studies to add in future safs,thrv, goldn, hypergen
pooled <- rbind(afccaf,afp,afvub,amish,aric,cfs,chs,copd,dhs,fhs,genestar,genoa,
                gensalt,jhs,mesa,mesafam,raw.HVH,raw.MGH,raw.VUdaw,sas,whi)

table(pooled$t2d,useNA='always') #-9 n=2834 NA n=1978
#pooled$t2d[pooled$t2d == -9] = 'NA'
table(pooled$t2d)
table(pooled$sex,useNA='always') # 595 missing 1 XXY
pooled2 <- pooled[is.na(pooled$t2d),] 
pooled2 = subset(pooled, sex == 'XXY')
table(pooled$CENTER,useNA='always')
#  baylor    broad illumina macrogen     nygc       uw     <NA> 
#  5308    29938      264     1283      851     7188        0 

##  !!  NEED TO CREATE THESE VARIABLES !!  ##
## t2d_ctrl t2d_superctrl
table(pooled$t2d,useNA='always')
pooled$t2d_ctrl[pooled$t2d == 2] = 1
pooled$t2d_ctrl[pooled$t2d == 1] = 0
pooled$t2d_ctrl[pooled$t2d == 0] = 0
pooled$t2d_ctrl[pooled$t2d == -9] = -9
with(pooled,table(t2d,t2d_ctrl,useNA='always'))

table(pooled$t2d,useNA='always')
pooled$t2d_nopre.ctrl[pooled$t2d == 2] = 1
pooled$t2d_nopre.ctrl[pooled$t2d == 1] = -9
pooled$t2d_nopre.ctrl[pooled$t2d == 0] = 0
pooled$t2d_nopre.ctrl[pooled$t2d == -9] = -9
with(pooled,table(t2d,t2d_nopre.ctrl,useNA='always'))

write.csv(pooled,row.names=F,quote=F,file='Pooled_MIXED_WesselJ_18JAN2017_T2D.csv')


#collapse AFib studies
table(pooled$study_ancestry,useNA='always')
pooled$study_ancestry_AFib[pooled$study_ancestry == "ARIC_EU"] = "Afib_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "CCAF_EU"] = "Afib_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "HVH_EU"] = "Afib_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "MGH_AF_EU"] = "Afib_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "Partners_EU"] = "Afib_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "VAFAR_EU"] = "Afib_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "VU_AF_EU"] = "Afib_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "Partners_AF"] = "Afib_AF"
pooled$study_ancestry_AFib[pooled$study_ancestry == "VU_AF_AF"] = "Afib_AF"
pooled$study_ancestry_AFib[pooled$study_ancestry == "CCAF_HS"] = "Afib_HS"
pooled$study_ancestry_AFib[pooled$study_ancestry == "MGH_AF_HS"] = "Afib_HS"
pooled$study_ancestry_AFib[pooled$study_ancestry == "VU_AF_HS"] = "Afib_HS"
pooled$study_ancestry_AFib[pooled$study_ancestry == "HVH_OTHER"] = "Afib_Other"
pooled$study_ancestry_AFib[pooled$study_ancestry == "Partners_OTHER"] = "Afib_Other"
pooled$study_ancestry_AFib[pooled$study_ancestry == "MGH_AF_AMR"] = "Afib_AMR"
pooled$study_ancestry_AFib[pooled$study_ancestry == "VU_AF_AMR"] = "Afib_AMR"
pooled$study_ancestry_AFib[pooled$study_ancestry == "MGH_AF_MIXED"] = "Afib_Mixed"
pooled$study_ancestry_AFib[pooled$study_ancestry == "Partners_AS"] = "Afib_AS"
pooled$study_ancestry_AFib[pooled$study_ancestry == "VU_AF_AS"] = "Afib_AS"
pooled$study_ancestry_AFib[pooled$study_ancestry == "Partners_NA"] = "Afib_NA"
pooled$study_ancestry_AFib[pooled$study_ancestry == "Amish_EU"] = "Amish_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "CFS_AF"] = "CFS_AF"
pooled$study_ancestry_AFib[pooled$study_ancestry == "CFS_EU"] = "CFS_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "COPDGene_AF"] = "COPDGene_AF"
pooled$study_ancestry_AFib[pooled$study_ancestry == "COPDGene_EU"] = "COPDGene_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "FHS_EU"] = "FHS_EU"
pooled$study_ancestry_AFib[pooled$study_ancestry == "JHS_AF"] = "JHS_AF"
pooled$study_ancestry_AFib[pooled$study_ancestry == "SAS_AS"] = "SAS_AS"
with(pooled,table(study_ancestry,study_ancestry_AFib,useNA='always'))

write.csv(pooled,row.names=F,quote=F,file='Pooled_MIXED_WesselJ_27AUG2017_T2D.csv')
save(pooled, file = "Pooled_MIXED_WesselJ_27AUG2017_T2D.RData")
