

####read duplicates file
d<-read.table("/restricted/projectnb/sequencing/achilleas/dbGaP-11723/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.5b/relatedness/freeze5b_duplicates.txt",header=T)
####read linker file
linker<-read.table("/restricted/projectnb/sequencing/achilleas/dbGaP-11723/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.5b/sample_annotation/sample_sets_2017-10-16/freeze5b_sample_annot_2017-10-16.txt",header=T)
t<-table(linker$study,linker$CENTER)
t<-prop.table(t,1)

###add center
linker<-linker[,c("sample.id","CENTER")]
names(linker)<-c("sample.id","center1")
d<-merge(d,linker,by.x="ID1",by.y="sample.id",all.x=T)
names(linker)<-c("sample.id","center2")
d<-merge(d,linker,by.x="ID2",by.y="sample.id",all.x=T)

####add call rate
q<-read.table("/restricted/projectnb/glycemic/achilleas/pass_minDP10_qc.tsv",header=T)

q<-q[,c("s","sa.qc.callRate")]

names(q)<-c("s","cr1")
d<-merge(d,q,by.x="ID1",by.y="s",all.x=T)

names(q)<-c("s","cr2")
d<-merge(d,q,by.x="ID2",by.y="s",all.x=T)

###add center percentage

for (j in 1:nrow(d))
{
  d$p1[j]<-ifelse(d$study1[j]%in%row.names(t),t[as.character(d$study1[j]),as.character(d$center1[j])],NA)
  d$p2[j]<-ifelse(d$study2[j]%in%row.names(t),t[as.character(d$study2[j]),as.character(d$center2[j])],NA)
  
  }

 ##remove controls
d<-d[which(!(d$study1%in%c("CONTROL"))&!(d$study2%in%c("CONTROL"))),]

for (j in 1:nrow(d))
{
  if (as.character(d$study1[j])==as.character(d$study2[j])) ###within same study
  {
    if (is.na(d$MZtwinID[j])==F){d$keep1[j]<-1;d$keep2[j]<-1} ###keep MZtwin
     else if (as.character(d$center1[j])!=as.character(d$center2[j])&is.na(d$p1[j])==F&is.na(d$p2[j])==F){
      d$keep1[j]<-ifelse(d$p1[j]>=d$p2[j],1,0)
      d$keep2[j]<-1-d$keep1[j]
     }
     else if (as.character(d$center1[j])==as.character(d$center2[j])){
       d$keep1[j]<-ifelse(d$cr1[j]>=d$cr2[j],1,0)
       d$keep2[j]<-1-d$keep1[j]
     }###Choose duplicate with highest sequencing call rate
  }

  #else {d$keep1[j]<-NA;d$keep2[j]<-NA} 
  #ARIC>     DHS, GOLDN, GENOA, HyperGEN, Mayo_VTE
   else if ((as.character(d$study1[j])=='ARIC'&as.character(d$study2[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE"))|
            (as.character(d$study2[j])=='ARIC'&as.character(d$study1[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='ARIC',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }#FHS>       COPDGene, MGH_AF, WHI
    else if ((as.character(d$study1[j])=='FHS'&as.character(d$study2[j])%in%c("COPDGene","MGH_AF","WHI"))|
             (as.character(d$study2[j])=='FHS'&as.character(d$study1[j])%in%c("COPDGene","MGH_AF","WHI")))
    {
      d$keep1[j]<-ifelse(as.character(d$study1[j])=='FHS',1,0)
      d$keep2[j]<-1-d$keep1[j]
    }#GeneSTAR>          COPDGene
   else if ((as.character(d$study1[j])=='GeneSTAR'&as.character(d$study2[j])%in%c("COPDGene"))|
            (as.character(d$study2[j])=='GeneSTAR'&as.character(d$study1[j])%in%c("COPDGene")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='GeneSTAR',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }#HyperGEN>          COPDGene
  else if ((as.character(d$study1[j])=='HyperGEN'&as.character(d$study2[j])%in%c("COPDGene"))|
           (as.character(d$study2[j])=='HyperGEN'&as.character(d$study1[j])%in%c("COPDGene")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='HyperGEN',1,0)
    d$keep2[j]<-1-d$keep1[j]
  }
  ##JHS>       GENOA
   else if ((as.character(d$study1[j])=='JHS'&as.character(d$study2[j])%in%c("GENOA"))|
            (as.character(d$study2[j])=='JHS'&as.character(d$study1[j])%in%c("GENOA")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='JHS',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }
  ##MESA>   COPDGene, DHS, GENOA, HyperGEN
  else if ((as.character(d$study1[j])=='MESA'&as.character(d$study2[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN"))|
           (as.character(d$study2[j])=='MESA'&as.character(d$study1[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='MESA',1,0)
    d$keep2[j]<-1-d$keep1[j]
  }
  ###VU_AF> WGHS
  else if ((as.character(d$study1[j])=='VU_AF'&as.character(d$study2[j])%in%c("WGHS"))|
           (as.character(d$study2[j])=='VU_AF'&as.character(d$study1[j])%in%c("WGHS")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='VU_AF',1,0)
    d$keep2[j]<-1-d$keep1[j]
  }

   ###WHI>     COPDGene, HyperGEN, Mayo_VTE, MESA ######Added in MESA######
  else if ((as.character(d$study1[j])=='WHI'&as.character(d$study2[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA"))|
           (as.character(d$study2[j])=='WHI'&as.character(d$study1[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='WHI',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
  ###ARIC>WHI
  else if ((as.character(d$study1[j])=='ARIC'&as.character(d$study2[j])%in%c("WHI"))|
           (as.character(d$study2[j])=='ARIC'&as.character(d$study1[j])%in%c("WHI")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='ARIC',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
  
  ###ARIC > MESA
   else if ((as.character(d$study1[j])=='ARIC'&as.character(d$study2[j])%in%c("MESA"))|
           (as.character(d$study2[j])=='ARIC'&as.character(d$study1[j])%in%c("MESA")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='ARIC',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
  ####GeneSTAR > MESA
   else if ((as.character(d$study1[j])=='GeneSTAR'&as.character(d$study2[j])%in%c("MESA"))|
           (as.character(d$study2[j])=='GeneSTAR'&as.character(d$study1[j])%in%c("MESA")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='GeneSTAR',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
   ####CHS > MESA, WHI ###Added in CHS ordering####
  else if ((as.character(d$study1[j])=='CHS'&as.character(d$study2[j])%in%c("MESA", "WHI"))|
          (as.character(d$study2[j])=='CHS'&as.character(d$study1[j])%in%c("MESA", "WHI")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='CHS',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
  else { d$keep1[j]<-NA; d$keep2[j]<-NA}
}


####add phenotype

p<-read.table("/restricted/projectnb/glycemic/peitao/phenotype_harmonization/pooled_analysis/Pooled_Glycemic_Traits_freeze5_duplicate_ID.ped",header=T)


####FG
##This adds keep column for fasting glucose and creates an output listing removed duplicates



p1<-p[,c("TOPMEDID","FastingGlucose")]
names(p1)<-c("TOPMEDID","FastingGlucose1")
d<-merge(d,p1,by.x = 'ID1',by.y = 'TOPMEDID',all.x=T)

p1<-p[,c("TOPMEDID","FastingGlucose")]
names(p1)<-c("TOPMEDID","FastingGlucose2")
d<-merge(d,p1,by.x = 'ID2',by.y = 'TOPMEDID',all.x=T)

d$keep1_fg<-d$keep1
d$keep2_fg<-d$keep2

missing_fg <-c()
keep_fg <- c()
cohort_ranking_fg <-c()
callrate_fg <-c()

for (j in 1:nrow(d)){
	
		
	if(is.na(d$FastingGlucose1[j])==F&is.na(d$FastingGlucose2[j])==F){
			d$keep1_fg[j]<-ifelse(is.na(d$FastingGlucose1[j])==F,d$keep1[j],ifelse(is.na(d$FastingGlucose1[j])==T,0,NA))
			d$keep2_fg[j]<-ifelse(is.na(d$FastingGlucose2[j])==F,d$keep2[j],ifelse(is.na(d$FastingGlucose2[j])==T,0,NA))
		##for additional output fgles
			if(d$keep1[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
				keep_fg <- c(keep_fg,d$ID1[j])
				cohort_ranking_fg <-c(cohort_ranking_fg, as.character(d$ID2[j]))
			}
			else if(d$keep2[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
				keep_fg <- c(keep_fg,as.character(d$ID2[j]))
				cohort_ranking_fg <-c(cohort_ranking_fg, as.character(d$ID1[j]))
			}
			else if(d$keep1[j]==1&d$keep2[j]==1){
				keep_fg <- c(keep_fg,as.character(d$ID1[j]), as.character(d$ID2[j]))
			}
			else if(d$keep1[j]==1&d$keep2[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
				keep_fg <- c(keep_fg, as.character(d$ID1[j]))
				callrate_fg <-c(callrate_fg, as.character(d$ID2[j]))
			}
			else if(d$keep2[j]==1&d$keep1[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
				keep_fg <- c(keep_fg,as.character(d$ID2[j]))
				callrate_fg <-c(callrate_fg, as.character(d$ID1[j]))
			}
	}
	else if(is.na(d$FastingGlucose1[j])==T&is.na(d$FastingGlucose2[j])==F)
	{
		d$keep1_fg[j]<-0
		d$keep2_fg[j]<-1
		missing_fg<-c(missing_fg,as.character(d$ID1[j]))
		keep_fg<-c(keep_fg, as.character(d$ID2[j]))
	} 
	else if(is.na(d$FastingGlucose1[j])==F&is.na(d$FastingGlucose2[j])==T)
	{
		d$keep1_fg[j]<-1
		d$keep2_fg[j]<-0
		missing_fg<-c(missing_fg,as.character(d$ID2[j]))
		keep_fg<-c(keep_fg, as.character(d$ID1[j]))
	}
	else if (is.na(d$FastingGlucose1[j])==T&is.na(d$FastingGlucose2[j])==T)
	{
		d$keep1_fg[j]<-0
		d$keep2_fg[j]<-0
		missing_fg<-c(missing_fg, as.character(d$ID1[j]),as.character(d$ID2[j]))
	}
	else{ 
		d$keep1_fg[j]<-NA
		d$keep2_fg[j]<-NA
	}
}

table(d$keep1_fg,d$keep2_fg,useNA='always')


names <- c("ID","reason")

missing_fg.df <- data.frame(missing_fg)
missing_fg.df$reason<-"missing trait data"
colnames(missing_fg.df) <-names

cohort_ranking_fg.df <-data.frame(cohort_ranking_fg)
cohort_ranking_fg.df$reason<-"used duplicate from older study"
colnames(cohort_ranking_fg.df) <-names

callrate_fg.df <-data.frame(callrate_fg)
callrate_fg.df$reason <-"used duplicate with higher call rate" 
colnames(callrate_fg.df) <-names

not_keep_fg <- rbind(missing_fg.df,cohort_ranking_fg.df, callrate_fg.df)
colnames(not_keep_fg) <-names



write.table(not_keep_fg,"/data4/dloesch/Duplicates/Test/not_keep_fg.txt",row.names=F,col.names=T,quote=F,sep='\t')






####a1c


p2<-p[,c("TOPMEDID","HbA1c")]
names(p2)<-c("TOPMEDID","HbA1c1")
d<-merge(d,p2,by.x = 'ID1',by.y = 'TOPMEDID',all.x=T)

p2<-p[,c("TOPMEDID","HbA1c")]
names(p2)<-c("TOPMEDID","HbA1c2")
d<-merge(d,p2,by.x = 'ID2',by.y = 'TOPMEDID',all.x=T)


d$keep1_a1c<-d$keep1
d$keep2_a1c<-d$keep2

#for building output file for listing removed duplicates
missing_a1c <-c()
keep_a1c <- c()
cohort_ranking_a1c <-c()
callrate_a1c <-c()

for (j in 1:nrow(d)){
	
		
	if(is.na(d$HbA1c1[j])==F&is.na(d$HbA1c2[j])==F){
			d$keep1_a1c[j]<-ifelse(is.na(d$HbA1c1[j])==F,d$keep1[j],ifelse(is.na(d$HbA1c1[j])==T,0,NA))
			d$keep2_a1c[j]<-ifelse(is.na(d$HbA1c2[j])==F,d$keep2[j],ifelse(is.na(d$HbA1c2[j])==T,0,NA))
		##for additional output a1cles
			if(d$keep1[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
				keep_a1c <- c(keep_a1c,d$ID1[j])
				cohort_ranking_a1c <-c(cohort_ranking_a1c, as.character(d$ID2[j]))
			}
			else if(d$keep2[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
				keep_a1c <- c(keep_a1c,as.character(d$ID2[j]))
				cohort_ranking_a1c <-c(cohort_ranking_a1c, as.character(d$ID1[j]))
			}
			else if(d$keep1[j]==1&d$keep2[j]==1){
				keep_a1c <- c(keep_a1c,as.character(d$ID1[j]), as.character(d$ID2[j]))
			}
			else if(d$keep1[j]==1&d$keep2[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
				keep_a1c <- c(keep_a1c, as.character(d$ID1[j]))
				callrate_a1c <-c(callrate_a1c, as.character(d$ID2[j]))
			}
			else if(d$keep2[j]==1&d$keep1[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
				keep_a1c <- c(keep_a1c,as.character(d$ID2[j]))
				callrate_a1c <-c(callrate_a1c, as.character(d$ID1[j]))
			}
	}
	else if(is.na(d$HbA1c1[j])==T&is.na(d$HbA1c2[j])==F)
	{
		d$keep1_a1c[j]<-0
		d$keep2_a1c[j]<-1
		missing_a1c<-c(missing_a1c,as.character(d$ID1[j]))
		keep_a1c<-c(keep_a1c, as.character(d$ID2[j]))
	} 
	else if(is.na(d$HbA1c1[j])==F&is.na(d$HbA1c2[j])==T)
	{
		d$keep1_a1c[j]<-1
		d$keep2_a1c[j]<-0
		missing_a1c<-c(missing_a1c,as.character(d$ID2[j]))
		keep_a1c<-c(keep_a1c, as.character(d$ID1[j]))
	}
	else if (is.na(d$HbA1c1[j])==T&is.na(d$HbA1c2[j])==T)
	{
		d$keep1_a1c[j]<-0
		d$keep2_a1c[j]<-0
		missing_a1c<-c(missing_a1c, as.character(d$ID1[j]),as.character(d$ID2[j]))
	}
	else{ 
		d$keep1_a1c[j]<-NA
		d$keep2_a1c[j]<-NA
	}
}

table(d$keep1_a1c,d$keep2_a1c,useNA='always')

names <- c("ID","reason")

missing_a1c.df <- data.frame(missing_a1c)
missing_a1c.df$reason<-"missing trait data"
colnames(missing_a1c.df) <-names

cohort_ranking_a1c.df <-data.frame(cohort_ranking_a1c)
cohort_ranking_a1c.df$reason<-"used duplicate from older study"
colnames(cohort_ranking_a1c.df) <-names

callrate_a1c.df <-data.frame(callrate_a1c)
callrate_a1c.df$reason <-"used duplicate with higher call rate" 
colnames(callrate_a1c.df) <-names

not_keep_a1c <- rbind(missing_a1c.df,cohort_ranking_a1c.df, callrate_a1c.df)
colnames(not_keep_a1c) <-names

#for building output file for listing removed duplicates
write.table(not_keep_a1c,"/data4/dloesch/Test/Duplicates/not_keep_a1c.txt",row.names=F,col.names=T,quote=F,sep='\t')





####Insulin
p3<-p[,c("TOPMEDID", "FastingInsulin")]
names(p3) <- c("TOPMEDID", "FastingInsulin1")
d<-merge(d,p3,by.x = 'ID1',by.y = 'TOPMEDID',all.x=T)

p3<-p[,c("TOPMEDID", "FastingInsulin")]
names(p3) <- c("TOPMEDID", "FastingInsulin2")
d<-merge(d,p3,by.x = 'ID2',by.y = 'TOPMEDID',all.x=T)

d$keep1_fi<-d$keep1
d$keep2_fi<-d$keep2

#for building output file for listing removed duplicates
missing_fi <-c()
keep_fi <- c()
cohort_ranking_fi <-c()
callrate_fi <-c()

for (j in 1:nrow(d)){
	
		
	if(is.na(d$FastingInsulin1[j])==F&is.na(d$FastingInsulin2[j])==F){
			d$keep1_fi[j]<-ifelse(is.na(d$FastingInsulin1[j])==F,d$keep1[j],ifelse(is.na(d$FastingInsulin1[j])==T,0,NA))
			d$keep2_fi[j]<-ifelse(is.na(d$FastingInsulin2[j])==F,d$keep2[j],ifelse(is.na(d$FastingInsulin2[j])==T,0,NA))
	## for building output file listing removed duplicates
			if(d$keep1[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
				keep_fi <- c(keep_fi,d$ID1[j])
				cohort_ranking_fi <-c(cohort_ranking_fi, as.character(d$ID2[j]))
			}
			else if(d$keep2[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
				keep_fi <- c(keep_fi,as.character(d$ID2[j]))
				cohort_ranking_fi <-c(cohort_ranking_fi, as.character(d$ID1[j]))
			}
			else if(d$keep1[j]==1&d$keep2[j]==1){
				keep_fi <- c(keep_fi,as.character(d$ID1[j]), as.character(d$ID2[j]))
			}
			else if(d$keep1[j]==1&d$keep2[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
				keep_fi <- c(keep_fi, as.character(d$ID1[j]))
				callrate_fi <-c(callrate_fi, as.character(d$ID2[j]))
			}
			else if(d$keep2[j]==1&d$keep1[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
				keep_fi <- c(keep_fi,as.character(d$ID2[j]))
				callrate_fi <-c(callrate_fi, as.character(d$ID1[j]))
			}
	}
	else if(is.na(d$FastingInsulin1[j])==T&is.na(d$FastingInsulin2[j])==F)
	{
		d$keep1_fi[j]<-0
		d$keep2_fi[j]<-1
		missing_fi<-c(missing_fi,as.character(d$ID1[j]))
		keep_fi<-c(keep_fi, as.character(d$ID2[j]))
	} 
	else if(is.na(d$FastingInsulin1[j])==F&is.na(d$FastingInsulin2[j])==T)
	{
		d$keep1_fi[j]<-1
		d$keep2_fi[j]<-0
		missing_fi<-c(missing_fi,as.character(d$ID2[j]))
		keep_fi<-c(keep_fi, as.character(d$ID1[j]))
	}
	else if (is.na(d$FastingInsulin1[j])==T&is.na(d$FastingInsulin2[j])==T)
	{
		d$keep1_fi[j]<-0
		d$keep2_fi[j]<-0
		missing_fi<-c(missing_fi, as.character(d$ID1[j]),as.character(d$ID2[j]))
	}
	else{ 
		d$keep1_fi[j]<-NA
		d$keep2_fi[j]<-NA
	}
}

table(d$keep1_fi,d$keep2_fi,useNA='always')

names <- c("ID","reason")

missing_fi.df <- data.frame(missing_fi)
missing_fi.df$reason<-"missing trait data"
colnames(missing_fi.df) <-names

cohort_ranking_fi.df <-data.frame(cohort_ranking_fi)
cohort_ranking_fi.df$reason<-"used duplicate from older study"
colnames(cohort_ranking_fi.df) <-names

callrate_fi.df <-data.frame(callrate_fi)
callrate_fi.df$reason <-"used duplicate with higher call rate" 
colnames(callrate_fi.df) <-names

not_keep_fi <- rbind(missing_fi.df,cohort_ranking_fi.df,callrate_fi.df)
colnames(not_keep_fi) <-names

#this output file lists ID and reason for removed duplicates
write.table(not_keep_fi,"/data4/dloesch/Duplicates/Test/not_keep_fi.txt",row.names=F,col.names=T,quote=F,sep='\t')



#### for subsetting phenotype data for fasting glucose (and removing duplicates)


fg <- subset(p, (!is.na(p$FastingGlucose)))


not_keep1<-subset(d, keep1_fg==0)
not_keep1<-not_keep1[c("ID1","keep1_fg")]


not_keep2<-subset(d, keep2_fg==0)
not_keep2<-not_keep2[c("ID2","keep2_fg")]

names<-c("ID","KEEP")
not_keep<-data.frame()
colnames(not_keep1)<-names
colnames(not_keep2)<-names
not_keep <-rbind(not_keep1, not_keep2)


fg<-merge(fg,not_keep,by.x="TOPMEDID",by.y="ID",all.x=T)


fg$KEEP[is.na(fg$KEEP)] <- 1

fg <-subset(fg, KEEP==1)


write.table(fg,"/data4/dloesch/Duplicates/Test/removed_duplicates_fg.ped",row.names=F,col.names=T,quote=F,sep='\t')

### for subsetting phenotype data for HbA1c (and removing duplicates)

HbA1c <- subset(p, (!is.na(p$HbA1c)))


not_keep1<-subset(d, keep1_a1c==0)
not_keep1<-not_keep1[c("ID1","keep1_a1c")]


not_keep2<-subset(d, keep2_a1c==0)
not_keep2<-not_keep2[c("ID2","keep2_a1c")]

names<-c("ID","KEEP")
not_keep<-data.frame()
colnames(not_keep1)<-names
colnames(not_keep2)<-names
not_keep <-rbind(not_keep1, not_keep2)


HbA1c<-merge(HbA1c,not_keep,by.x="TOPMEDID",by.y="ID",all.x=T)


HbA1c$KEEP[is.na(HbA1c$KEEP)] <- 1

HbA1c <-subset(HbA1c, KEEP==1)


write.table(HbA1c,"/data4/dloesch/Duplicates/Test/removed_duplicates_HbA1c.ped",row.names=F,col.names=T,quote=F,sep='\t')



####to subset for fasting insulin and to remove duplicates from phenotype file 

fi <- subset(p, (!is.na(p$FastingInsulin)))


not_keep1<-subset(d, keep1_fi==0)
not_keep1<-not_keep1[c("ID1","keep1_fi")]


not_keep2<-subset(d, keep2_fi==0)
not_keep2<-not_keep2[c("ID2","keep2_fi")]

names<-c("ID","KEEP")
not_keep<-data.frame()
colnames(not_keep1)<-names
colnames(not_keep2)<-names
not_keep <-rbind(not_keep1, not_keep2)


fi<-merge(fi,not_keep,by.x="TOPMEDID",by.y="ID",all.x=T)


fi$KEEP[is.na(fi$KEEP)] <- 1

fi <-subset(fi, KEEP==1)


write.table(fi,"/data4/dloesch/Duplicates/Test/removed_duplicates_fi.ped",row.names=F,col.names=T,quote=F,sep='\t')




### for reordering so ID1 and ID2 are no longer flipped (ID2 was listed before ID1, but everything else was in the correct order). 

d <- d[c("ID1", "ID2", "study1", "study2", "MZtwinID", "center1", "center2", "cr1", "cr2", "p1", "p2", "keep1", "keep2","FastingGlucose1", "FastingGlucose2","keep1_fg",
"keep2_fg","HbA1c1","HbA1c2","keep1_a1c", "keep2_a1c","FastingInsulin1", "FastingInsulin2", "keep1_fi", "keep2_fi")]


###for the duplicates.txt file with keep columns
write.table(d,"/restricted/projectnb/glycemic/peitao/phenotype_harmonization/pooled_analysis/duplicates.txt",row.names=F,col.names=F,quote=F,sep='\t')