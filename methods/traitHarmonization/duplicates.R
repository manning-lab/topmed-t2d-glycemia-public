

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
  #ARIC>     DHS, GOLDN, GENOA, HyperGEN, Mayo_VTE
   else if ((as.character(d$study1[j])=='ARIC'&as.character(d$study2[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE"))|
            (as.character(d$study2[j])=='ARIC'&as.character(d$study1[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='ARIC',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }#FHS>       COPDGene, MGH_AF, WHI
    else if ((as.character(d$study1[j])=='FHS'&as.character(d$study2[j])%in%c("COPDGene","MGH_AF","WHI"))|
             (as.character(d$study2[j])=='FHS'&as.character(d$study1[j])%in%c("COPDGene","MGH_AF","WHI")))
    {
      d$keep1[j]<-ifelse(as.character(d$study1[j])=='FHS',1,0)
      d$keep2[j]<-1-d$keep1[j]
    }#GeneSTAR>          COPDGene
   else if ((as.character(d$study1[j])=='GeneSTAR'&as.character(d$study2[j])%in%c("COPDGene"))|
            (as.character(d$study2[j])=='GeneSTAR'&as.character(d$study1[j])%in%c("COPDGene")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='GeneSTAR',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }#HyperGEN>          COPDGene
  else if ((as.character(d$study1[j])=='HyperGEN'&as.character(d$study2[j])%in%c("COPDGene"))|
           (as.character(d$study2[j])=='HyperGEN'&as.character(d$study1[j])%in%c("COPDGene")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='HyperGEN',1,0)
    d$keep2[j]<-1-d$keep1[j]
  }
  ##JHS>       GENOA
   else if ((as.character(d$study1[j])=='JHS'&as.character(d$study2[j])%in%c("GENOA"))|
            (as.character(d$study2[j])=='JHS'&as.character(d$study1[j])%in%c("GENOA")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='JHS',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }
  ##MESA>   COPDGene, DHS, GENOA, HyperGEN
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

d<-d[which(!(d$study1%in%c("CONTROL"))&!(d$study2%in%c("CONTROL"))),]
####add phenotype

p<-read.table("/restricted/projectnb/glycemic/peitao/phenotype_harmonization/pooled_analysis/Pooled_Glycemic_Traits_freeze5_duplicate_ID.ped",header=T)


####FG
p1<-p[,c("TOPMEDID","FastingGlucose")]
names(p1)<-c("TOPMEDID","FastingGlucose1")
d<-merge(d,p1,by.x = 'ID1',by.y = 'TOPMEDID',all.x=T)

p1<-p[,c("TOPMEDID","FastingGlucose")]
names(p1)<-c("TOPMEDID","FastingGlucose2")
d<-merge(d,p1,by.x = 'ID2',by.y = 'TOPMEDID',all.x=T)

d$keep1_fg<-d$keep1
d$keep2_fg<-d$keep2

for (j in 1:nrow(d)){
  if (is.na(d$keep1[j])==T&is.na(d$keep2[j])==T) {
d$keep1_fg[j]<-ifelse(is.na(d$FastingGlucose1[j])==F&is.na(d$FastingGlucose2[j])==T,1,ifelse(is.na(d$FastingGlucose1[j])==T&is.na(d$FastingGlucose2[j])==F,0,NA))
d$keep2_fg[j]<-1-d$keep1_fg[j]
}
}

for (j in 1:nrow(d)){
  if (is.na(d$keep1[j])==T&is.na(d$keep2[j])==T) {
    if(is.na(d$FastingGlucose1[j])==T&is.na(d$FastingGlucose2[j])==T){
      d$keep1_fg[j]<-0
      d$keep2_fg[j]<-0
    }
  }
}
#d_fg<-d[which(is.na(d$keep1)==T&is.na(d$keep2)==T),]






####a1c


p2<-p[,c("TOPMEDID","HbA1c")]
names(p2)<-c("TOPMEDID","HbA1c1")
d<-merge(d,p2,by.x = 'ID1',by.y = 'TOPMEDID',all.x=T)

p2<-p[,c("TOPMEDID","HbA1c")]
names(p2)<-c("TOPMEDID","HbA1c2")
d<-merge(d,p2,by.x = 'ID2',by.y = 'TOPMEDID',all.x=T)


d$keep1_a1c<-d$keep1
d$keep2_a1c<-d$keep2

for (j in 1:nrow(d)){
  if (is.na(d$keep1[j])==T&is.na(d$keep2[j])==T) {
    d$keep1_a1c[j]<-ifelse(is.na(d$HbA1c1[j])==F&is.na(d$HbA1c2[j])==T,1,ifelse(is.na(d$HbA1c1[j])==T&is.na(d$HbA1c2[j])==F,0,NA))
    d$keep2_a1c[j]<-1-d$keep1_a1c[j]
  }
}

for (j in 1:nrow(d)){
  if (is.na(d$keep1[j])==T&is.na(d$keep2[j])==T) {
    if(is.na(d$HbA1c1[j])==T&is.na(d$HbA1c2[j])==T){
      d$keep1_a1c[j]<-0
      d$keep2_a1c[j]<-0
    }
  }
}


table(d$keep1_fg,d$keep2_fg,useNA='always')
table(d$keep1_a1c,d$keep2_a1c,useNA='always')



####Insulin: different method for selecting for a trait.





p3<-p[,c("TOPMEDID", "FastingInsulin")]
names(p3) <- c("TOPMEDID", "FastingInsulin1")
d<-merge(d,p3,by.x = 'ID1',by.y = 'TOPMEDID',all.x=T)

p3<-p[,c("TOPMEDID", "FastingInsulin")]
names(p3) <- c("TOPMEDID", "FastingInsulin2")
d<-merge(d,p3,by.x = 'ID2',by.y = 'TOPMEDID',all.x=T)

d$keep1_fi<-d$keep1
d$keep2_fi<-d$keep2

for (j in 1:nrow(d)){
	
		
	if(is.na(d$FastingInsulin1[j])==F&is.na(d$FastingInsulin2[j])==F){
			d$keep1_fi[j]<-ifelse(is.na(d$FastingInsulin1[j])==F,d$keep1[j],ifelse(is.na(d$FastingInsulin1[j])==T,0,NA))
			d$keep2_fi[j]<-ifelse(is.na(d$FastingInsulin2[j])==F,d$keep2[j],ifelse(is.na(d$FastingInsulin2[j])==T,0,NA))
	}
	else if(is.na(d$FastingInsulin1[j])==T&is.na(d$FastingInsulin2[j])==F)
	{
		d$keep1_fi[j]<-0
		d$keep2_fi[j]<-1
	} 
	else if(is.na(d$FastingInsulin1[j])==F&is.na(d$FastingInsulin2[j])==T)
	{
		d$keep1_fi[j]<-1
		d$keep2_fi[j]<-0
	}
	else if (is.na(d$FastingInsulin1[j])==T&is.na(d$FastingInsulin2[j])==T)
	{
		d$keep1_fi[j]<-0
		d$keep2_fi[j]<-0
	}
	else{ 
		d$keep1_fi[j]<-NA
		d$keep2_fi[j]<-NA
	}
}

table(d$keep1_fi,d$keep2_fi,useNA='always')


write.table(d,"/restricted/projectnb/glycemic/peitao/phenotype_harmonization/pooled_analysis/duplicates.txt",row.names=F,col.names=F,quote=F,sep='\t')