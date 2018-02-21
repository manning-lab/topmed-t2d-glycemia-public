args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
######INPUTS#######
####read duplicates file
d<-read.table("freeze5b_duplicates.txt",header=T)
####read linker file
linker<-read.table("freeze5b_sample_annot_2017-12-01.txt",header=T)
t<-table(linker$study,linker$CENTER)
t<-prop.table(t,1)

###add center
linker<-linker[,c("sample.id","CENTER")]
names(linker)<-c("sample.id","center1")
d<-merge(d,linker,by.x="ID1",by.y="sample.id",all.x=T)
names(linker)<-c("sample.id","center2")
d<-merge(d,linker,by.x="ID2",by.y="sample.id",all.x=T)

####add call rate
q<-read.table("freeze_5b_pass_minDP10_qc.tsv",header=T)

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

 ##remove controls####
d<-d[which(!(d$study1%in%c("CONTROL"))&!(d$study2%in%c("CONTROL"))),]

#####FOR WITHIN THE SAME STUDY: Keeps monozygotic twins and chooses duplicates with highest sequencing call rate
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

  ###FOR ACROSS STUDIES: ranks duplicates according to cohort type, such as population-based or with longest follow up#####
  #ARIC>  DHS, GOLDN, GENOA, HyperGEN, Mayo_VTE, MESA, WHI, JHS, GeneSTAR
   else if ((as.character(d$study1[j])=='ARIC'&as.character(d$study2[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE", "MESA", "WHI", "JHS", "GeneSTAR"))|
            (as.character(d$study2[j])=='ARIC'&as.character(d$study1[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE", "MESA", "WHI", "JHS", "GeneSTAR")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='ARIC',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }#FHS>   COPDGene, MGH_AF, WHI
    else if ((as.character(d$study1[j])=='FHS'&as.character(d$study2[j])%in%c("COPDGene","MGH_AF","WHI"))|
             (as.character(d$study2[j])=='FHS'&as.character(d$study1[j])%in%c("COPDGene","MGH_AF","WHI")))
    {
      d$keep1[j]<-ifelse(as.character(d$study1[j])=='FHS',1,0)
      d$keep2[j]<-1-d$keep1[j]
    }#GeneSTAR>  COPDGene
   else if ((as.character(d$study1[j])=='GeneSTAR'&as.character(d$study2[j])%in%c("COPDGene"))|
            (as.character(d$study2[j])=='GeneSTAR'&as.character(d$study1[j])%in%c("COPDGene")))
   {
     d$keep1[j]<-ifelse(as.character(d$study1[j])=='GeneSTAR',1,0)
     d$keep2[j]<-1-d$keep1[j]
   }#HyperGEN>          COPDGene, DHS
  else if ((as.character(d$study1[j])=='HyperGEN'&as.character(d$study2[j])%in%c("COPDGene", "DHS"))|
           (as.character(d$study2[j])=='HyperGEN'&as.character(d$study1[j])%in%c("COPDGene", "DHS")))
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
  ##MESA>   COPDGene, DHS, GENOA, HyperGEN, GeneSTAR
  else if ((as.character(d$study1[j])=='MESA'&as.character(d$study2[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN", "GeneSTAR"))|
           (as.character(d$study2[j])=='MESA'&as.character(d$study1[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN", "GeneSTAR")))
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
   ###WHI>     COPDGene, HyperGEN, Mayo_VTE, MESA
  else if ((as.character(d$study1[j])=='WHI'&as.character(d$study2[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA"))|
           (as.character(d$study2[j])=='WHI'&as.character(d$study1[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='WHI',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
   ####CHS > MESA, WHI 
  else if ((as.character(d$study1[j])=='CHS'&as.character(d$study2[j])%in%c("MESA", "WHI"))|
          (as.character(d$study2[j])=='CHS'&as.character(d$study1[j])%in%c("MESA", "WHI")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='CHS',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
  ####MGH_AF > VU_AF, CCAF, Partners
   else if ((as.character(d$study1[j])=='MGH_AF'&as.character(d$study2[j])%in%c("CCAF", "VU_AF", "Partners"))|
          (as.character(d$study2[j])=='MGH_AF'&as.character(d$study1[j])%in%c("CCAF", "VU_AF", "Partners")))
  {
    d$keep1[j]<-ifelse(as.character(d$study1[j])=='MGH_AF',1,0)
    d$keep2[j]<-1-d$keep1[j]
  } 
  else { d$keep1[j]<-NA; d$keep2[j]<-NA}
}


d <- d[c("ID1", "ID2", "study1", "study2", "MZtwinID", "center1", "center2", "cr1", "cr2", "p1", "p2", "keep1", "keep2")]
### for reordering ID1 and ID2. Used column names for clarity but can use column numbers for simplicity. 


write.table(d,"duplicates.txt",row.names=F,col.names=T,quote=F,sep='\t')

}

if(length(args) >1) {
	stop ("Specify only one trait at a time, without spaces.")
}

if(length(args) ==1) {

p<-read.table("Pooled_Glycemic_Traits_freeze5_duplicate_ID_20180116.ped",header=T)

d<- read.table("duplicates.txt",header=T)

TEST <- as.character(args[1])

	if(TEST %in% colnames(p)==FALSE) {
	print(TEST)
	stop("Trait entered must match column name in phenotype file")
	}
	
	if(TEST %in% colnames(p)==TRUE){
		####ACROSS STUDIES: indicates which duplicate to remove based on missing trait data and cohort type 

		####TRAIT
		p3<-p[,c("TOPMEDID", "TRAIT")]
		names(p3) <- c("TOPMEDID", "TRAIT1")
		d<-merge(d,p3,by.x = 'ID1',by.y = 'TOPMEDID',all.x=T)

		p3<-p[,c("TOPMEDID", "TRAIT")]
		names(p3) <- c("TOPMEDID", "TRAIT2")
		d<-merge(d,p3,by.x = 'ID2',by.y = 'TOPMEDID',all.x=T)

		d$keep1_TRAIT<-d$keep1
		d$keep2_TRAIT<-d$keep2

		d$reason_TRAIT <- NA
	
		
		for (j in 1:nrow(d)){
			
				
			if(is.na(d$TRAIT1[j])==F&is.na(d$TRAIT2[j])==F&d$TRAIT1[j]>=0&d$TRAIT2[j]>=0){
					d$keep1_TRAIT[j]<-ifelse(is.na(d$TRAIT1[j])==F,d$keep1[j],ifelse(is.na(d$TRAIT1[j])==T,0,NA))
					d$keep2_TRAIT[j]<-ifelse(is.na(d$TRAIT2[j])==F,d$keep2[j],ifelse(is.na(d$TRAIT2[j])==T,0,NA))
				##for additional output files
				if(is.na(d$keep1[j])==F&is.na(d$keep2[j])==F){
					if(d$keep1[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
						d$reason_TRAIT[j] <- "used duplicate from older study"
					}
					else if(d$keep2[j]==1&as.character(d$study1[j])!=as.character(d$study2[j])){
						d$reason_TRAIT[j] <- "used duplicate from population-based and/or older study"
					}
					else if(d$keep1[j]==1&d$keep2[j]==1){
						d$reason_TRAIT[j] <- "kept both; monozygotic twins" 
					}
					else if(d$keep1[j]==1&d$keep2[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
						d$reason_TRAIT[j] <- "kept duplicate with higher callrate"
					}
					else if(d$keep2[j]==1&d$keep1[j]==0&as.character(d$study1[j])==as.character(d$study2[j])){
						d$reason_TRAIT[j] <- "kept duplicate with higher callrate"
					}
				}
			}
			else if(is.na(d$TRAIT1[j])==T&is.na(d$TRAIT2[j])==F)
			{
				d$keep1_TRAIT[j]<-0
				d$keep2_TRAIT[j]<-1
				d$reason_TRAIT[j] <- "missing trait data"
			} 
			else if(is.na(d$TRAIT1[j])==F&is.na(d$TRAIT2[j])==T)
			{
				d$keep1_TRAIT[j]<-1
				d$keep2_TRAIT[j]<-0
				d$reason_TRAIT[j] <- "missing trait data"
			}
			else if (is.na(d$TRAIT1[j])==T&is.na(d$TRAIT2[j])==T)
			{
				d$keep1_TRAIT[j]<-0
				d$keep2_TRAIT[j]<-0
				d$reason_TRAIT[j] <- "missing trait data"
			}
			else if (d$TRAIT1[j]<0&d$TRAIT2[j]>=0) ##for binary trait if negative numbers are used for NA
			{
				d$keep1_TRAIT[j]<-0
				d$keep2_TRAIT[j]<-1
				d$reason_TRAIT[j] <- "missing trait data"
			}
			else if (d$TRAIT1[j]>=0&d$TRAIT2[j]<0) ##for binary trait if negative numbers are used for NA
			{
				d$keep1_TRAIT[j]<-1
				d$keep2_TRAIT[j]<-0
				d$reason_TRAIT[j] <- "missing trait data"
			}
			else if(d$TRAIT1[j]<0&d$TRAIT2[j]<0)  ##for binary trait if negative numbers are used for NA
			{
				d$keep1_TRAIT[j]<-0
				d$keep2_TRAIT[j]<-0
				d$reason_TRAIT[j] <- "missing trait data"
			}
			else{ 
				d$keep1_TRAIT[j]<-NA
				d$keep2_TRAIT[j]<-NA
			}
		}

		print(table(d$keep1_TRAIT,d$keep2_TRAIT,useNA='always'))
		print("Error if NAs occur")

		
		####OUTPUT: to subset for TRAIT and to remove duplicates as determined by algorithm####

		not_keep1<-subset(d, keep1_TRAIT==0)
		not_keep1<-not_keep1[c("ID1","keep1_TRAIT")]

		not_keep2<-subset(d, keep2_TRAIT==0)
		not_keep2<-not_keep2[c("ID2","keep2_TRAIT")]

		names<-c("ID","KEEP_TRAIT")
		not_keep<-data.frame()
		colnames(not_keep1)<-names
		colnames(not_keep2)<-names
		not_keep <-rbind(not_keep1, not_keep2)

	
		p<-merge(p,not_keep,by.x="TOPMEDID",by.y="ID",all.x=T)

		p$KEEP_TRAIT[is.na(p$KEEP_TRAIT)&is.na(p$TRAIT)] <- 0
		p$KEEP_TRAIT[is.na(p$KEEP_TRAIT)&p$TRAIT<0] <- 0 ##if trait is binary and uses negative number for NA
		p$KEEP_TRAIT[is.na(p$KEEP_TRAIT)&!is.na(p$TRAIT)] <- 1

		final <- subset(p, (!is.na(p$TRAIT)))
		final <-subset(final, KEEP_TRAIT==1)

		write.table(final,"TRAIT/removed_duplicates_TRAIT.ped",row.names=F,col.names=T,quote=F,sep='\t')


		## for pooled phenotypes file with additional columns (for inspection purposes)

		write.table(p,"./TRAIT/all_TRAIT.ped",row.names=F,col.names=T,quote=F,sep='\t')

		###for the duplicates.txt file with additional trait-specific keep columns, if desired: 
		d <- d[c("ID1", "ID2", "study1", "study2", "MZtwinID", "center1", "center2", "cr1", "cr2", "p1", "p2", "keep1", "keep2", "TRAIT1", "TRAIT2", "keep1_TRAIT", "keep2_TRAIT", "reason_TRAIT")]
		
		write.table(d,"TRAIT/duplicates_TRAIT.txt",row.names=F,col.names=T,quote=F,sep='\t')

		#### for lists of duplicates that were removed####

		names <- c("ID", "reason")

		removed1_TRAIT <- subset(d, keep1_TRAIT==0)
		removed1_TRAIT <- removed1_TRAIT[c("ID1", "reason_TRAIT")]
		colnames(removed1_TRAIT) <- names

		removed2_TRAIT <-subset(d, keep2_TRAIT==0)
		removed2_TRAIT <- removed2_TRAIT[c("ID2", "reason_TRAIT")]
		colnames(removed2_TRAIT) <- names

		removed_TRAIT <- rbind(removed1_TRAIT, removed2_TRAIT)

		write.table(removed_TRAIT, "TRAIT/removed_ID_TRAIT.txt", row.names=F, col.names=T, quote=F, sep='\t')
	}
}
