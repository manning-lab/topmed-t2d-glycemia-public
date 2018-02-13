# T2D Harmonization readme

## Prereq's

1. Put all of the study-specific phenotype files, below, in a single directory. Note the absolutely path to that directory.

```bash
freeze5b_sample_annot_2017-12-01.txt
MGH_Parent_dbGaP_SubjectPhenotypesDS_v3.txt
MGH_Parent_PedigreeDS_v1.txt
phs001032.v3.pht005675.v3.p2.c1.TOPMed_WGS_VUH_AF_Subject_Phenotypes.GRU-IRB.txt
phs001032.v3.pht007135.v1.p2.TOPMed_WGS_VUH_AF_Pedigree.MULTI.txt
phs001024.v2.pht005693.v1.p1.c1.TOPMed_WGS_Partners_AFGen_Subject_Phenotypes.HMB.txtsep=/)
VUH_Ablation_dbGaP_SubjectPhenotypesDS_v2.txt
CCAF_dbGaP_SubjectPhenotypesDS_v5.txt
HVH_FLOYD_20170818_T2DPhenotype.PED
ARIC_AA_NMAYA_20161220_T2D.ped
ARIC_EU_NMAYA_20161220_T2D.ped
phs000179.v5.pht002239.v4.p2.c1.COPDGene_Subject_Phenotypes.HMB.txt
phs000179.v5.pht002239.v4.p2.c2.COPDGene_Subject_Phenotypes.DS-CS-RD.txt
Amish_HuichunXu_May152017_T2D.ped
CFS_12May17_T2D_Freeze4.ped
SAS_LIULIN_20170212_T2D.ped
JHS_T2Dincidence_16Aug2017.ped
FHS_T2D_20170214_CS.ped
WHI_EU_MP_20170726_T2D.ped
WHI_AA_MP_20170726_T2D.ped
WHI_AS_MP_20170726_T2D.ped
WHI_HA_MP_20170726_T2D.ped
GenSalt_EA_XuenanMi_20170614_T2D.ped
GeneSTAR_AA_YANEK_20170720_T2D.ped
GeneSTAR_EU_YANEK_20170720_T2D.ped
MESA_AA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv
MESA_EU_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv
MESA_HA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv
MESA_SA_ABIGAILBALDRIDGE_04JAN17_T2D_sidno.csv
MESA_Family_AA_YAO_20171031_T2D_sidno.txt
CHS_FLOYD_201708018_T2DPhenotype.PED
GENOA_AA_BIELAK_20171030_T2D.ped
2a_dbGaP_SubjectPhenotypesDS_AACAC.txt
6a_dbGaP_PedigreeDS_AACAC_revised072817.txt
```

2. Make sure your github repo is up-to-date and that you're on the correct branch.

```bash
cd /my/local/git/repo/path/topmed-t2d-glycemia-public/ 
git stash 
git fetch
git checkout glycemic_traits_harmonization
git pull origin glycemic_traits_harmonization
```

3. Note the absolute path to the file `Harmonization.19JAN2017.GitHub.SourceFiles.R`

## How to Run

### From RStudio

Open R studio and the script `Harmonization.19JAN2017.GitHub.R`. Set these two variables:
```R
f.dir <- /phenotype/file/directory/ # from step 1 in Prereq's
source.file <- /my/local/git/repo/path/topmed-t2d-glycemia-public/methods/traitHarmonization/Harmonization.19JAN2017.GitHub.SourceFiles.R # from step 3 of Prereq's
```

*You may also need to comment out lines 15 and 16* 

### From the terminal

```bash
R --vanilla --args /phenotype/file/directory/ /my/local/git/repo/path/topmed-t2d-glycemia-public/methods/traitHarmonization/Harmonization.19JAN2017.GitHub.SourceFiles.R < Harmonization.19JAN2017.GitHub.R`
```