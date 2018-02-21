# Glycemic Traits Harmonization readme

## Prereq's

1. Put all of the study-specific phenotype files, below, in a single directory. Note the absolutely path to that directory.

```bash
freeze5b_sample_annot_2017-10-16.txt
FHS_glycemicTraits_20170214_CS.ped
share_unrel_comb_01012017.csv
JHS_glycemictraits_29Mar2017.ped
SAS_LIULIN_20170212_glycemic_traits.ped
CFS_12May17_glycemic_traits_Freeze4.ped
Amish_HuichunXu_May15th2017_glycemic_traits.ped
GenSalt_EA_XuenanMi_20170614_glycemic_traits.ped
ARIC_EU_NDAYA_20161220_glycemic_traits_ped.csv
ARIC_AA_NDAYA_20161220_glycemic_traits_ped.csv
WHI_HA_MP_20170726_glycemic_traits.ped
WHI_EU_MP_20170726_glycemic_traits.ped
WHI_AS_MP_20170726_glycemic_traits.ped
WHI_AA_MP_20170726_glycemic_traits.ped
GeneSTAR_AA_YANEK_20170720_glycemic_traits.ped
GeneSTAR_EU_YANEK_20170720_glycemic_traits.ped
CHS_FLOYD_201708018_glycemic_traits.PED
MESA_HA_ABIGAILBALDRIDGE_04JAN17_glycemic_traits_sidno.csv
freeze5b_sample_annot_2017-10-16.txt
MESA_SA_ABIGAILBALDRIDGE_04JAN17_glycemic_traits_sidno.csv
MESA_AA_ABIGAILBALDRIDGE_04JAN17_glycemic_traits_sidno.csv
GENOA_AA_BIELAK_20171124_glycemic_traits.ped
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