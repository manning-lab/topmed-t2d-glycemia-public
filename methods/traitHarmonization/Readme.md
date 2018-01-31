
# Algorithm for duplicates.R file

For every genotype file freeze (for example, freeze 5b), we have the following algorithm for determining subject_id / NWD mapping for each paper (T2D, FG/FI, HbA1c)

Inputs:
- File with NWD ID sequencing call rate and sequencing center
- File with NWD and outcome data <- trait file from which to determine if individual has missing trait data
- string specifying which trait column to use in the trait file
- Duplicates file released by TOPMed DCC
- Linker file released by TOPMed DCC

Algorithm
- Within Study
  - Keep monozygotic twins
  - If duplicates are at different sequencing centers, choose NWD by sequencing center for cohort
  - Choose duplicate with highest sequencing call rate
- Across Studies
  - Remove missing trait data, and identify remaining duplicates
  - Choose duplicates according to cohort type:
  - Population cohort > ascertained cohort
  - Choose duplicate from cohort with longest follow-up period (# years follow-up for entire cohort;  not by individual)
  - If there are any remaining duplicates, add a rule to the algorithm

iii)     Remove intentional duplicates

Output:
- A subject_id -> NWD mapping file for each trait file

## Cohort Rules:

•	ARIC > DHS, GOLDN, GENOA, HyperGEN, Mayo_VTE, MESA, WHI, JHS, GeneSTAR
•	FHS > COPDGene, MGH_AF, WHI
•	GeneSTAR > COPDGene
•	HyperGEN > COPDGene, DHS
•	JHS > GENOA
•	MESA >COPDGene, DHS, GENOA, HyperGEN, GeneSTAR
•	VU_AF > WGHS
•	WH I> COPDGene, HyperGEN, Mayo_VTE, MESA
•	CHS > MESA, WHI 
•	MGH_AF > VU_AF, CCAF, Partners

Name		Type				Year started
ARIC		Population-based	1988
FHS			Population-based	1948
GeneSTAR	Family-based (heart disease and stroke)	1982
HyperGEN	2000 subjects enriched for hypertension, 800 random subjects	1995
JHS			Population-based and family-based	1998
MESA		Population-based	2000
VU_AF		Case-ascertained 	2010
WHI			Population-based	1991
CHS			Population-based	1989
MGH_AF		Case-ascertained	2001
DHS			Family-based, diabetes enriched	2007
GOLDN		Population based	~2010
GENOA		Family-based	1995
Mayo_VTE	cases with VTE from other studies, i.e. ARIC	2013
COPDGene	Only smokers 	1994
WGHS		Case-ascertained	1992
Partners	Biorepository for Case-ascertained cohorts	
CCAF		Case-ascertained	2005


