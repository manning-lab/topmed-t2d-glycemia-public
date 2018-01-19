# topmed-t2d-glycemia-public
Scripts, docker images and pipelines for analysis by the TOPMed Diabetes working group


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


Table of Cohort Rules:
[[ FILL THIS IN ]]

