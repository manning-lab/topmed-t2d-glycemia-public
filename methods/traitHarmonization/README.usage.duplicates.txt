READ ME for duplicates.sh and duplicates.r

DESCRIPTION

Duplicates.sh is a shell script that manages and runs the duplicates.r script. 

Duplicates.r is an r script that performs two tasks:

1. selects duplicates for removal based on an algorithm (in Readme.me file)
2. when a trait is specified, subsets a phenotype file by removing the duplicates (based on trait data and selections in step1) and by removing subjects with missing trait data. 

*NOTE: trait data ultimately determines which duplicate will be kept. If only one of the duplicates have trait data, then that one will be kept.
Cohort characteristics and call rates are used to break "ties" when both duplicates have trait data.  

INPUTS
List of duplicates = freeze5b_duplicates.txt
Center call rates = freeze_5b_pass_minDP10_qc.tsv
Linker file = freeze5b_sample_annot_2017-12-01.txt
Phenotype file = Pooled_Glycemic_Traits_freeze5_duplicate_ID_20180116.ped
*NOTE: inputs should be in same directory as scripts. 

USAGE

-Both scripts must be in the same directory

-Duplicates.sh takes a trait as an argument. If no argument is specified, it will only run the first part of duplicates.r script. 
If you already ran the duplicates.r script and do not specify a trait, than the script will exit.

-The trait entered as an argument must match how it is written in the column header in the phenotype file. The script will tell you if there is a typo. 

-If you specify a trait and never ran the duplicates.r script, it will run both parts, giving as outputs:
   In working directory:
	1. A duplicates.txt file with columns indicating which ID to keep based on call rate and cohort characteristics 
   In trait subdirectory of working directory (created by script):
	2. A phenotype file with duplicates and subjects missing trait data removed
	3. A text file indicating which duplicates were removed a why
	4. an r script for that trait for future use 
     For inspection/record keeping:
	5. a phenotype file with all subjects (including duplicates), with a "keep" column indicating which subjects are flagged for removal
	6. a duplicates_trait.txt file with columns to keep based on missing trait data, call rate, and cohort characteristics.

If you specify a trait and already ran the duplicates.r script, then it will return the same outputs except #1 (since it already exits)

The duplicates.r script can be run on its own without running the shell script. It will only give you output#1 and will not subset for trait data.  

