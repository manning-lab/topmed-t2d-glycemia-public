# Single Variant Association -- Linear Mixed Models ###

## Description 

This workflow performs a single variant association analysis of genotype data with a single phenotype using linear mixed models. The primary code is written in R using the GENESIS package for model fitting. See the accompanying README for full description and input file format.

## Main Functions

### fitNull

This function generates a null model to be used in association testing using Genesis. 

Inputs:
* genotype_files : genotype data for all samples (array of VCF or GDS file)
* phenotype_file : phenotype data for all samples to be included in analysis (CSV or TSV file)
* outcome : the outcome to be tested (string)
* outcome_type : the type of outcome being tested (dichotomous or continuous)
* covariates_string : covariates to condition on in linear mixed modeling (comma separated string, default = '')
* sample_file : a file containing a list of sample ids (matching the genotype and phenotype files) to be included, one per line (.txt)
* label : prefix for output filename (string)
* kinship_file : relatedness measures for all samples (CSV or TSV file)
* id_col : column name of id column in phenotype file (string)

Outputs:
* model : generated null model (.RDa)
