# Fit Null Model using Genesis

## Description 

This workflow generates a null model for use in association testing. The primary code is written in R using the GENESIS package for model fitting.

### Authors

This workflow is produced and maintained by the [Manning Lab](https://manning-lab.github.io/). It has been heavily influenced by code written by Jen Brody for the DNANexus platform. Contributing authors include:

* Tim Majarian (tmajaria@broadinstitute.org)
* Alisa Manning (amanning@broadinstitute.org)

## Dependencies

### Workflow execution

* [WDL](https://software.broadinstitute.org/wdl/documentation/quickstart)
* [Chromwell](http://cromwell.readthedocs.io/en/develop/)

### R packages

* [GENESIS](https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html)
* [GWASTools](https://www.bioconductor.org/packages/release/bioc/html/GWASTools.html)
* [SeqArray](https://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)
* [SeqVarTools](https://www.bioconductor.org/packages/release/bioc/html/SeqVarTools.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)

## Main Functions

### fitNull

This function generates a null model to be used in association testing in Genesis.

Inputs:
* genotype_files : genotype data for all samples (array of VCF or GDS file)
* phenotype_file : phenotype data for all samples to be included in analysis (CSV or TSV file)
* outcome_name : the outcome to be tested (string)
* outcome_type : the type of outcome being tested (dichotomous or continuous)
* covariates_string : covariates to condition on in linear mixed modeling (comma separated string, default = '')
* sample_file : a file containing a list of sample ids (matching the genotype and phenotype files) to be included, one per line (.txt)
* label : prefix for output filename (string)
* kinship_file : relatedness measures for all samples (CSV or TSV file)
* id_col : column name of id column in phenotype file (string)

Outputs:
* model : generated null model (.RDa)

### summary

This function generates summary data for the input phenotype file and covariates. Plots and statistics for covariates subset by cohort are produced.

Inputs:
* phenotype_file : phenotype data for all samples to be included in analysis (CSV or TSV file)
* outcome_name : the outcome to be tested (string)
* covariates_string : covariates to condition on in linear mixed modeling (comma separated string, default = '')
* label : prefix for output filename (string)
* cohort_column : column header in phenotype file for cohort annotations (string)

Outputs:
* plots : a pdf of box and whisker plots for all covariates subset by cohort (.pdf)
* stats : a csv of mean, median, standard deviation, min, and max for each covariate subset by cohort and outcome (.csv)

## Other workflow inputs

* this_memory : amount of memory in GB for each execution of a task (int)
* this_disk : amount of disk space in GB to allot for each execution of a task (int)



