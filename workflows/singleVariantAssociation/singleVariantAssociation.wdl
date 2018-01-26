task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/issue_33/workflows/singleVariantAssociation/association.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/fitNull/genesis_nullmodel.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/singleVariantAssociation/summary.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/issue_33/workflows/singleVariantAssociation/preprocess_conditional.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File assoc_script = "association.R"
		File null_script = "genesis_nullmodel.R"
		File summary_script = "summary.R"
		File conditional_script = "preprocess_conditional.R"
	}
}

task conditionalPhenotype {
	Array[File] genotype_files
	File? phenotype_file
	String? id_col
	File? sample_file
	String? snps
	String label

	File script

	Int disk

	command {
		R --vanilla --args ${sep="," genotype_files} ${phenotype_file} ${id_col} ${sample_file} ${snps} ${label} < ${script}
	}

	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk ${disk} SSD"
		memory: "10G"
	}

	output {
		File new_phenotype_file = "${label}_phenotypes.csv"
		File alt_ref = "${label}_alleles.txt"
	}
}

task fitNull {
	File genotype_file
	File phenotype_file
	String outcome_name
	String outcome_type
	String covariates_string
	String? conditional_string
	File sample_file
	String label
	File kinship_matrix
	String id_col
	File script

	Int memory
	Int disk

	command {
		R --vanilla --args ${genotype_file} ${phenotype_file} ${outcome_name} ${outcome_type} ${covariates_string} ${default="NA" conditional_string} ${sample_file} ${label} ${kinship_matrix} ${id_col} < ${script}
	}

	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File model = "${label}_null.RDa"
	}
}

task assocTest {

	File gds_file
	File? null_file
	String label
	String? test
	Int? mac
	File script

	Int memory
	Int disk

	command {
		R --vanilla --args ${gds_file} ${null_file} ${label} ${default="Score" test} ${default="5" mac} < ${script} 
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File assoc = "${label}.assoc.RData"
	}
}

task summary {
	String pval
	Float? pval_threshold
	String label
	Array[File] assoc
	File script

	Int memory
	Int disk

	command {
		R --vanilla --args ${pval} ${default="0.0001" pval_threshold} ${label} ${sep = ',' assoc} < ${script}	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
  	    disks: "local-disk ${disk} SSD"
        memory: "${memory}G"
	}

	output {
		File allassoccsv = "${label}.assoc.csv"
		File topassoccsv = "${label}.topassoc.csv"
		File plots = "${label}_association_plots.png"
	}
}

workflow w_assocTest {
	# conditionalPhenotype inputs
	String? these_snps


	# fitNull inputs
	Array[File] these_genotype_files
	File? this_phenotype_file
	String? this_outcome_name
	String? this_outcome_type
	String? this_covariates_string
	File? this_sample_file
	String this_label
	File? this_kinship_matrix
	String? this_id_col
	
	# assocTest inputs
	File? this_null_file
	String? this_test
	Int? this_mac

	# summary inputs
	String this_pval
	Float? this_pval_threshold	

	# inputs to all
	Int this_memory
	Int this_disk

	call getScript

	File null_genotype_file = these_genotype_files[0]

	Boolean need_null = defined(this_null_file)

	if(defined(these_snps)) {

		call conditionalPhenotype {
			input: genotype_files = these_genotype_files, phenotype_file = this_phenotype_file, id_col = this_id_col, sample_file = this_sample_file, snps = these_snps, label = this_label, script = getScript.conditional_script, disk = this_disk
		}
		
		call fitNull as fitNullConditional {
			input: genotype_file = null_genotype_file, phenotype_file = conditionalPhenotype.new_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, conditional_string = these_snps, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, script = getScript.null_script, memory = this_memory, disk = this_disk
		}

		scatter(this_genotype_file in these_genotype_files) {
		
			call assocTest as assocTestConditional {
				input: gds_file = this_genotype_file, null_file = fitNullConditional.model, label = this_label, test = this_test, mac = this_mac, script = getScript.assoc_script, memory = this_memory, disk = this_disk
			}
		}

		call summary as summaryConditional {
			input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = assocTestConditional.assoc, script = getScript.summary_script, memory = this_memory, disk = this_disk
		}
	}

	if (!defined(these_snps)) {
	
		if(!need_null) {
			
			call fitNull {
				input: genotype_file = null_genotype_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, script = getScript.null_script, memory = this_memory, disk = this_disk
			}

			scatter(this_genotype_file in these_genotype_files) {
			
				call assocTest {
					input: gds_file = this_genotype_file, null_file = fitNull.model, label = this_label, test = this_test, mac = this_mac, script = getScript.assoc_script, memory = this_memory, disk = this_disk
				}
			}

			call summary {
				input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = assocTest.assoc, script = getScript.summary_script, memory = this_memory, disk = this_disk
			}

		} 

		if(need_null) {

			scatter(this_genotype_file in these_genotype_files) {
			
				call assocTest as assocNull {
					input: gds_file = this_genotype_file, null_file = this_null_file, label = this_label, test = this_test, mac = this_mac, script = getScript.assoc_script, memory = this_memory, disk = this_disk
				}
			}

			call summary as summaryNull {
				input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = assocNull.assoc, script = getScript.summary_script, memory = this_memory, disk = this_disk
			}
		}
	}
}
