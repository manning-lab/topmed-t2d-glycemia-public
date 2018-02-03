task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/singleVariantAssociation/association.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/fitNull/genesis_nullmodel.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/singleVariantAssociation/summary.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File assoc_script = "association.R"
		File null_script = "genesis_nullmodel.R"
		File summary_script = "summary.R"
	}
}

task fitNull {
	File genotype_file
	File? phenotype_file
	String? outcome_name
	String? outcome_type
	String? covariates_string
	File? sample_file
	String label
	File? kinship_matrix
	String? id_col
	File script

	Int memory
	Int disk

	command {
		R --vanilla --args ${genotype_file} ${phenotype_file} ${outcome_name} ${outcome_type} ${covariates_string} ${sample_file} ${label} ${kinship_matrix} ${id_col} < ${script}
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
	
	if(!defined(this_null_file)) {

		File null_genotype_file = these_genotype_files[0]
		
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

	if(defined(this_null_file)) {

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
