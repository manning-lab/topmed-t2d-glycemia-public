task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/fitNull/genesis_nullmodel.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/agg_assoc/workflows/aggregateAssociation/aggregateAssociation.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/agg_assoc/workflows/aggregateAssociation/aggregateSummary.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File null_script = "genesis_nullmodel.R"
		File assoc_script = "aggregateAssociation.R"
		File summary_script = "aggregateSummary.R"
	}
}

task fitNull {
	File genotype_file
	File? phenotype_file
	String? outcome_name
	String? outcome_type
	String? covariates_string
	String? ivars_string
	File? sample_file
	String label
	File? kinship_matrix
	String? id_col
	File script

	Int memory
	Int disk

	command {
		R --vanilla --args ${genotype_file} ${phenotype_file} ${outcome_name} ${outcome_type} ${default="NA" covariates_string} "NA" ${default="NA" ivars_string} ${sample_file} ${label} ${kinship_matrix} ${id_col} < ${script}
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

task aggAssocTest {
	File gds_file
	File null_file
	File group_file
	String label
	String? test
	String? pval
	String? weights
	
	
	File script

	Int memory
	Int disk

	command {
		R --vanilla --args ${gds_file} ${null_file} ${group_file} ${label} ${default="SKAT" test} ${default="kuonen" pval} ${default="1,25" weights} < ${script} 
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
	String label
	Array[File] assoc

	File summaryScript

	Int memory
	Int disk

	command {
		R --vanilla --args ${label} ${sep = ',' assoc} < ${summaryScript}	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
  	    disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File mhplot = "${label}_association_plots.png"
		File assoc_res = "${label}.groupAssoc.csv"
	}
}

workflow group_assoc_wf {
	# fitNull inputs
	File? this_phenotype_file
	String? this_outcome_name
	String? this_outcome_type
	String? this_covariates_string
	String? this_ivars_string
	File? this_sample_file
	String this_label
	File? this_kinship_matrix
	String? this_id_col

	# aggAssocTest inputs
	Array[File] these_gds_files
	File? this_null_file
	Array[File] these_group_files
	String? this_test
	String? this_pval
	String? this_weights
	
	# other inputs
	Int this_memory
	Int this_disk
	
	# gds and group files must be in the same order, one group file per gds file
	Array[Pair[File,File]] these_gds_groups = zip(these_gds_files, these_group_files)

	call getScript

	File null_genotype_file = these_gds_files[0]
	Boolean have_null = defined(this_null_file)

	if (!have_null) {

		call fitNull {
				input: genotype_file = null_genotype_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, ivars_string = this_ivars_string, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, script = getScript.null_script, memory = this_memory, disk = this_disk
			}

		scatter(this_gds_group in these_gds_groups) {
			
			call aggAssocTest {
				input: gds_file = this_gds_group.left, null_file = fitNull.model, group_file = this_gds_group.right, label = this_label, test = this_test, pval = this_pval, weights = this_weights, memory = this_memory, disk = this_disk, script = getScript.assoc_script
			}
		}
	

		call summary {
			input: assoc = aggAssocTest.assoc, label = this_label, memory = this_memory, disk = this_disk, summaryScript = getScript.summary_script
		}
	}

	if (have_null) {

		scatter(this_gds_group in these_gds_groups) {
			
			call aggAssocTest as aggAssocTest_null_in {
				input: gds_file = this_gds_group.left, null_file = this_null_file, group_file = this_gds_group.right, label = this_label, test = this_test, pval = this_pval, weights = this_weights, memory = this_memory, disk = this_disk, script = getScript.assoc_script
			}
		}
	

		call summary as summary_null_in {
			input: assoc = aggAssocTest.assoc, label = this_label, memory = this_memory, disk = this_disk, summaryScript = getScript.summary_script
		}
	}
}