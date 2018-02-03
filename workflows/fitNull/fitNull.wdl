task getScript {
	command {
		wget https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/fitNull/genesis_nullmodel.R
		wget https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/methods/phenotypeSummary/phenotypeSummary.R
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "genesis_nullmodel.R"
		File summary_script = "phenotypeSummary.R"
	}
}

task fitNull {
	File genotype_file
	File phenotype_file
	String outcome_name
	String outcome_type
	String covariates_string
	File sample_file
	String label
	File kinship_matrix
	String id_col
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

task summary {
	File phenotype_file
	String outcome_name
	String covariates_string
	String label
	String cohort_column

	File script
	Int memory
	Int disk

	command {
		R --vanilla --args ${phenotype_file} ${outcome_name} ${covariates_string} ${label} ${cohort_column} < ${script}
	}

	runtime {
		docker: "tmajarian/r-ggplot@sha256:e4a9a53a49faf7a4d0318e56db14b85284b8d1a1bb5ceee8b04a9979efa8d4ca"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File plots = "${label}_plots.pdf"
		File stats = "${label}_stats.csv"
	}
}

workflow nullModel {
	# fitNull inputs
	Array[File] these_genotype_file
	File this_phenotype_file
	String this_outcome_name
	String this_outcome_type
	String this_covariates_string
	File this_sample_file
	String this_label
	File this_kinship_matrix
	String this_id_col
	
	# summary inputs
	String this_cohort_column

	# other workflow inputs
	Int this_memory
	Int this_disk
	
	File this_genotype_file = these_genotype_file[0]

	call getScript
	
	call fitNull {
            input: genotype_file = this_genotype_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, script = getScript.script, memory = this_memory, disk = this_disk
	}

	call summary {
		input: phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, covariates_string = this_covariates_string, label = this_label, cohort_column = this_cohort_column, script = getScript.summary_script, memory = this_memory, disk = this_disk
	}

	output {
		File out_file = fitNull.model
		File out_plots = summary.plots
		File out_stats = summary.stats
	}
}
