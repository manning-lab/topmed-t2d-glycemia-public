task getScript {
	command {
		wget https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/fitNull/genesis_nullmodel.R
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File null_script = "genesis_nullmodel.R"
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


workflow nullModel {
	Array[File] this_genotype_file
	File this_phenotype_file
	String this_outcome_name
	String this_outcome_type
	String this_covariates_string
	File this_sample_file
	String this_label
	File this_kinship_matrix
	String this_id_col
	
	Int this_memory
	Int this_disk
	

	call getScript
	
	call fitNull {
            input: genotype_file = this_genotype_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, script = getScript.script, memory = this_memory, disk = this_disk
	}

	output {
		File out_file = fitNull.model
	}
}