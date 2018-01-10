# This is a ported over copy of the DNANexus app genesis_v0.7 written by Jen Brody 

task getScript {
	command {
		wget https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/dnanexus/genesis_nullmodel_v0.2/genesis_nullmodel.R
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File null_script = "genesis_nullmodel.R"
	}
}

task fitNull {
	File phenofile
	String outcomename
	String outcometype
	String covariates
	File genotype
	String output_file
	File kinshipmatrix
	String phenoid
	String? test_stat
	String? conditional 
	String? het_vars

	File script
	Int disk
	Int memory

	command {
		R --vanilla --args ${phenofile} ${outcomename} ${outcometype} ${covariates} ${genotype} ${output_file} ${kinshipmatrix} ${phenoid} ${default="Score" test_stat} ${default="NA" conditional} ${default="NA" het_vars} < ${script}
	}

	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File model = "${output_file}.RDa"
	}
}


workflow nullModel {
	File this_phenofile
	String this_outcomename
	String this_outcometype
	String this_covariates
	File this_genotype
	String this_output_file
	File this_kinshipmatrix
	String this_phenoid
	String? this_test_stat
	String? this_conditional
	String? this_het_vars

	Int this_disk
	Int this_memory

	call getScript
	
	call fitNull {
            input: phenofile = this_phenofile, outcomename = this_outcomename, outcometype = this_outcometype, covariates = this_covariates, genotype = this_genotype, output_file = this_output_file, kinshipmatrix = this_kinshipmatrix, phenoid = this_phenoid, test_stat = this_test_stat, conditional = this_conditional, het_vars = this_het_vars, script=getScript.null_script, disk = this_disk, memory=this_memory
	}
}