task getScript {
	command {
		wget https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/fitNull/genesis_nullmodel.R && head genesis_nullmodel.R
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
	File genotypefile
	String label

	File kinshipmatrix
	String phenoid # column name in other scripts

	File script

	command {
		R --vanilla --args ${phenofile} ${outcomename} ${outcometype} ${covariates} ${genotypefile} ${label} ${kinshipmatrix} ${phenoid} < ${script}
	}

	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk 20 SSD"
		memory: "20G"
	}

	output {
		File model = "${label}.Rda"
	}
}


workflow nullModel {
	File this_phenofile
	String this_outcomename
	String this_outcometype
	String this_covariates
	Array[File] this_genotypefile
	String this_label

	File this_kinshipmatrix
	String this_phenoid

	call getScript
	
	call fitNull {
            input: phenofile=this_phenofile, outcomename=this_outcomename, outcometype=this_outcometype, covariates=this_covariates, genotypefile=this_genotypefile[0], label=this_label, kinshipmatrix=this_kinshipmatrix, phenoid=this_phenoid, script=getScript.null_script
	}
}