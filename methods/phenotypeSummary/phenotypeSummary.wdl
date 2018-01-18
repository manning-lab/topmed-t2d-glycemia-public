task getScript {
	command {
		wget https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/phenotype_summary/methods/phenotypeSummary/phenotypeSummary.R
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "phenotypeSummary.R"
	}
}

task summary {
	File ped_file
	String outcome 
	String covars 
	String label 
	String cohort_column 

	File script
	Int memory
	Int disk

	command {
		R --vanilla --args ${ped_file} ${outcome} ${covars} ${label} ${cohort_column} < ${script}
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


workflow phenotypeSummary {
	File this_ped_file
	String this_outcome
	String this_covars
	String this_label
	String this_cohort_column

	Int this_memory
	Int this_disk

	call getScript
	
	call summary {
            input: ped_file = this_ped_file, outcome = this_outcome, covars = this_covars, label = this_label, cohort_column = this_cohort_column, script = getScript.script, memory = this_memory, disk = this_disk
	}
}