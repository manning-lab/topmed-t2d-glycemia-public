# getPassVariants
task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/methods/passVariants/getPassVariants.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "getPassVariants.R"
	}
}


task getPass {
	File gds_file
	File script

	command {
		R --vanilla --args ${gds_file} < ${script}
	}

	runtime {
		docker: "manninglab/singlevariantassociation:latest"
		disks: "local-disk 100 SSD"
		memory: "20G"
		bootDiskSizeGb: 20
	}

	output {
		File varlist = select_first(glob("freeze.5b.passing.variants.*.minDP10.csv"))
	}
}

task combine {
	Array[File] ids

	command {
		cat ${sep=" " ids} >> freeze.5b.passing.variants.all.minDP10.csv
	}

	runtime {
		docker: "manninglab/singlevariantassociation:latest"
		disks: "local-disk 50 SSD"
		memory: "5G"
		bootDiskSizeGb: 20
	}

	output {
		File allvar = "freeze.5b.passing.variants.all.minDP10.csv"
	}
	
}

workflow passVar {
	Array[File] these_gds_files

	call getScript

	scatter(this_gds_file in these_gds_files) {
		call getPass {
			input: gds_file = this_gds_file, script = getScript.script
		}
	}

	call combine {
		input: ids = getPass.varlist
	}

	output {
		File all_variants = combine.allvar
	}
}