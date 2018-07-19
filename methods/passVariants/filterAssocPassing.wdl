# getPassVariants
task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/methods/passVariants/filterAssocPassing.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "filterAssocPassing.R"
	}
}


task filterPass {
	File res_file
	File pass_file
	File script

	command {
		R --vanilla --args ${res_file} ${pass_file} < ${script}
	}

	runtime {
		docker: "manninglab/singlevariantassociation:latest"
		disks: "local-disk 100 SSD"
		memory: "15G"
		bootDiskSizeGb: 20
	}

	output {
		File varlist = "${res_file}"
	}
}

workflow filterVar {
	File this_res_file
	File this_pass_file

	call getScript
	
	call filterPass {
		input: res_file = this_res_file, pass_file = this_pass_file, script = getScript.script
	}
}