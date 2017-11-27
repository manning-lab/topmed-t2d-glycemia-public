task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/methods/dataModel/cohortToSamples.py"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File outscript = "cohortToSamples.py"
	}
}

task runScript {
	String cohorts
	String label
	File script

	command {
		python ${script} --cohorts ${cohorts} --outfile_pref ${label}
	}

	runtime {
		docker: "broadgdac/fiss@sha256:a65324c8cf1edc769bee3195c798468defacefece3a3d551143706cd412e4c39"
		disks: "local-disk 10 SSD"
		memory: "2G"
	}

	output { 
		File out_file = "${label}.txt"
	}
}

workflow w {
	String these_cohorts
	String this_label
	
	call getScript

	call subset { 
		input: cohorts=this_cohorts, label=this_label, script=getScript.outscript
	}

	output {
		File sample_list = runScript.out_file
	}
}