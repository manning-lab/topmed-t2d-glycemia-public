task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/methods/fileSubset/subsetGdsBySample.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File outscript = "subsetGdsBySample.R"
	}
}

task runScript {
	File gds_in
	Array[String] subset_ids

	String outbase = basename(gds_in,".gds")
	String outfile = outbase + "_subset.gds"

	File script
	Float memory
	Int disksize
	
	command {
		R --vanilla --args ${gds_in} ${sep= "," subset_ids} ${outfile} < ${script}
	}
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:96d6b7b1ab1ec89d6413c09f7fef3fdaac5e9d4293b259492ddda4cf5883d354"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}
	output { 
		File out_file = outfile
	}
}

workflow w {
	Array[File] this_gds_in_arr
	Array[String] this_subset_ids
	Float this_memory
	Int this_disksize
	
	call getScript

	scatter(this_gds_in in this_gds_in_arr) {
		call runScript{ input: gds_in=this_gds_in, subset_ids=this_subset_ids, script=getScript.outscript, memory=this_memory, disksize=this_disksize}
	}

	output {
		Array[File] this_gds_out = runScript.out_file
	}
}