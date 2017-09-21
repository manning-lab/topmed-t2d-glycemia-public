task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/methods/fileSubset/subsetGdsByVariant.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File outscript = "subsetGdsByVariant.R"
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

task combine {
	Array[File] in_files 
	String? label = "subset"
	String outfile = label + ".gds"
	Int disksize
	Float memory

	command {
	awk '
    	FNR==1 && NR!=1 { while (/^<header>/) getline; }
    	1 {print}
	' *.gds >${outfile}
	}

	runtime {
		docker: "alpine@sha256:1072e499f3f655a032e88542330cf75b02e7bdf673278f701d7ba61629ee3ebe"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}

	output {
		File out = outfile
	}
}

workflow w {
	Array[File] this_gds_in_arr
	Array[String] this_subset_ids
	Float this_memory
	Int this_disksize
	
	call getScript
	scatter(this_gds in this_gds_in_arr){
		call runScript{ input: gds_in=this_gds, subset_ids=this_subset_ids, script=getScript.outscript, memory=this_memory, disksize=this_disksize}
	}

	call combine {input: runScript.out_file}

	output {
		File this_gds_out = combines.out
	}
}