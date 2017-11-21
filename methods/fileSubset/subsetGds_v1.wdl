task subset {
	File gds
	File? samples = "none"
	File? variants = "none"
	String label

	File script
	Float? memory
	Int? disksize
	
	command {
		R --vanilla --args ${gds} ${samples} ${variants} ${label} < ${script}
	}
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:96d6b7b1ab1ec89d6413c09f7fef3fdaac5e9d4293b259492ddda4cf5883d354"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}
	output { 
		File out_file = "${label}.gds"
	}
}

workflow w {
	Array[File] these_gds
	File these_samples
	File these_variants
	String this_label
	File script
	
	Float? this_memory = 10.0
	Int? this_disksize = 20
	
	scatter(this_gds in these_gds) {
		call subset { 
			input: gds=this_gds, samples=these_samples, variants=these_variants, label=this_label, script=script, memory=this_memory, disksize=this_disksize}
	}

	output {
		Array[File] this_gds_out = subset.out_file
	}
}