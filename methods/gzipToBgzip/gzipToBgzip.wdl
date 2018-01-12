task gzipToBgzip {
	File infile
	Int disksize
	Float memory

	String base = basename(infile)
	String outfile = sub(base, ".gz", ".bgz")

	command {
		echo ${base} && echo ${outfile} && gunzip -c ${infile} | bgzip  > ${outfile}
	}

	runtime {
		docker: "biowardrobe2/samtools@sha256:e4dad5f38c1b782d3f1608410c07e8dc47fb7b92bc427175a160dfa0813c48d8"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}

	output {
		File out_file = outfile
	}
}


workflow w {
	Array[File] infiles
	Int thisDisksize
	Float thisMemory

	scatter(this_file in infiles) {
		call gzipToBgzip { input: infile = this_file, disksize = thisDisksize, memory = thisMemory }
	}

	output {
		Array[File] outfiles = gzipToBgzip.out_file
	}
}