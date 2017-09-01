# adding a task that takes the vcf file and creates a bgzipped vcf and a tabix index file of the vcf.
# can be optional output

task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/methods/bcfToGds/vcfToGds.R"
	}

	runtime {
		docker: "tmajarian/wget@sha256:2808366a79e76c01fdfa007413f400c811918060f0c25d1c7b28d5e19c504215"
	}

	output {
		File outscript = "vcfToGds.R"
	}
}

task runScript {
	File in
	File indat
	Int disksize
	Float memory
	
	command {
		R --vanilla --args ${indat} < ${in}
	}
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:96d6b7b1ab1ec89d6413c09f7fef3fdaac5e9d4293b259492ddda4cf5883d354"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}
	output { 
		File out = read_string("output.txt")
	}
}

workflow w {
	File infile

	call getScript
	call runScript {input: in=getScript.outscript, indat=infile}

	output {
		File outfile = runScript.out
	}
}