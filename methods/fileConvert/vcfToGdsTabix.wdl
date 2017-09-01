# adding a task that takes the vcf file and creates a bgzipped vcf and a tabix index file of the vcf.
# can be optional output

task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/methods/fileConvert/vcfToGds.R"
	}

	runtime {
		docker: "tmajarian/wget@sha256:2808366a79e76c01fdfa007413f400c811918060f0c25d1c7b28d5e19c504215"
	}

	output {
		File outscript = "vcfToGds.R"
	}
}

task vcfToGds {
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

task vcfToBgvcf {
	File vcf_in
	Int diskSize
	Float Memory
	File bgfile = vcf_in + ".gz"

	command {
		bgzip -c ${vcf_in} > ${bgfile}
		tabix -p vcf ${bgfile}
	}

	runtime {
		docker: "biowardrobe2/samtools@sha256:e4dad5f38c1b782d3f1608410c07e8dc47fb7b92bc427175a160dfa0813c48d8"
		disks: "local-disk ${diskSize} SSD"
		memory: "${Memory}G"
	}

	output {
		File zip = "${vcf_in}.gz"
		File ind = "${vcf_in}.gz.tbi"
	}
}

workflow w {
	File infile
	Int thisDisksize
	Float thisMemory
	Boolean gzip
	Boolean gds

	if(gds) {
		call getScript
		call vcfToGds {input: in=getScript.outscript, indat=infile, disksize=thisDisksize, memory=thisMemory}
	}

	if(gzip) {
		call vcfToBgvcf {input: vcf_in=infile, diskSize=thisDisksize, Memory=thisMemory}
	}
}