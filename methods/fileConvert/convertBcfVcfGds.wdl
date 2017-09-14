task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/methods/fileConvert/convertBcfVcfGds.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File outscript = "convertBcfVcfGds.R"
	}
}

task runScript {
	File script
	File indat
	Boolean? gds_flag = false
	Boolean? vcf_flag = false
	Int disksize
	Float memory
	
	command {
		R --vanilla --args ${indat} ${gds_flag} ${vcf_flag} < ${script}
	}
	runtime {
		docker: "akmanning/r-mkl-bioconductor-bcftools@sha256:432300b876907590fc0778dd214f35718477a670acc0501427dcb5a30486ae4a"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}
	output { 
		File? gds_out = read_string("output_gds.txt")
		File? vcf_out = read_string("output_vcf.txt")
		File? tbi_out = read_string("output_tbi.txt")
	}
}

workflow w {
	Array[File] infiles
	Boolean this_gds_flag
	Boolean this_vcf_flag
	Int thisDisksize
	Float thisMemory

	call getScript
	scatter(this_file in infiles) {
		call runScript {input: script=getScript.outscript, indat=this_file, gds_flag=this_gds_flag, vcf_flag=this_vcf_flag, disksize=thisDisksize, memory=thisMemory}
	}

	output {
		Array[File]? gdsOut = select_all(runScript.gds_out)
		Array[File]? vcfzip = select_all(runScript.vcf_out)
		Array[File]? vcfind = select_all(runScript.tbi_out)
	}
}