task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/methods/fileConvert/convertBcfVcfGds.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File outscript = "convertBcfVcfGds.R"
	}
}

task runGds {
	File script
	File indat
	Int disksize
	Float memory

	String no_path = basename(indat)
	String no_vcf = sub(no_path,".vcf","")
	String no_ext = sub(no_vcf,".bcf","")
	String outfile = no_ext + ".gds"

	String out_type = "gds"

	command {
		R --vanilla --args ${indat} ${out_type} ${outfile} < ${script}
	}

	runtime {
		docker: "akmanning/r-mkl-bioconductor-bcftools@sha256:432300b876907590fc0778dd214f35718477a670acc0501427dcb5a30486ae4a"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}

	output {
		File out_file = outfile
	}
}

task runVcf {
	File script
	File indat
	Int disksize
	Float memory

	String no_path = basename(indat)
	String no_vcf = sub(no_path,".vcf","")
	String no_ext = sub(no_vcf,".bcf","")
	String outfile = no_ext + ".vcf.gz"
	String outfile_tbi = no_ext + ".vcf.gz.tbi"

	String out_type = "vcf"

	command {
		R --vanilla --args ${indat} ${out_type} ${outfile} ${outfile_tbi} < ${script}
	}

	runtime {
		docker: "akmanning/r-mkl-bioconductor-bcftools@sha256:432300b876907590fc0778dd214f35718477a670acc0501427dcb5a30486ae4a"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}

	output {
		File out_file = outfile
		File out_file_tbi = outfile_tbi
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
		if (this_gds_flag) {
			call runGds { input: script=getScript.outscript, indat=this_file, disksize=thisDisksize, memory=thisMemory}
		}

		if (this_vcf_flag) {
			call runVcf { input: script=getScript.outscript, indat=this_file, disksize=thisDisksize, memory=thisMemory}
		}
	}

	output {
		Array[File?]? gds_out = select_all(runGds.out_file)
		Array[File?]? vcf_out = select_all(runVcf.out_file)
		Array[File?]? tbi_out = select_all(runVcf.out_file_tbi)
	}
}