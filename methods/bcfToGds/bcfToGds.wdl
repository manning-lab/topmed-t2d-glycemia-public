task convertBCFtoVCF {
	File bcf
	Int disksize
	Int memory
	String base = basename(bcf)
	String vcf_out_file = sub(base, "\\.bcf$", ".vcf")

	
	command {
		bcftools convert --output-type z ${bcf} -o ${vcf_out_file}
	}
	
	meta {
		author: "Alisa Manning, T Majarian"
		email: "amanning@broadinstitute.org, tmajaria@broadinstitute.org"
	}

	runtime {
		# docker: "vandhanak/bcftools:1.3.1"
		docker: "vandhanak/bcftools@sha256:2d15d24e39085c3b8a0c716a2dcd7aadd48f33f213661432fe6f279eaf8b293d"
		disks: "local-disk ${disksize} SSD"
        memory: "${memory}G"
	}

	output {
		File vcfOut = vcf_out_file
	}
}

task convertVCFtoGDS {
	File vcf_in
	String method
	Int disksize
	Int memory
	String base = basename(vcf_in)
	String gds_out_file = sub(base, "\\.vcf$", ".gds")

	command {
		R --vanilla --args ${vcf_in} ${method} < 
	}

	meta {
		author: "jasen jackson, T Majarian"
		email: "jasenjackson97@gmail.com, tmajaria@broadinstitute.org"
	}
	
	runtime {
		# docker: "robbyjo/r-mkl-bioconductor:3.4.1"
		docker: "robbyjo/r-mkl-bioconductor@sha256:96d6b7b1ab1ec89d6413c09f7fef3fdaac5e9d4293b259492ddda4cf5883d354"
		disks: "local-disk ${disksize} SSD"
        memory: "${memory}G"
	}

	output {
		File gdsOut = gds_out_file
	}
}

workflow bcfToGds {
	File bcf_file
	File this_script
	String this_method
	
	call convertBCFtoVCF {
		input: bcf=bcf_file
	}

	call convertVCFtoGDS {
		input: vcf_in=convertBCFtoVCF.vcfOut, method=this_method
	}

	output {
		File out = convertVCFtoGDS.gdsOut
	}
}
