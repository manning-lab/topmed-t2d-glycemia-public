task subset {
	File vcf
	File? samples = "none"
	Boolean filtersamples
	String label
	String out = "${label}.vcf.bgz"

	Float? memory
	Int? disksize
	
	command {
		bcftools view ${true='-S ' false='' filtersamples}${samples} -c 1 ${vcf} | bgzip > ${out}
		tabix -p vcf ${out}

	}
	runtime {
		docker: "vandhanak/bcftools@sha256:2d15d24e39085c3b8a0c716a2dcd7aadd48f33f213661432fe6f279eaf8b293d"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}
	output { 
		File out_file = "${label}.vcf.bgz"
		File out_tbi = "${label}.vcf.bgz.tbi"
	}
}

workflow w {
	Array[File] these_vcf
	File? these_samples
	Boolean this_filtersamples
	String this_label
	
	Float? this_memory = 30.0
	Int? this_disksize = 200
	
	scatter(this_vcf in these_vcf) {
		call subset { 
			input: vcf=this_vcf, samples=these_samples, filtersamples=this_filtersamples, label=this_label, memory=this_memory, disksize=this_disksize}
	}

	output {
		Array[File] this_vcf_out = subset.out_file
		Array[File] this_tbi_out = subset.out_tbi
	}
}