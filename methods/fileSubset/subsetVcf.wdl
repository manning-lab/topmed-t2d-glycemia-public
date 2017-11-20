task subset {
	File vcf
	File? samples = "none"
	Boolean filtersamples
	File? variants = "none"
	Boolean filtervariants
	String label
	String out = "${label}.vcf.bgz"

	File script
	Float? memory
	Int? disksize
	
	command {
		vcftools --gzvcf ${vcf} ${true='--snps ' false='' filtervariants}${variants} ${true='--indv ' false='' filtersamples}${samples} --recode --recode-INFO-all --stdout | bgzip > ${out}
		tabix -p vcf ${out}
	}
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:96d6b7b1ab1ec89d6413c09f7fef3fdaac5e9d4293b259492ddda4cf5883d354"
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
	File these_samples
	File these_variants
	String this_label
	File script
	
	Float? this_memory = 10.0
	Int? this_disksize = 20
	
	scatter(this_vcf in these_vcf) {
		call subset { 
			input: vcf=this_vcf, samples=these_samples, variants=these_variants, label=this_label, script=script, memory=this_memory, disksize=this_disksize}
	}

	output {
		Array[File] this_vcf_out = subset.out_file
		Array[File] this_tbi_out = subset.out_tbi
	}
}