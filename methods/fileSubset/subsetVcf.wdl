task subset {
	File vcf
	File? samples = "none"
	Boolean filtersamples
	File? variants = "none"
	Boolean filtervariants
	String label
	String out = "${label}.vcf.bgz"

	Float? memory
	Int? disksize
	
	command {
		vcftools --gzvcf ${vcf} ${true='--snps ' false='' filtervariants}${variants} ${true='--indv ' false='' filtersamples}${samples} --recode --recode-INFO-all --stdout | bgzip > ${out}
		tabix -p vcf ${out}
	}
	runtime {
		docker: "biocontainers/vcftools@sha256:2e9dbb86adc69833634e35512a0791bc7ce25e87e1a0f4604a6df15d46ed4f7f"
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
	Boolean this_filtersamples
	File these_variants
	Boolean this_filtervariants
	String this_label
	
	Float? this_memory = 10.0
	Int? this_disksize = 20
	
	scatter(this_vcf in these_vcf) {
		call subset { 
			input: vcf=this_vcf, samples=these_samples, filtersamples=this_filtersamples, variants=these_variants, filtervariants=this_filtervariants, label=this_label, memory=this_memory, disksize=this_disksize}
	}

	output {
		Array[File] this_vcf_out = subset.out_file
		Array[File] this_tbi_out = subset.out_tbi
	}
}