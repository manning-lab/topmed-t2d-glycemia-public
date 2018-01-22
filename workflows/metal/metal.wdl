task makeScript {
	Array[File] result_files
	String marker_column
	String weight_column
	String allele_effect_column
	String allele_non_effect_column
	String freq_column
	String pval_column
	String out_pref
	String out_file

	command {
		echo "# this is your metal script" > ${out_file}
		echo "MARKER ${marker_column}" >> ${out_file}
		echo "WEIGHT ${weight_column}" >> ${out_file}
		echo "ALLELE ${allele_effect_column} ${allele_non_effect_column}" >> ${out_file}
		echo "FREQ ${freq_column}" >> ${out_file}
		echo "PVAL ${pval_column}" >> ${out_file}
		echo "PROCESS ${sep = "\nPROCESS " result_files}" >> ${out_file}
		echo "OUTFILE ${out_pref} .TBL" >> ${out_file}
		echo "SEPARATOR ${default= "COMMA" separator}" >> ${out_file}
		echo "ANALYZE" >> ${out_file}
	}

	runtime {
		docker: "ubuntu@sha256:d3fdf5b1f8e8a155c17d5786280af1f5a04c10e95145a515279cf17abdf0191f"
	}

	output {
		File out = out_file
	}
}

task runMetal {
	Array[File] result_files
	String marker_column
	String weight_column
	String allele_effect_column
	String allele_non_effect_column
	String freq_column
	String pval_column
	String out_pref
	String? separator
	String out_file

	Int memory
	Int disk

	command {
		echo "# this is your metal out_file" > ${out_file}
		echo "MARKER ${marker_column}" >> ${out_file}
		echo "WEIGHT ${weight_column}" >> ${out_file}
		echo "ALLELE ${allele_effect_column} ${allele_non_effect_column}" >> ${out_file}
		echo "FREQ ${freq_column}" >> ${out_file}
		echo "PVAL ${pval_column}" >> ${out_file}
		echo "PROCESS ${sep = "\nPROCESS " result_files}" >> ${out_file}
		echo "OUTFILE ${out_pref} .TBL" >> ${out_file}
		echo "SEPARATOR ${default= "COMMA" separator}" >> ${out_file}
		echo "ANALYZE" >> ${out_file}
		metal ${out_file} > "${out_pref}.log"
	}

	runtime {
		docker: "tmajarian/metal@sha256:27e2c3189cff9c974d814a3ea3968330160d068b98e7d53dec85c098e5c57c1a"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		Array[File] out_files = glob("${out_pref}*")
	}
}

workflow w_metal {
	Array[File] these_result_files
	String this_marker_column
	String this_weight_column
	String this_allele_effect_column
	String this_allele_non_effect_column
	String this_freq_column
	String this_pval_column
	String this_out_pref
	String? this_separator
	String this_out_file

	Int this_memory
	Int this_disk

	

	call runMetal {
		input: result_files = these_result_files, marker_column = this_marker_column, weight_column = this_weight_column, allele_effect_column = this_allele_effect_column, allele_non_effect_column = this_allele_non_effect_column, freq_column = this_freq_column, pval_column = this_pval_column, out_pref = this_out_pref, separator = this_separator, out_file = this_out_file, memory = this_memory, disk = this_disk
		
	}
}