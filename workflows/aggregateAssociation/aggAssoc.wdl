task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/agg_assoc/workflows/aggregateAssociation/aggregateAssociation.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/agg_assoc/workflows/aggregateAssociation/aggregateSummary.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File assoc_script = "aggregateAssociation.R"
		File summary_script = "aggregateSummary.R"
	}
}

task aggAssocTest {
	File gds_file
	File null_file
	File group_file
	String label
	String? test
	String? pval
	String? weights
	
	
	File script

	Int memory
	Int disk

	command {
		R --vanilla --args ${gds_file} ${null_file} ${group_file} ${label} ${default="SKAT" test} ${default="kuonen" pval} ${default="1,25" weights} < ${script} 
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File assoc = "${label}.assoc.RData"
	}
}

task summary {
	String label
	Array[File] assoc

	File summaryScript

	Int? memory = 10
	Int? disk = 50

	command {
		R --vanilla --args ${label} ${sep = ' ' assoc} < ${summaryScript}	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
  	    disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File mhplot = "${label}.mhplot.pdf"
		File qqplot = "${label}.qqplot.pdf"
		File assoc_res = "${label}.groupAssoc.csv"
	}
}

workflow group_assoc_wf {
	Array[File] these_gds_files
	File this_null_file
	Array[File] these_group_files
	String this_label
	String this_test
	String this_pval
	String this_weights
	Int this_memory
	Int this_disk
	
	Array[Pair[File,File]] these_gds_groups = zip(these_gds_files, these_group_files)

	call getScript

	scatter(this_gds_group in these_gds_groups) {
		
		call aggAssocTest {
			input: gds_file = this_gds_group.left, null_file = this_null_file, groups = this_gds_group.right, label=this_label, test = this_test, pval = this_pval, weights = this_weights, memory = this_memory, disk = this_disk, script = getScript.assoc_script
		}
	}

	call summary {
		input: assoc = aggAssocTest.assoc, label=this_label, memory = this_memory, disk = this_disk, summaryScript = getScript.summary_script
	}

}