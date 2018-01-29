task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/aggregateAssociation/aggregateAssociation.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/aggregateAssociation/aggregateSummary.R"
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
	File gds
	String label
	String test
	String pval
	File groups
	File model_file
	
	File assocTestScript

	Int? memory = 10
	Int? disk = 50

	command {
		R --vanilla --args ${gds} ${label} ${test} ${pval} ${groups} ${model_file} < ${assocTestScript} 
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
	Array[Pair[File,File]] these_gds_groups
	File this_model
	String this_label
	String this_test
	String this_pval
	Int? this_memory
	Int? this_disk
	
	call getScript

	scatter(this_gds in these_gds_groups) {
		
		call aggAssocTest {
			input: gds = this_gds.left, groups = this_gds.right, model_file = this_model, label=this_label, test = this_test, pval = this_pval, memory = this_memory, disk = this_disk, assocTestScript = getScript.assoc_script
		}
	}

	call summary {
		input: assoc = aggAssocTest.assoc, label=this_label, memory = this_memory, disk=this_disk, summaryScript = getScript.summary_script
	}

}
