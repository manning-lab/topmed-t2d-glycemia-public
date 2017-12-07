task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/singleVariantAssociation/association_null_in.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/mdev/workflows/singleVariantAssociation/summary_v2.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File assoc_script = "association_null_in.R"
		File summary_script = "summary_v2.R"
	}
}

task assocTest {
	File gds
	File null_file
	String label
	String test
	File script

	command {
		R --vanilla --args ${gds} ${null_file} ${label} ${test} < ${script} 
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk 100 SSD"
		memory: "40G"
	}

	output {
		File assoc = "${label}.assoc.RData"
	}
}

task summary {
	Array[File] assoc
	String pval
	String label
	File script

	command {
		R --vanilla --args ${pval} ${label} ${sep = ' ' assoc} < ${script}	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
  	    disks: "local-disk 100 SSD"
        memory: "30G"
	}

	output {
		File allqqplot = "${label}.all.qqplot.png"
		File commonqqplot = "${label}.common.qqplot.png"
		File uncommonqqplot = "${label}.uncommon.qqplot.png"

		File allmhplot = "${label}.all.mhplot.png"
		File commonmhplot = "${label}.common.mhplot.png"
		File uncommonmhplot = "${label}.uncommon.mhplot.png"

		File topassoccsv = "${label}.topassoc.csv"
		File allassoccsv = "${label}.assoc.csv"
	}
}

workflow w_assocTest {
	Array[File] gdsFiles
	File this_null
	String this_label
	String this_test
	String this_pval

	call getScript
	
	scatter(oneFile in gdsFiles) {
		
		call assocTest {
			input: gds = oneFile, null_file = this_null, label = this_label, test = this_test, script = getScript.assoc_script
		}
	}

	call summary {
		input: assoc = assocTest.assoc, pval=this_pval, label=this_label, script = getScript.summary_script
	}

}