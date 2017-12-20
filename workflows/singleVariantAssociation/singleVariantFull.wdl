task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/singleVariantAssociation/association.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/singleVariantAssociation/summary.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File assoc_script = "association.R"
		File summary_script = "summary.R"
	}
}

task assocTest {

	File gds
	File ped
	File GRM
	File commonIDs
	String colname
	String label
	String outcome
	String outcomeType
	Int minmac
	String? covariates
	File assocTestScript

	Int memory

	command {
		R --vanilla --args ${gds} ${ped} ${GRM} ${commonIDs} ${colname} ${label} ${outcome} ${outcomeType} ${minmac} ${covariates} < ${assocTestScript} 
	}

	meta {
		author: "jasen jackson; Alisa Manning, Tim Majarian"
		email: "jasenjackson97@gmail.com; amanning@broadinstitute.org, tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk 100 SSD"
		memory: "${memory}G"
	}

	output {
		File assoc = "${label}.assoc.RData"
	}
}

task summary {
	Array[File] assoc
	String test
	String label
	File script

	Int memory

	command {
		R --vanilla --args ${test} ${label} ${sep = ' ' assoc} < ${script}	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
  	    disks: "local-disk 100 SSD"
        memory: "${memory}G"
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
	File this_ped
	File this_kinshipGDS
	File this_sampleids
	String this_colname
	String this_label
	String this_outcome
	String this_outcomeType
	Int this_minmac
	String? this_covariates
	String this_test

	call getScript
		
	scatter(oneFile in gdsFiles) {
		
		call assocTest {
			input: gds = oneFile, ped = this_ped, GRM = this_kinshipGDS, commonIDs = this_sampleids, colname = this_colname, label=this_label, outcome = this_outcome, outcomeType = this_outcomeType, minmac = this_minmac, covariates = this_covariates, assocTestScript = getScript.assoc_script
		}
	}

	call summary {
		input: assoc = assocTest.assoc, test=this_test, label=this_label, script = getScript.summary_script
	}

}