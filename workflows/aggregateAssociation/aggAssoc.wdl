task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/aggregateAssociation/aggregateAssociationGroup.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/singleVariantAssociation/commonID.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/aggregateAssociation/aggregateSummary.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File assoc_script = "aggregateAssociationGroup.R"
		File commonID_script = "commonID.R"
		File summary_script = "aggregateSummary.R"
	}
}

task common_ID {
# get only those sample IDs that are in both gds and ped files
# commonID.R

        File gds
        File ped
        File script
        String idcol
        String label

        command {
                R --vanilla --args ${gds} ${ped} ${idcol} ${label} < ${script}
        }

        meta {
                author: "jasen jackson"
                email: "jasenjackson97@gmail.com"
        }

        runtime {
    	   docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		   disks: "local-disk 100 SSD"
		   memory: "3G"
        }

        output {
                File commonIDstxt = "${label}.commonIDs.txt"
                File commonIDsRData = "${label}.commonIDs.RData"
        }
}

task aggAssocTest {
# perform association test

	File gds
	File ped
	File GRM
	File commonIDs
	File groups
	String colname
	String label
	String outcome
	String outcomeType
	String test
	String pval
	String? covariates
	File assocTestScript

	command {
		R --vanilla --args ${gds} ${ped} ${GRM} ${commonIDs} ${groups} ${colname} ${label} ${outcome} ${outcomeType} ${covariates} < ${assocTestScript} 
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk 50 SSD"
		memory: "20G"
	}

	output {
		File assoc = "${label}.assoc.RData"
	}
}

task summary {
	File assoc
	String label
	File summaryScript

	command {
		R --vanilla --args ${assoc} ${label} < ${summaryScript}	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
  	    disks: "local-disk 100 SSD"
        memory: "30G"
	}

	output {
		File mhplot = "${label}.mhplot.pdf"
		File qqplot = "${label}.qqplot.pdf"
		File assoc_res = "${label}.groupAssoc.csv"
	}
}

workflow w_assocTest {
	File this_gds
	File this_ped
	File this_GRM
	File this_groups
	String this_colname
	String this_label
	String this_outcome
	String this_outcomeType
	String this_test
	String this_pval
	String? this_covariates
	
	call getScript
	
	call common_ID {
            input: gds=this_gds, ped=this_ped, script=getScript.commonID_script, idcol=this_colname, label=this_label
	}
		
	call aggAssocTest {
		input: gds = this_gds, ped = this_ped, GRM = this_GRM, commonIDs = common_ID.commonIDsRData, groups = this_groups, colname = this_colname, label=this_label, outcome = this_outcome, outcomeType = this_outcomeType, test = this_test, pval = this_pval, covariates = this_covariates, assocTestScript = getScript.assoc_script
	}


	call summary {
		input: assoc = aggAssocTest.assoc, label=this_label, summaryScript = getScript.summary_script
	}

}