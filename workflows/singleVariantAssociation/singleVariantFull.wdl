task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/singleVariantAssociation/association.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/singleVariantAssociation/commonID.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/singleVariantAssociation/summary.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File assoc_script = "association_opt.R"
		File commonID_script = "commonID.R"
		File summary_script = "summary.R"
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

task assocTest {
# perform association test
# assocSingleChr.R

	File gds
	File ped
	File GRM
	File commonIDs
	String label
	String colname
	String outcome
	String outcomeType
	String? covariates
	File assocTestScript

	command {
		R --vanilla --args ${gds} ${ped} ${GRM} ${commonIDs} ${colname} ${label} ${outcome} ${outcomeType} ${covariates} < ${assocTestScript} 
	}

	meta {
		author: "jasen jackson; Alisa Manning, Tim Majarian"
		email: "jasenjackson97@gmail.com; amanning@broadinstitute.org, tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk 100 SSD"
		memory: "30G"
	}

	output {
		File assoc = "${label}.assoc.RData"
	}
}

task summary {
# summarizes results from association test
# saves qq and manhattan plots
# writes tables with association values

	Array[File] assoc
	String pval
	String label
	String title
	File summaryScript

	command {
		R --vanilla --args ${pval} ${label} ${title} ${sep = ' ' assoc} < ${summaryScript}	
	}
	
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
  	    disks: "local-disk 100 SSD"
        memory: "30G"
	}

	output {
		File mhplot = "${label}.mhplot.png"
		File qqplot = "${label}.qqplot.png"
		File topassoccsv = "${label}.topassoc.csv"
		File allassoccsv = "${label}.assoc.csv"
	}
}

workflow w_assocTest {
	File this_ped
	File this_kinshipGDS
	String this_label
	String this_colname
	String this_outcome
	String this_outcomeType
	String? this_covariates
	Array[File] gdsFiles
	String this_pval 
	String this_title

	call getScript
	
	call common_ID {
            input: gds=gdsFiles[0], ped=this_ped, script=getScript.commonID_script, idcol=this_colname, label=this_label
	}
		
	scatter(oneFile in gdsFiles) {
		
		call assocTest {
			input: gds = oneFile, ped = this_ped, GRM = this_kinshipGDS, commonIDs = common_ID.commonIDsRData, colname = this_colname, outcome = this_outcome, outcomeType = this_outcomeType, covariates = this_covariates, assocTestScript = getScript.assoc_script, label=this_label
		}
	}

	call summary {
		input: assoc = assocTest.assoc, pval=this_pval, label=this_label, title=this_title, summaryScript = getScript.summary_script
	}

}
