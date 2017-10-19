task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/aggregateAssociation/parse_wgsa.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "parse_wgsa.R"
	}
}

task parse {
	File script
	File anno
	String anno_unzip = basename(anno,".gz")

	command {
        	R --vanilla --args ${anno} ${anno_unzip} < ${script}
    }

    meta {
            author: "TM"
            email: "tmajaria@broadinstitute.org"
    }

    runtime {
    	   docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		   disks: "local-disk 100 SSD"
		   memory: "30G"
    }

    output {
    	File anno_out = "${anno_unzip}.csv"
    }	

}

workflow wf {
	Array[File] these_anno

	call getScript
	
	scatter(this_anno in these_anno) {
		call parse {
				input: script=getScript.script, anno=this_anno
		}
	}
}