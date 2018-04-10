task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/make_groups/methods/makeGroups/make_mask3.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/make_groups/methods/makeGroups/combineGroups.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "make_mask3.R"
		File comb_script = "combineGroups.R"
	}
}

task makeGroups {
	File gds
	File expr
	File exon
	File tfbs
	File states
	File ptv
	File genh
	String outpref
	File script

	Int disk
	Int mem 

	command {
        	R --vanilla --args ${gds} ${expr} ${exon} ${tfbs} ${states} ${ptv} ${genh} ${outpref} < ${script}
    }

    meta {
            author: "TM"
            email: "tmajaria@broadinstitute.org"
    }

    runtime {
    	   docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		   disks: "local-disk ${disk} SSD"
		   memory: "${mem}G"
    }

    output {
    	File rdata = "${outpref}.mask3.RData"
    	File csv = "${outpref}.mask3.csv"
    }	

}

task combineGroups {
	Array[File] groups
	String outpref
	File script

	command {
		R --vanilla --args ${sep="," groups} < ${script}
	}

	runtime {
    	   docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		   disks: "local-disk 100 SSD"
		   memory: "10G"
    }

    output {
    	File all_groups = "${outpref}.all.groups.mask3.csv"
    }	
}

workflow wf {
	Array[File] these_gds
	File this_expr
	File this_exon
	File this_tfbs
	File this_states
	File this_ptv
	File this_genh
	String this_outpref

	Int this_disk
	Int this_mem 



	call getScript
	
	scatter(cur_gds in these_gds) {
		call makeGroups {
				input: gds = cur_gds, expr = this_expr, exon = this_exon, tfbs = this_tfbs, states = this_states, ptv = this_ptv, genh = this_genh, outpref = this_outpref, disk = this_disk, mem = this_mem, script = getScript.script
		}
	}

	call combineGroups {
		input: groups = makeGroups.csv, outpref = this_outpref, script = getScript.comb_script
	}
}
