task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/make_groups/methods/makeGroups/make_mask3.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "make_mask3.R"
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
    	File rdata = "${outpref}.${chr}.mask3.RData"
    	File csv = "${outpref}.${chr}.mask3.csv"
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
				input: gds = cur_gds, expr = this_expr, tfbs = this_tfbs, states = this_states, ptv = this_ptv, genh = this_genh, outpref = this_outpref, disk = this_disk, mem = this_mem, script = getScript.script
		}
	}
}
