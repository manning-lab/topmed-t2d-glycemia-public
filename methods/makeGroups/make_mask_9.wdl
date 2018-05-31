task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/make_groups/methods/makeGroups/make_mask_9.R"
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/make_groups/methods/makeGroups/combineGroups.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "make_mask_9.R"
		File comb_script = "combineGroups.R"
	}
}

task makeGroups {
	File gds
	File hub
	File reg
	String chr
	File script

	Int disk
	Int mem 

	command {
        	R --vanilla --args ${gds} ${hub} ${reg} ${chr} < ${script}
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
    	File csv = "freeze5b_dp10.${chr}.mask9.csv"
    }	

}

task combineGroups {
	Array[File] groups
	String outpref
	File script

	command {
		R --vanilla --args ${sep="," groups} ${outpref} < ${script}
	}

	runtime {
    	   docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		   disks: "local-disk 100 SSD"
		   memory: "10G"
    }

    output {
    	File all_groups = "${outpref}.all.groups.mask9.csv"
    }	
}

workflow wf {
	Array[File] these_gds
	File this_hub
	File this_reg
	String this_outpref

	Int this_disk
	Int this_mem

	Array[String] chrs = ["chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9"]
	

	Array[Pair[String,File]] chrs_gds =  zip(chrs, these_gds)

	call getScript
	
	scatter(cur_gds in chrs_gds) {
		call makeGroups {
				input: gds = cur_gds.right, hub = this_hub, reg = this_reg, chr = cur_gds.left, disk = this_disk, mem = this_mem, script = getScript.script
		}
	}

	call combineGroups {
		input: groups = makeGroups.csv, outpref = this_outpref, script = getScript.comb_script
	}
}
