task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/aggregateAssociation/groupByGene_tfbs.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File group_script = "groupByGene_tfbs.R"
	}
}

task get_groups {
	File gds
	File genes
	File anno
	File tfbs
	File state
	String label

	File groupScript

	Int? memory = 10
	Int? disk = 50

	command {
		R --vanilla --args ${gds} ${genes} ${anno} ${tfbs} ${state} ${label} < ${groupScript}
    }

    meta {
		author: "TM"
		email: "tmajaria@broadinstitute.org"
    }

    runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk ${disk} SSD"
	    memory: "${memory}G"
    }

    output {
    	Pair[File,File] out_groups = (gds,"${label}_groups.RData")
    }	

}

workflow wf {
	Map[Int,File] chr_gds
	Map[Int,File] chr_anno
	Map[Int,File] chr_tfbs
	File this_genes
	File this_state
	String this_label

	Array[Int] chrs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

	Int? this_memory
	Int? this_disk 
	

	call getScript
	
	scatter(chr in chrs) {
		call get_groups {
				input: gds=chr_gds[chr], genes=this_genes, anno=chr_anno[chr], tfbs=chr_tfbs[chr], state=this_state, label=this_label, memory=this_memory, disk=this_disk, groupScript=getScript.group_script
		}
	}
}