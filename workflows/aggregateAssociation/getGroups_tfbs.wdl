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
	String annotation_values
	String anno_field
	File tfbs
	File state
	String chromatin_states
	String label
	Float min_maf
	Int pad
	String tfbs_states
	File groupScript

	Int? memory = 10
	Int? disk = 50

	command {
		R --vanilla --args ${gds} ${genes} ${anno} ${annotation_values} ${anno_field} ${tfbs} ${state} ${chromatin_states} ${label} ${min_maf} ${pad} ${tfbs_states} < ${groupScript}
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
	File this_genes
	Map[Int,File] chr_anno
	String this_anno_val
	String this_anno_field
	Map[Int,File] chr_tfbs
	File this_state
	String this_chr_states
	String this_label
	Float this_maf
	Int this_pad
	String this_tfbs_states

	Array[Int] chrs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

	Int? this_memory
	Int? this_disk 
	

	call getScript
	
	scatter(chr in chrs) {
		call get_groups {
				input: gds=chr_gds[chr], genes=this_genes, anno=chr_anno[chr], annotation_values=this_anno_val, anno_field=this_anno_field, tfbs=chr_tfbs[chr], state=this_state, chromatin_states=this_chr_states, label=this_label, min_maf=this_maf, pad=this_pad, tfbs_states=this_tfbs_states, memory=this_memory, disk=this_disk, groupScript=getScript.group_script
		}
	}
}