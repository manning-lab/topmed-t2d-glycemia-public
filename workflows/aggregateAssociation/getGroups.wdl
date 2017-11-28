task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/aggregateAssociation/groupByGene.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File group_script = "groupByGene.R"
	}
}

task get_groups {
	File gds
	File genes
	File anno
	String annotation_values
	String anno_field
	File state
	String chromatin_states
	String label
	Float min_maf
	Int pad
	File groupScript

	Int? memory = 10
	Int? disk = 50

	command {
		R --vanilla --args ${gds} ${genes} ${anno} ${annotation_values} ${anno_field} ${state} ${chromatin_states} ${label} ${min_maf} ${pad} < ${groupScript}
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
	Array[File] these_gds
	File this_genes
	Array[File] these_anno
	String these_annotation_values
	String this_anno_field
	File this_state
	String this_chromatin_states
	String this_label

	Float this_min_maf
	Int this_pad

	Int? this_memory
	Int? this_disk 
	
	Array[Pair[File,File]] these_gds_anno = zip(these_gds,these_anno)

	call getScript
	
	scatter(this_gds_anno in these_gds_anno) {
		call get_groups {
				input: gds=this_gds_anno.left, genes=this_genes, anno=this_gds_anno.right, annotation_values=these_annotation_values, anno_field=this_anno_field, state=this_state, chromatin_states=this_chromatin_states, label=this_label, min_maf=this_min_maf, pad=this_pad, memory=this_memory, disk=this_disk, groupScript=getScript.group_script
		}
	}
}