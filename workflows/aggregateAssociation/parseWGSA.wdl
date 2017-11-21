task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/workflows/aggregateAssociation/wgsa_parser.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File script = "wgsa_parser.R"
	}
}

task parse {
	File anno
	String label
	String desired_columns
	String columns_to_split
	File script

	Int? disk = 100
	Int? mem = 10

	command {
        	R --vanilla --args ${anno} ${label} ${desired_columns} ${columns_to_split} < ${script}
    }

    meta {
            author: "TM"
            email: "tmajaria@broadinstitute.org"
    }

    runtime {
    	   docker: "tmajarian/wgsa_parser"
		   disks: "local-disk ${disk} SSD"
		   memory: "${mem}G"
    }

    output {
    	File anno_out = "${label}.tsv"
    }	

}

workflow wf {
	Map[String,File] these_chr_anno
	String this_label
	String these_cols
	String these_split

	Int? this_disk
	Int? this_mem

	Array[String] chrs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X]

	call getScript
	
	scatter(this_chr in chrs) {
		call parse {
				input: anno=these_chr_anno[this_chr], label=this_label, desired_columns=these_cols, columns_to_split=these_split, script=getScript.script, disk=this_disk, mem=this_mem
		}
	}
}