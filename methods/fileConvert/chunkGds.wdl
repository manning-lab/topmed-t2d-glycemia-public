# adding a task that takes the vcf file and creates a bgzipped vcf and a tabix index file of the vcf.
# can be optional output

task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/methods/fileConvert/chunkGds.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File outscript = "chunkGds.R"
	}
}

task chunk {
	File script
	File gds_in
	Int disksize
	Float memory
	
	command {
		R --vanilla --args ${gds_in} < ${script}
	}

	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:96d6b7b1ab1ec89d6413c09f7fef3fdaac5e9d4293b259492ddda4cf5883d354"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}

	output { 
		Array[File] out = glob("*.gds")
	}

}

task flattenArray {
    Array[Array[String]] arr

    command{
    	echo "${sep='\n' arr}" > raw_array
    	sed -e 's/\n//g' -e 's/\[//g' -e $'s/\]//g' -e $'s/ /\\\n/g' -e 's/,//g' raw_array > file_of_filenames
    }

    runtime {
    	docker: "frolvlad/alpine-bash@sha256:1c9aa5fb2c5feaacff1ea25de537116235726430c18f9061585c0d1283bcf0dd"
    }

    output{
        Array[String] filenames = read_lines("file_of_filenames")
    }
}

workflow w {
	Array[File] infiles
	# File infile
	Int thisDisksize
	Float thisMemory

	call getScript

	scatter(this_file in infiles) {
		call chunk {input: script=getScript.outscript, gds_in=this_file, disksize=thisDisksize, memory=thisMemory}
	}

	call flattenArray { input: arr=chunk.out }

	output { Array[File] chunks = flattenArray.filenames }

}
