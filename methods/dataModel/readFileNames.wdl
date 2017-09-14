task read {
	File list

	command {
		ls
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		Array[File] filepaths = read_lines(list)
	}
}

workflow w {
	File list_file
	call read {input: list=list_file}
}