task readMap {
	File chrgds
	File chranno
	File chrtfbs

	command {
		ls
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		Map[String,String] ctog = read_map(chrgds)
		Map[String,String] ctoa = read_map(chranno)
		Map[String,String] ctot = read_map(chrtfbs)
	}
}

workflow w {
	File gds
	File anno
	File tfbs
	call readMap {input: chrgds = gds, chranno = anno, chrtfbs = tfbs}
}