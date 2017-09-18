workflow regionalsummary {
	call taskregionalsummary
}

task taskregionalsummary {
        String study
        File regionalscript
        File results
	String label
	File inputfilescript
	File functionsfile
	File gtexadiposesubfile
	File gtexadiposevisfile
	File gtexpancreasfile
	File gtexmusclefile
	File gtexliverfile
	String pvalthreshold
	String mafthreshold
	

command {
ls -la && Rscript -e 'rmarkdown::render("${regionalscript}",output_file="${label}_regional.html",output_dir="/cromwell_root/")' ${results} ${study} ${label} ${inputfilescript} ${functionsfile} ${gtexadiposesubfile} ${gtexadiposevisfile} ${gtexpancreasfile} ${gtexmusclefile} ${gtexliverfile} ${pvalthreshold} ${mafthreshold} && ls -la && ls -la /cromwell_root/ 
}

        output {
                File regional_html = "${label}_regional.nb.html"
        }
        runtime {
        	docker: "bioconductor/release_core2"
        	memory: "40GB"
    	 	disks: "local-disk 100 SSD"

    }


}


