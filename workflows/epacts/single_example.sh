
# From Epacts wiki, example single variant test (https://genome.sph.umich.edu/wiki/EPACTS#Getting_Started_With_Examples)
# The command above will perform single variant association test using a dummy case-control phenotype file and a subset of 1000 genomes exome VCF file (chr20) using score test statistic for all variants over 1% of higher MAF using 2 parallel runs.
# You will see the 4 output files as the main outcome of the analysis

# since files are copied simultaneously in DNA nexus, our inputs no longer satisfy the EPACTS requirement of .tbi files being newer than .vcg.gz files. 
# the following line is a work around and should be included whenever you are working with indexed vcf and EPACTS in DNA nexus
# use the format : touch <name_of_vcf_tbi.vcf.gz.tbi>

touch $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz.tbi

# now we can run the epacts command

$runepacts single --vcf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz --ped  $epactsdata/1000G_dummy_pheno.ped  --min-maf 0.001 --chr 20 --pheno DISEASE --cov AGE --cov SEX --test b.score --anno --out test --run 2

########## expected outputs ##########
# /test.epacts.top5000
# /test.phe
# /test.epacts.R
# /test.epacts.conf
# /test.epacts.OK
# /test.ind
# /test.epacts.qq.pdf
# /test.cov
# /test.epacts.gz
# /test.Makefile
# /test.epacts.mh.pdf
# /test.epacts.gz.tbi

touch $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz.tbi