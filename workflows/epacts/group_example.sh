
# From Epacts wiki, example group test (https://genome.sph.umich.edu/wiki/EPACTS#Getting_Started_With_Examples)

# since files are copied simultaneously in DNA nexus, our inputs no longer satisfy the EPACTS requirement of .tbi files being newer than .vcg.gz files. 
# the following line is a work around and should be included whenever you are working with indexed vcf and EPACTS in DNA nexus
# use the format : touch <name_of_vcf_tbi.vcf.gz.tbi>

touch $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz.tbi

# If the VCF is not annotated, 'epacts makegroup' cannot be used. In order to annotate VCF, one can use the example VCF using ANNOVAR as follows:
# The epacts anno script will add "ANNO=[function]:[genename]" entry into the INFO field based on gencodeV7 (default) or refGene database.
# It is important to check whether the VCF file is already annotated or not in order to avoid no or redundant annotation.

$runepacts anno --in $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz --out $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.vcf.gz

# Make the group file

$runepacts make-group --vcf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.vcf.gz --out $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.grp --format epacts --nonsyn

 # To perform a groupwise burden test on the example VCF (annotated as above), run the following command

$runepacts group --vcf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.vcf.gz --groupf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.grp --out test.gene.skat --ped $epactsdata/1000G_dummy_pheno.ped --max-maf 0.05  --pheno QT --cov AGE --cov SEX --test skat --skat-o --run 2

########## expected outputs ##########
# test.gene.skat.phe
# test.gene.skat.ind
# test.gene.skat.cov
# test.gene.skat.Makefile
# test.gene.skat.epacts
# test.gene.skat.epacts.R
# test.gene.skat.epacts.conf
# test.gene.skat.epacts.qq.pdf
# test.gene.skat.epacts.mh.pdf
# test.gene.skat.epacts.top5000
# test.gene.skat.epacts.OK