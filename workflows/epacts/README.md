# Running EPACTS in DNAnexus with the *Swiss Army Knife* app

Here is a quick tuturial on running EPACTS in DNAnexus. We use custom shell scripts, called from the *Swiss Army Knife* command line to run EPACTS commands (just like you would locally) together with a docker image containing all of the EPACTS code. This is entirely based on the examples provided with the [EPACTS documentation](https://genome.sph.umich.edu/wiki/EPACTS#Getting_Started_With_Examples).

## Getting Started

First, you'll need to login to DNAnexus and create a new folder for your EPACTS scripts. 

### Prerequisites

Docker image available on docker hub: [tmajarian/epacts_rmkl:1.1](https://hub.docker.com/r/tmajarian/epacts_rmkl/)


Tutorial scripts available on [GitHub](https://github.com/manning-lab/topmed-t2d-glycemia-public/tree/dev/workflows/epacts)

You do not need to download the docker image to your local machine, although it can be useful for reproducing the compute environment in DNAnexus while testing. Download the tutorial scripts and upload to your EPACTS scripts folder.

## Writing scripts for DNAnexus

To run EPACTS, we need to write shell scripts containing calls to EPACTS that we would normally use from the command line. (Note: within the docker image, you can use $runepacts for your epacts command.)

```bash
# test_run_epacts.sh
/EPACTS_install/EPACTS-3.2.6/bin/test_run_epacts.sh
```

## Running the tutorial scripts

### Tutorial 1 -- Testing

1. From your EPACTS scripts folder, create a new workflow and rename to something informative like _Epacts test run_.
2. Choose _Swiss Army Knife_ under the _Add a Step_ menu. 
3. Add the first script *test_run.sh* to the workflow input files. This will call the *test_run_epacts.sh* script found on the [EPACTS wiki](https://genome.sph.umich.edu/wiki/EPACTS#Getting_Started_With_Examples).
4. Click _Swiss Army Knife_ to configure the workflow.
   1. Set the instance type to something very minimal; these scripts do not take much power.
   2. For the command line call, add *sh test_run_epacts.sh*
   3. For the docker image, add *tmajarian/epacts_rmkl:1.1*
5. Run workflow and wait to complete.
6. Results files will default to workflow home directory.

### Tutorial 2 --  Single Variant analysis

For single variant analysis using zipped vcf (with index) files, we need to add an extra step when importing the input files. EPACTS requires index files to be newer than the vcf. Since DNAnexus copies files simultaneously, this will not be true when the command is run. To get around this, we can add a **touch** command to the script being called.

```bash
# single_example.sh
touch $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz.tbi
$runepacts single --vcf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz --ped  $epactsdata/1000G_dummy_pheno.ped  --min-maf 0.001 --chr 20 --pheno DISEASE --cov AGE --cov SEX --test b.score --anno --out test --run 2
```

1. Make a new workflow, choosing the Swiss Army knife again.
2. Configure the tool as before, except for the command line call.
3. Add *sh single_example.sh* for the command line.
4. Run the workflow and wait for out (located in ./singleVar/)

### Tutorial 3 -- Groupwise analysis

For group testing (burden/skat), we have a few more steps in the script. We annotate the vcf file, create a group file, and run the test (in this case skat-o).

```bash
# group_example.sh
touch $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz.tbi

$runepacts anno --in $epactsdata/1000G_exome_chr20_example_softFiltered.calls.vcf.gz --out $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.vcf.gz

$runepacts make-group --vcf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.vcf.gz --out $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.grp --format epacts --nonsyn

$runepacts group --vcf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.vcf.gz --groupf $epactsdata/1000G_exome_chr20_example_softFiltered.calls.anno.grp --out test.gene.skat --ped $epactsdata/1000G_dummy_pheno.ped --max-maf 0.05  --pheno QT --cov AGE --cov SEX --test skat --skat-o --run 2
```

1. Make a new workflow and configure as before.
2. For the command line, use *sh group_example.sh*
