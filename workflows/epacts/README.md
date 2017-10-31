# Running EPACTS in DNAnexus

Here is a quick tuturial on running EPACTS in DNAnexus. We use custom shell scripts, called from the *Swiss Army Knife* command line to run EPACTS commands (just like you would locally) together with a docker image containing all of the EPACTS code. This is entirely based on the examples provided with the [EPACTS documentation](https://genome.sph.umich.edu/wiki/EPACTS#Getting_Started_With_Examples).

We can also use the *cloud-workstation app* to run EPACTS from the command line within DNAnexus.

## Getting Started -- Swiss Army Knife

First, you'll need to login to DNAnexus and create a new folder for your EPACTS scripts. 

### Prerequisites

Docker image available on docker hub: [tmajarian/epacts_rmkl:1.1](https://hub.docker.com/r/tmajarian/epacts_rmkl/)

Tutorial scripts available on [GitHub](https://github.com/manning-lab/topmed-t2d-glycemia-public/tree/dev/workflows/epacts)

You do not need to download the docker image to your local machine, although it can be useful for reproducing the compute environment in DNAnexus while testing. Download the tutorial scripts and upload to your EPACTS scripts folder.

### Writing scripts for DNAnexus

Running EPACTS on DNAnexus boils down to copying your normal command line input to a shell script. The script is then interpreted within a docker image, and the EPACTS commands executed. For more information on running EPACTS from the commend line, see their [documentation](https://genome.sph.umich.edu/wiki/EPACTS). From here, we'll assume that you know the basics of using EPACTS. 

Here is a very short example that requires no inputs and just mimics the call from the EPACTS wiki:

```bash
# test_run_epacts.sh
/EPACTS_install/EPACTS-3.2.6/bin/test_run_epacts.sh
```

A more complex example (where we need to specify inputs that live in some folder in DNAnexus) is below. Notice that we have to know the names of the input files prior to writing our shell script. Since DNAnexus imports files to the working directory inside the docker image, we do not need to prepend any paths to input files. (Note: within the docker image, you can use $runepacts for your epacts command.)

```bash
# run_firth_chr10small.sh
touch freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi
/EPACTS_install/EPACTS-3.2.6/bin/epacts single --vcf freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz --ped AFEU.Firth.24AUG2017.ped --pheno t2d_ctrl --cov age --cov sex --cov PC1 --test b.firth --out testsmallminmaf001 --chr 10 -run 4 --min-maf 0.001

```

In DNAnexus, our input files would then be:
* run_firth_chr10small.sh
* freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz
* freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi
* AFEU.Firth.24AUG2017.ped

## Running the tutorial scripts from the GUI

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

## Running jobs from the command line

If you have yet to do so, follow these [instructions](https://wiki.dnanexus.com/Command-Line-Client/Quickstart) for DNAnexus command line tools installation.

### Scripts with inputs

Let's use the script shown in the __Writing scripts for DNAnexus__ section, reproduced below.

```bash
# run_firth_chr10small.sh
touch freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi
/EPACTS_install/EPACTS-3.2.6/bin/epacts single --vcf freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz --ped AFEU.Firth.24AUG2017.ped --pheno t2d_ctrl --cov age --cov sex --cov PC1 --test b.firth --out testsmallminmaf001 --chr 10 -run 4 --min-maf 0.001

```

First, we need to upload our shell script to our DNAnexus project. See this [documentation](https://wiki.dnanexus.com/Command-Line-Client/Quickstart#Upload) for command line uploading.

Once the script is uploaded, cd to the correct directory using the *dx cd* command. Next, we want to see if the Swiss Army Knife app is installed.

`dx find apps --installed`

If the Swiss Army Knife app is not within the list printed to the console:

`dx install swiss-army-knife`

Now lets check out the app:

`dx run swiss-army-knife -h`

```
usage: dx run swiss-army-knife [-iINPUT_NAME=VALUE ...]

App: Swiss Army Knife

A multi-purpose tool for all your basic analysis needs

See the app page for more information:
  https://platform.dnanexus.com/app/swiss-army-knife

Inputs:
  Input files: [-iin=(file) [-iin=... [...]]]

  Command line: -icmd=(string)

  Optional Docker image identifier: [-iimage=(string)]
        Instead of using the default Ubuntu 14.04 environment, the input command will be run using
        the specified Docker image as it would be when running 'docker run image cmd'. Example
        images identifiers are 'ubuntu:16.04', 'quay.io/ucsc_cgl/samtools'.

Outputs:
	Output files: [out (array:file)]
```

So we need three things to run our job:

1. a list of input files
2. a command
3. a docker image identifier

For the job that we want to run we have:

1. run_firth_chr10small.sh, freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz, freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi, AFEU.Firth.24AUG2017.ped
2. sh run_firth_chr10small.sh
3. tmajarian/epacts_rmkl:1.1

Then our call to DNAnexus is (assuming that all of our inputs are in our current working directory, otherwise use full paths relative to the DNAnexus project):

`dx run swiss-army-knife -iin=run_firth_chr10small.sh -iin=freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz -iin=freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi -iin=AFEU.Firth.24AUG2017.ped -icmd=sh\ run_firth_chr10small.sh -iimage=tmajarian/epacts_rmkl:1.1`

```
dx run swiss-army-knife -iin=run_firth_chr10small.sh -iin=freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz -iin=freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi -iin=AFEU.Firth.24AUG2017.ped -icmd=sh\ run_firth_chr10small.sh -iimage=tmajarian/epacts_rmkl:1.1

Using input JSON:
{
    "cmd": "sh run_firth_chr10small.sh",
    "image": "tmajarian/epacts_rmkl:1.1",
    "in": [
        {
            "$dnanexus_link": {
                "project": "project-F292q6001JK51p1PJ64F8fb0",
                "id": "file-F716F6001JKFBvbVBpV6xjvG"
            }
        },
        {
            "$dnanexus_link": {
                "project": "project-F292q6001JK51p1PJ64F8fb0",
                "id": "file-F70qbJ806fXJZ3bg6j7pyV29"
            }
        },
        {
            "$dnanexus_link": {
                "project": "project-F292q6001JK51p1PJ64F8fb0",
                "id": "file-F70qbJ806fXJZ3O66j7dfD17"
            }
        },
        {
            "$dnanexus_link": {
                "project": "project-F292q6001JK51p1PJ64F8fb0",
                "id": "file-F70qbJ806fX55P4q6g9vzXPk"
            }
        }
    ]
}

Confirm running the executable with this input [Y/n]:
```

Confirm and that job will be submitted. Now we just wait to gather our results.


## Running EPACTS with the Cloud Workstation App

**_Note that most of this information is adapted from the dnanexus [documentation](https://wiki.dnanexus.com/developer-tutorials/cloud-workstations)._**

### Motivation

Most EPACTS jobs can easily be run from the command line. To avoid downloading and uploading many large files, a user may want to run EPACTS from within the DNAnexus infrustructure directly. The cloud-workstation app allows us access to this infrustructure from a familiar command line interface without running anything locally -- so no downloading files or hogging all of your local computer's memory or storage. As long as the files that we want to work on exist within a DNAnexus project, we can run EPACTS commands.

### Prerequisites

#### Install the app

First, install the cloud workstation app

`dx install cloud_workstation`

#### Configure SSH

First, configure SSH for your account. 

`dx ssh_config`

Follow the prompts and look [here](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#ssh_config) for more info.

### Running the app

#### Checkout what inputs we might need to run the app

`dx run cloud_workstation -h`

```bash
usage: dx run cloud_workstation [-iINPUT_NAME=VALUE ...]

App: Cloud Workstation

This app sets up a cloud workstation which you can access by running the applet with the --ssh or
--allow-ssh flags

See the app page for more information:
  https://platform.dnanexus.com/app/cloud_workstation

Inputs:
  Maximum Session Length (suffixes allowed: s, m, h, d, w, M, y): [-imax_session_length=(string, default="1h")]
        The maximum length of time to keep the workstation running.  Value should include units of
        either s, m, h, d, w, M, y for seconds, minutes, hours, days, weeks, months, or years
        respectively.

  Files: [-ifids=(file) [-ifids=... [...]]]
        An optional list of files to download to the cloud workstation on startup.

  Snapshot: [-isnapshot=(file)]
        An optional list of snapshots to restore to the workshop.

Outputs: <none>
```

Since we would like this to be an interactive session (where we can directly access the command line from within the DNAnexus framework), use the *-ssh* flag

`dx run cloud_workstation --ssh`

```bash
Select an optional parameter to set by its # (^D or <ENTER> to finish):

 [0] Maximum Session Length (suffixes allowed: s, m, h, d, w, M, y) (max_session_length) [default="1h"]
 [1] Files (fids)
 [2] Snapshot (snapshot)

 ```

 Here is where you would set the job time limit.

```bash
Select an optional parameter to set by its # (^D or <ENTER> to finish):

 [0] Maximum Session Length (suffixes allowed: s, m, h, d, w, M, y) (max_session_length) [default="1h"]
 [1] Files (fids)
 [2] Snapshot (snapshot)

Optional param #: 0

Input:   Maximum Session Length (suffixes allowed: s, m, h, d, w, M, y) (max_session_length)
Class:   string

Enter string value ('?' for more options)
max_session_length: 1h

Select an optional parameter to set by its # (^D or <ENTER> to finish):

 [0] Maximum Session Length (suffixes allowed: s, m, h, d, w, M, y) (max_session_length) [="1h"]
 [1] Files (fids)
 [2] Snapshot (snapshot)

Optional param #:

Using input JSON:
{
    "max_session_length": "1h"
}

Confirm running the executable with this input [Y/n]:
```

After conformation, we'll be connected to a virtual machine within DNAnexus 

```bash
Welcome to DNAnexus!

This is the DNAnexus Execution Environment, running job-F79kYj801JK9fKP4KF3B2vB0.
Job: Cloud Workstation
App: cloud_workstation:main
Instance type: mem1_ssd1_x4
Project: cdbg_t2d (project-F292q6001JK51p1PJ64F8fb0)
Workspace: container-F79kYjQ034ybKvpb61q9yVJb
Running since: Tue Oct  3 14:13:27 UTC 2017
Running for: 0:00:43
The public address of this instance is ec2-34-228-254-97.compute-1.amazonaws.com.
You are running byobu, a terminal session manager.
If you get disconnected from this instance, you can log in again; your work will be saved as long as the job is running.
For more information on byobu, press F1.
The job is running in terminal 1. To switch to it, use the F4 key (fn+F4 on Macs; press F4 again to switch back to this terminal).
Use sudo to run administrative commands.
From this window, you can:
 - Use the DNAnexus API with dx
 - Monitor processes on the worker with htop
 - Install packages with apt-get install
 - Use this instance as a general-purpose Linux workstation
OS version: Ubuntu 14.04.5 LTS (GNU/Linux 3.13.0-117-generic x86_64)
dnanexus@job-F79kYj801JK9fKP4KF3B2vB0:~⟫
```

#### Setting up EPACTS environment

We can use docker to setup our compute environment. This will import all of the EPACTS code and allow us to run some commands. First pull the image from the [remote repository](https://hub.docker.com/r/tmajarian/epacts_rmkl/).

`dx-docker pull tmajairan/epacts_rmkl:1.1`

Next, we'll want to import any files from our DNAnexus project that we might want to use. If we don't know the exact path to a file, we can view the file structure of our project

`dx ls $DX_PROJECT_CONTEXT_ID`

```bash
Annotation/
Chr10_Pilot_Rare/
freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz
freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi
```

Now we can import the files we need

```bash
dx download $DX_PROJECT_CONTEXT_ID:/freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz
dx download $DX_PROJECT_CONTEXT_ID:/freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi
dx download $DX_PROJECT_CONTEXT_ID:/AFEU.Firth.24AUG2017.ped
```

Now that we have the correct input files downloaded, we'll want to run our docker image. Be sure to include all input files that you might need when running the *dx-docker run* command. We are essentially creating a self contained virtual machine (within our current virtual machine) that has no outside access. Any files not imported will not be visible from within the docker container. More information on this can be found [here](https://docs.docker.com/engine/admin/volumes/volumes/)

```bash
# first create and move the downloaded files to their own folder
mkdir data
mv freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz ./data/
mv freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi ./data/
mv AFEU.Firth.24AUG2017.ped ./data/
dx-docker run -v data:/data tmajarian/epacts_rmkl:1.1
```

Once we're within the docker container, we can *ls* and see our local file system. 
(highlight data)
```bash
EPACTS_install  bin  boot  data  dev  etc  gs_install  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
```

Notice that the folder `data` is exactly what we imported using the `-v` command above.

`ls ./data`

```bash
freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi AFEU.Firth.24AUG2017.ped
```

#### Running an epacts command

Just like our single variant workflow script (above), we can run the command directly from the command line in the docker container.

```bash
cd data/
touch freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz.tbi
/EPACTS_install/EPACTS-3.2.6/bin/epacts single --vcf freeze4.chr10.pass.gtonly.minDP10.genotypes.small.vcf.gz --ped AFEU.Firth.24AUG2017.ped --pheno t2d_ctrl --cov age --cov sex --cov PC1 --test b.score --out testsmallminmaf001 --chr 10 -run 4 --min-maf 0.001
```

#### Pushing result files back to the DNAnexus project

To finish our workflow, we still need to make sure our generated results end up back in our project in DNAnexus. We'll need to upload them back to the project directory.

`dx upload --path "$DX_PROJECT_CONTEXT_ID:" /data/out/*.png`

Check that the files were actuall uploaded.

`dx ls`

### Ending the job

Make sure to end your session when everything is finished.

`dx terminate $DX_JOB_ID`

## Other resources

For more information on the cloud-workstation app, see [this page](https://wiki.dnanexus.com/developer-tutorials/cloud-workstations).

For information on running docker commands, see [this page](https://docs.docker.com/)






