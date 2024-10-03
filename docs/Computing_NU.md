# Northeastern University Computing Resources 

A high-performance computing resource for the Northeastern research community, the Discovery cluster is located in the [Massachusetts Green High Performance Computing Center](https://www.mghpcc.org/) in Holyoke, Massachusetts.

The Discovery cluster provides Northeastern researchers with access to more than 45,000 CPU cores and more than 400 GPUs. Connected to the university network over 10 Gbps Ethernet (GbE) for high-speed data transfer, Discovery provides 5 PB of available storage on a high-performance GPFS parallel file system. Compute nodes are connected with either 10 GbE or a high-performance HDR100 InfiniBand (IB) interconnect running at 200 Gbps, supporting all types and scales of computational workloads.

As GMGI's researchers, we have access to Northeastern's HPC resources through an MOU established in Fall 2023. 

Read the HPC resource documentation prior to getting started: [Research Computing - NURC RTD (northeastern.edu)](https://rc-docs.northeastern.edu/en/latest/welcome/index.html).

![](https://github.com/GMGI-Fisheries/resources/blob/master/img/NU_computing_resources.png?raw=true)

## Logging in

Before creating an account with the NU Discovery Cluster, claim your Northeastern email and sponsored account.

1. Log into (with the northeastern email and pw previously claimed): [Home - Northeastern Tech Service Portal](https://service.northeastern.edu/tech?id=tech_index_home).  
2. Navigate to [High Performance Computing - Northeastern Tech Service Portal](https://service.northeastern.edu/tech?id=sc_category&sys_id=43a3aef7db45cdd0ca10819b13961998).    
3. Request an account ([Research Computing Access Request - Northeastern Tech Service Portal](https://service.northeastern.edu/tech?id=sc_cat_item&sys_id=0ae24596db535fc075892f17d496199c)). Fill out the form with your Northeastern email, and following options:  
- Select: I do not have access to Discovery. I am requesting a new account.  
- Affiliation with Northeastern University: Visiting Researcher (Greg says this answer doesn't really matter).   
- University Sponsor: Geoffrey Trusell  
- Gaussian: No (This is a specific program, if you don't know what it is, you don't need it).  
- Select: By clicking here, I acknowledge that I have read and agree to the following.  
4. This triggers an email to Geoff Trusell to approve your account. Send an email to Geoff letting him know that you are activating your account that needs his approval.  
5. Once Geoff approves the account sponsorship, then Greg and the computing team will finish setting up your account.

You can operate on the Discovery Cluster in two ways:  
1. via Linux operating system on your computer or ssh client  
2. NU's Open On Demand ([Open OnDemand (OOD) - RC RTD (northeastern.edu)](https://rc-docs.northeastern.edu/en/latest/using-ood/index.html)) GUI. [Accessing Open OnDemand - RC RTD (northeastern.edu)](https://rc-docs.northeastern.edu/en/latest/using-ood/accessingood.html#access-ood). With OOD interface, you can access plug-ins (i.e. RStudio) and launch the server from the web. Serena recommended using incognito window b/c OOD usually works better without the caching.

```
ssh username@login.discovery.neu.edu
```

## Server structure

We have a 30 TB maximum in our working space `/work/gmgi`:  
- Each lab has their own subfolder that serves as their storage and working space (e.g., `work/gmgi/Fisheries/`).  
- `databases/`: Shared folder for common databases for amplicon sequencing (i.e., 12S, 16S, 18S, COI) and NCBI nt database. View the README file for databases sources.
- `containers/`: Shared folder for custom built containers.  
- `packages/`: Shared folder for programs downloaded for all users.  
- `miniconda3/` and `Miniconda3-latest-Linux-x86_64.sh`: Shared resource for building conda environments. Environments built here can be used by everyone. *Do not edit.*  
- `check_storage.sh`: Bash script built to calculate TB usage from each lab's folder and gmgi's work and output is `storage_summary_2024-09-23.txt` with the date calculated.

```
[e.strand@login-00 ~]$ cd /work/gmgi
[e.strand@login-00 gmgi]$ ls
check_storage.sh  containers  databases  ecosystem-diversity  Fisheries  miniconda3  Miniconda3-latest-Linux-x86_64.sh  packages  storage_summary_2024-09-23.txt
```

Fisheries folder (`work/gmgi/Fisheries/`) is split by the type of project. `reference_genomes` includes .fasta reference files for organisms rather than a database (`Haddock_ref.fasta`). 

```
[e.strand@login-00 Fisheries]$ ls
eDNA  epiage  reference_genomes
```

## Storage rules

While analyzing data, NOT just at the end of a project!

Raw data files are backed up on AWS services and on GMGI RHEL Gadus immediately upon receiving data. If working on NU cluster, once user is happy with data analysis pipeline, raw and final data is to be removed from NU cluster and only kept on AWS services or GMGI's RHEL server. 

Compressing files:  
- Gzip all fastq files (e.g., raw data, trimmed data), .fasta/.fa files (e.g., reference genomes), and large .txt files (e.g., intermediate files created during analysis): `gzip *.fastq` or create a slurm array with a sbatch script.    
- [Genozip](https://www.genozip.com/standard) all .bam, .sam, .vcf files (e.g., intermediate files created during analysis). Genozip has been downloaded in /work/gmgi/packages for general use.   

Space-related commands:  
- List all files within a directory and human-readable sizes (folder size is not total size of folder contents): `ls -lha`  
- Calculate total storage taken up by one directory (change path as needed): `du -shc .[^.]* /work/gmgi/fisheries`    
- In /work/gmgi/, there is a `check_storage.sh` bash script that will use the above commands to create a summary .txt file with the storage use of each team.  

## NU Contacts and Research Computing Help

GMGI's two main contacts are: Greg Shomo (g.shomo@northeastern.edu) and Serena Caplins (s.caplins@northeastern.edu). 

If you need assistance, NU's support team is available at rchelp@northeastern.edu or consult the [Frequently Asked Questions (FAQs)](https://rc-docs.northeastern.edu/en/latest/faq.html#faq). Emailing the rchelp team will create a help ticket one of NU's team members will be assigned to your case. When creating a help ticket, CC Geoff Trussell (g.trussell@northeastern.edu) and Jon Grabowski (j.grabowski@northeastern.edu). 

RC Help hosts office hours to connect with RC staff and Graduate Research Assistants (GRAs) to ask questions about any RC-related questions that you might have (e.g., Setting up conda environments, Installing software, Optimizing the runtime of your sbatch scripts, Effectively using the Open onDemand website):  
- Wednesdays 3 pm - 4 pm: [zoom link](https://url2.mailanyone.net/scanner?m=1rzNQw-0008Ns-3f&d=4%7Cmail%2F90%2F1713906600%2F1rzNQw-0008Ns-3f%7Cin2e%7C57e1b682%7C28509242%7C14152682%7C6628242AF22ECB27DE9BAA81E1B29580&o=%2Fphto%3A%2Fntseertnstrhaso.zj.u%2Fom20%2F951142466&s=ruDKyxUNCNd62B7o2lO1X9v5I5U)  
- Thursdays 11 am - 12 pm: [zoom link](https://url2.mailanyone.net/scanner?m=1rzNQw-0008Ns-3f&d=4%7Cmail%2F90%2F1713906600%2F1rzNQw-0008Ns-3f%7Cin2e%7C57e1b682%7C28509242%7C14152682%7C6628242AF22ECB27DE9BAA81E1B29580&o=%2Fphto%3A%2Fntseertnstrhaso.zj.u%2Fom63%2F914254083&s=feG6PL708qP4uNB_somz55Vucz0)  

Ask 2-3 others at GMGI for assistance and/or attend office hours prior to creating a help ticket. We do not want to overwhelm the rchelp desk if we can troubleshoot internally first. GMGI has a #bioinformatics slack channel for this purpose.

## Running a bioinformatic script 

Read through: [Running Jobs - RC RTD (northeastern.edu)](https://rc-docs.northeastern.edu/en/latest/runningjobs/index.html) before starting. 

Jobs are run either through:  
1. Interactive mode (immediate execution and feedback): `srun --pty bash` to claim a node and then utilize `bash scriptname.sh` to run a script.    
2. Batch jobs: using scripts to manage longer-running jobs: `sbatch scriptname.sh` to run a script.  

Interactive mode would be equivalent to running a job directly in your terminal window without starting a tmux session. Using batch jobs and shell scripts would be similar to tmux session where you can turn off wifi, walk away, etc. Interactive requires you to stay connected. 

Introduction to Slurm scripts:    
- [Slurm - RC RTD (northeastern.edu)](https://rc-docs.northeastern.edu/en/latest/slurmguide/index.html)     
- [Best Practices - RC RTD (northeastern.edu)](https://rc-docs.northeastern.edu/en/latest/best-practices/index.html)    

## Packages and modules 

NU has some modules downloaded that are accessible for all users. Otherwise, larger packages should be installed in a conda environment by the user or simple packages can be downloaded to the `/work/gmgi/packages/` folder. Instructions for downloading to conda env: [https://rc-docs.northeastern.edu/en/latest/software/index.html]

Common commands:   
- To find already installed programs: `module avail`      
- To get information about a module: `module help [module/version]` or `module whatis [module/version]`. "help" will provide what the module is, package information including version and install date, and a link to the documentation/github. "whatis" will provide a short, one line description of the program.    
- To load a module: `module load [module/version]` (e.g., `module load bamUtil/v1.0.15`). Loading a module will put all the necessary executables and dependencies for that program in your path so you can call the commands from any location (i.e. your working directory).   