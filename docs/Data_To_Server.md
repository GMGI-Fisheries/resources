# Downloading sequencing data to servers

Goal: Download .fastq files from sequencing center to HPC and/or move data between HPCs and personal computers. 

#### Table of Contents
- [GMGI in-house sequencing to HPC](#illumina-basespace-to-nu-discovery-cluster-or-gmgi-in-house-hpc)  
- [External sequencing to HPC](#globus-to-nu-discovery-cluster-or-gmgi-in-house-hpc)  
- [HPC to HPC or Personal Computer](#hpc-to-hpc-or-hpc-to-personal-computer)
- [HPC to AWS back-up](#aws-back-up)

## Illumina BaseSpace to NU Discovery Cluster or GMGI in-house HPC

[Illumina BaseSpace CLI instructions](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview#InstallBaseSpaceSequenceHubCLI)

Connecting your user to Illumina BaseSpace:  

*Each user only needs to complete this once to set-up. If completed for previous projects, skip to downloading data steps.* 

1. Create folder called bin: `mkdir $HOME/bin`  
2. Download BaseSpace CLI: `wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O $HOME/bin/bs`  
3. Change the file permissions to make the downloaded binary executable: `chmod u+x $HOME/bin/bs`  
4. Authenticate your account: `bs auth`  
5. Navigate to the webpage provided and authenticate use of BaseSpace CLI.  

Download data from each run to desired output path: 

6. Find the Run ID of desired download: Within Sequence Hub, navigate to Runs and select the desired run. The Run ID is in the webpage handle (e.g., https://basespace.illumina.com/run/123456789/details). 

![](https://github.com/GMGI-Fisheries/resources/blob/master/img/Data_To_Server_basespace_runsPage.png?raw=true)

7. Navigate to `cd $HOME/bin` and download dataset: `bs download run -n run_name --extension=fastq.gz -o /local/output/path`. Replace `run_name` with the exact name of the run on BaseSpace.  
8. Navigate to the output path `cd /local/output/path` and move all files out of subdirectories: `mv */* .` 

## Globus to NU Discovery Cluster or GMGI in-house HPC 

External sequencing centers (e.g., UConn) will share data via Globus. Instructions from NU on [transfering data](https://rc-docs.northeastern.edu/en/latest/datamanagement/transferringdata.html) and using [Globus](https://rc-docs.northeastern.edu/en/latest/datamanagement/globus.html#using-globus). Globus works by transferring data between 'endpoints'. NU's endpoint is called Discovery Cluster which is searchable but our in-house GMGI endpoint needs to be created for each user. 

Globus instructions: [https://docs.globus.org/globus-connect-server/v5.4/quickstart/]. Create a Globus account prior to instructions below. 

### GMGI in-house server 

GMGI endpoint set-up (only need to do this once):  
1. Navigate to the globusconnectpersonal-3.2.2 module that is already downloaded on GMGI's in-house server: `cd /data/resources/app_modules/globusconnectpersonal-3.2.2`.  
2. Set-up an endpoint: `./globusconnectpersonal -setup --no-gui`  
3. This will then ask you to click on a log-in link. Once logged in, you receive a authorization code. Paste that in your terminal window where it asked for this code.
4. Name your endpoint with your own user (change this to your first and last name): `user.name`  
5. If successfully, globus will output: 

```
Input a value for the Endpoint Name: user.name

registered new endpoint, id: [unique ID to you]

setup completed successfully
```

Start Globus transfer:  
1. Navigate to the globusconnectpersonal-3.2.2 module that is already downloaded on GMGI's in-house server: `cd /data/resources/app_modules/globusconnectpersonal-3.2.2`.
2. Activate personal endpoint: `./globusconnectpersonal -start &`  
3. Your `user.name` endpoint will now appear as an option on the Globus online interface.   
4. Log into Globus and Navigate to 'Collections' on the left hand panel. Confirm that your GMGI endpoint is activated (green icon):

![]()

### NU Discovery Cluster 





## HPC to HPC or HPC to personal computer




## AWS Back-up

AWS is our Amazon Web Services option for long-term data storage and back-up. GMGI uploads data from our in-house server to AWS.

What should be backed up:  
- Raw data files such as fastq files directly from the sequencer  
- Final result data files (i.e. count table, assemblies, etc.)  

1. Make sure all files are compressed by:  
- Gzip all fastq files (e.g., raw data, trimmed data), .fasta/.fa files (e.g., reference genomes), and large .txt files (e.g., intermediate files created during analysis): `gzip *.fastq` or create a slurm array with a sbatch script.  
- [Genozip](https://www.genozip.com/standard) all .bam, .sam, .vcf files (e.g., intermediate files created during analysis). Genozip program is downloaded NU in `/work/gmgi/packages/` for general use.  

**Prior to AWS back-up, check with Tim or Emma for approval of files and compression.** 

2. Create a new screen session called AWS_tar (user can change this name to desired): `tmux new -s AWS_tar`  
3. Create a txt file with file sizes of all desired input: `ls -l *gz > file_size.txt`  
4. Edit this file to be only file sizes and names: `awk '{print $5,$9}' file_size.txt > file_size_edited.txt`  
5. View this edited file: `head file_size_edited.txt`

```
11400971821 Mae-263_S1_R1_001.fastq.gz

12253428145 Mae-263_S1_R2_001.fastq.gz

11962611469 Mae-266_S2_R1_001.fastq.gz

12839131166 Mae-266_S2_R2_001.fastq.gz

9691926610 Mae-274_S3_R1_001.fastq.gz
```

6. Sum the first column: `awk '{sum += $1} END {print sum}' file_size_edited.txt`. This value is in Byte (B), but convert to MB or TB for a more helpful value to work with. It's important to know the size of the data you are working for storage and cost purposes. Backing up to AWS costs $$/monthly based on TBs stored and our HPC systems have max TB storage limitations.  
7. Tar all desired data to result in on zipped file: `tar -czvf HaddockEpiAge1_rawData_20231113.tar.gz ./*fastq.gz`. Tar file name needs to be a unique identifier that includes project name, date, and description of the files included.  
8. Detach from a session: Press Ctrl+B, release, and then press D. This tar function will take awhile especially for large datasets.  
9. Reopen/attach a detached session: `tmux attach-session -t AWS_tar`.  
10. Check the tar file size: `ls -lha`. Example output: `-rw-rw-r--. 1 estrand science 1736366676523 Nov 15 17:54 HaddockEpiAge1_rawData_20231113.tar.gz`. This 1736366676523 value is the size in B which should match exactly the sum of all input file sizes calculated previously.   
11. Once this tar function is complete, end a tmux session (forever - not just detached): In the attached session, type exit and press enter. Or press Ctrl+D.  
12. Move the packed tar archive file to AWS transfer folder: `mv HaddockEpiAge1_rawData_20231113.tar.gz /data/prj/AWS-transfers`.  
13. Notify Jen that there is a transfer waiting so she can move this to AWS services and then remove from the /AWS-transfers folder. 