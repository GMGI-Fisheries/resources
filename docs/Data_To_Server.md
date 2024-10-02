# Downloading sequencing data to servers

Goal: Download .fastq files from sequencing center to HPC.

#### Table of Contents
- [GMGI in-house sequencing to HPC](#Illumina BaseSpace to NU Discovery Cluster or GMGI in-house HPC)
- [External sequencing to HPC](#Globus to NU Discovery Cluster or GMGI in-house HPC-content)
- [HPC to HPC](#GMGI in-house to NU HPC)

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

GMGI endpoint set-up (if using NU cluster, skip to next section):  
1. 


## GMGI in-house to NU HPC 

