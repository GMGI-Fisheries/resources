# eDNA Metabarcoding for COI target region

The COI region is commonly used for metabarcoding practices and consequently there are many primer options to choose from. The Fisheries team at GMGI has optimized the Leray Geller set (outlined in red box below). Citation: [Leray et al 2013](https://link.springer.com/article/10.1186/1742-9994-10-34).

We primarily use this set for invertebrate targets and 12S for vertebrate communities. 

![](https://github.com/GMGI-Fisheries/resources/blob/master/img/eDNA_meta_COI_primerset.png?raw=true)

Workflow done on HPC. Scripts to run: 

1. Confirm conda environment is available      
2. Assess quality of raw data (00-fastqc.sh)        
3. Visualize quality of raw data (00-multiqc.sh)     
4. Set-up Mothur program (01-Mothur-setup.sh)  
5. QC'ing seqs via Mothur (02-Mothur-QC.sh)      
6. Determining and counting unique sequences with Mothur (03-Mothur-unique.sh)   
7. Taxonomic Assignment with Mothur (04-Mothur-tax.sh)  

*Descriptions of Mothur steps are from [A. Huffmyer](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/_posts/2022-01-12-16S-Analysis-in-Mothr-Part-1.md#makecontigs) and [E. Strand](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2022-01-24-KBay-Bleached-Pairs-16S-Analysis-Mothur.md) Mothur for 16S notebook posts.* 

[Mothur program page](https://mothur.org/)

## Step 1: Confirm conda environment is available 

The conda environment is started within each slurm script, but to activate conda environment outside of the slurm script to update packages or check what is installed:

```
# Activate conda
source /work/gmgi/miniconda3/bin/activate

# Activate fisheries eDNA conda environment 
conda activate fisheries_eDNA
conda activate eDNA_COI

# List all available environments 
conda env list 

# List all packages installed in fisheries_eDNA
conda list 

# Update a package
conda update [package name]
``` 

## Step 2: Assess quality of raw data  

Background information on [FASTQC](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/05_qc_running_fastqc_interactively.html). 

`00-fastqc.sh`: 

```
#!/bin/bash
#SBATCH --error=output/fastqc_output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/fastqc_output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=fastqc
#SBATCH --mem=3GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for your project

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate fisheries_eDNA

## SET PATHS 
raw_path=""
out_dir=""

## CREATE SAMPLE LIST FOR SLURM ARRAY
### 1. Create list of all .gz files in raw data path
ls -d ${raw_path}/*.gz > ${raw_path}/rawdata

### 2. Create a list of filenames based on that list created in step 1
mapfile -t FILENAMES < ${raw_path}/rawdata

### 3. Create variable i that will assign each row of FILENAMES to a task ID
i=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

## RUN FASTQC PROGRAM 
fastqc ${i} --outdir ${out_dir}
```

To run:    
- Start slurm array (e.g., with 138 files) = `sbatch --array=0-137 00-fastqc.sh`.

Notes:  

- This is going to output *many* error and output files. After job completes, use `cat *output.* > ../fastqc_output.txt` to create one file with all the output and `cat *error.* > ../fastqc_error.txt` to create one file with all of the error message outputs. 
- Within the `out_dir` output folder, use `ls *html | wc` to count the number of html output files (1st/2nd column values). This should be equal to the --array range used and the number of raw data files. If not, the script missed some input files so address this before moving on.  

## Step 3: Visualize quality of raw data  

Background information on [MULTIQC](https://multiqc.info/docs/#:~:text=MultiQC%20is%20a%20reporting%20tool%20that%20parses%20results,experiments%20containing%20multiple%20samples%20and%20multiple%20analysis%20steps).

`00-multiqc.sh` 

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=multiqc
#SBATCH --mem=8GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for your project
## 2. Optional: change file name (multiqc_raw.html) as desired

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate fisheries_eDNA

## SET PATHS 
## fastqc_output = output from 00-fastqc.sh; fastqc program
fastqc_output="" 
multiqc_dir="" 

## RUN MULTIQC 
multiqc --interactive ${fastqc_output} -o ${multiqc_dir} --filename multiqc_raw.html
```

To run:  
- `sbatch 00-multiqc.sh` 

Notes:  

- Depending on the number of files per project, multiqc can be quick to run without a slurm script. To do this, run each line separately in the command line after activating the conda environment.  


## Step 4: Set-up Mothur program 

#### Create primer file

Project set-up create primer seq file, **only needed once, if file is already there do not create a new one**:  
1. Navigate to COI databases folder: `cd /work/gmgi/databases/COI`  
2. Create new file for Leray Gellar primer set: `nano mothur_oligos_LG`  
3. Copy and paste the below content:

```
## Leray Geller 2013 primers 
## format = primer F R

primer GGWACWGGWTGAACWGTWTAYCCYCC TAIACYTCIGGRTGICCRAARAAYCA
```

#### Mothur program set-up project 

Create directory: `mkdir Mothur_data`

`make.file()`: tells mothur to look for fastq files in the data directory (`rawdata_dir`) and identify forward and reverse reads. Put type=gz to look for gz files. If there are .fasta or .fastq files you can use type=fasta or type=fastq. The prefix gives the output file a name - this can be a project name. This creates a file called proj_name.files with sample name, R1, R2 listed in columns.  
- Output: `proj.paired.files` within the `rawdata_dir`  

`make.contigs()`: makes contigs of forward and reverse reads for each sample. This will assemble contigs for each pair of files and use R1 and R2 quality information to correct any sequencing errors in locations where there is overlap and the quality of the complementary sequence is sufficient. Because we added an oligos file, make.contigs will remove the primers on these sequences.   
- `pdiffs=5`: primer differences=5. Try various parameters from 2-10 (based on other pipeline I've used or seen in the past).  
- `checkorient=t`: Degenerate primers can be used in our oligos file. We are also using check orient to allow Mothur to "flip" the reverse primer if it is not found.  
- `trimoverlap=F`: TRUE is needed if sequencing length was longer than amplicon length. The COI region is 313 bp which is longer than 2x250 bp sequencing.  

`summary.seqs()`: generates summary information about the sequences from the files we made above.

Run time: ~30-45 min for 30-40 files (from previous project; update values with GMGI COI info)

`01-Mothur-setup.sh`

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=mothur_setup
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=3

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate eDNA_COI

# Set Paths 
oligo_file="/work/gmgi/databases/COI/mothur_oligos_LG"

## User specific paths
rawdata_dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/raw_data"
proj_name="OSW_2023_invert" 
output_dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/Mothur_data"

## Make symlinks from raw data with new name (_ instead of -)
for file in ${rawdata_dir}/*.gz; do
    filename=$(basename "$file")
    newname=${filename//-/_}
    ln -s "$file" "${output_dir}/${newname}"
done

cd ${output_dir} 

## Make file with raw fastq files
mothur "#make.file(inputdir=., type=gz, prefix=${proj_name})"   

## Make contig file 
mothur "#make.contigs(inputdir=., outputdir=., file=${proj_name}.paired.files, trimoverlap=F, oligos=${oligo_file}, pdiffs=5, checkorient=T)"

## Create a summary file with trimmed contigs 
mothur "#summary.seqs(fasta=${proj_name}.paired.trim.contigs.fasta)"

## Count the number of sequences that were removed and the number that were kept by counting sequences in each fasta file
## This will be in the output file created by the .sh script 
grep -c "^>" ${proj_name}.paired.trim.contigs.fasta
grep -c "^>" ${proj_name}.paired.scrap.contigs.fasta 
```

The logfiles will be named with the run ID so I used `mv` to change these to names that reflect the step:  
- `mothur.make.contigs.logfile`  
- `mothur.summaryseqs.trimmed.logfile`  

Output:    
- proj_name.contigs.groups                 
- proj_name.contigs.report                                        
- proj_name.scrap.contigs.fasta            
- proj_name.trim.contigs.fasta 
- proj_name.trim.contigs.summary  

Descriptions of contig files:   
- Trim file = sequences that were "good".  
- Scrap file = sequences that were "bad".  
- Groups file = what group each sequence belongs to map sequence to each sample from the trimmed sequence file.  
- Contigs report file = information on sequences that were aligned and paired together.  

View the output file (`less mothur_setup.output*`) to check how many sequences were kept/lost with certain parameters. Troubleshoot this as needed. Shift+G to reach the end of the output file and view summary table (example below): 

```
## OSW data example - room for troubleshooting here 
10,430,038 reads in the trimmed fasta 
2,876,813 reads in the scrapped fasta 

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       124     124     0       3       1
2.5%-tile:      1       251     251     0       4       260751
25%-tile:       1       365     365     0       5       2607510
Median:         1       365     365     0       6       5215020
75%-tile:       1       365     365     1       6       7822529
97.5%-tile:     1       404     404     8       8       10169288
Maximum:        1       502     502     191     233     10430038
Mean:   1       363     363     1       5
# of Seqs:      10430038
```

This table shows quantile values about the distribution of sequences for a few things:   
- Start position: All at 1 now, will start at different point after some QC.  
- End position: We see that there are some sequences that are very short and we may need to remove those later.  
- Number of bases: length. we see most are in expected range here, but one is super long! This might tell us there is no overlap so they are butted up against each other. We will remove things like this.  
- Ambigs: Number of ambiguous calls in sequences. Here there are a few that have ambiguous base calls. We will remove any sequence with an ambiguous call or any longer than we would expect for V4 region.  
- Polymer: Length of polymer repeats.  
- NumSeqs: Number of sequences.  

View contigs report for quality assessment (`less proj_name.contigs.report`)  
 
**Conduct quality assessment before proceeding with next steps.** 

## Step 5: QC'ing seqs via Mothur 

Additional QC steps besides removing primers. 

Check the 25%-tile, median, and 75%-tile values from the table above. That will determine the max cut-offs below along with the expected length of the COI region (~313 bp but in OSW data, the majority are >350 bp). 

`screen.seqs()`: specify the fasta file of the contigs generated in the previous step and remove any sequence with an ambiguous call ("N"). We will also remove sequences >350 nt. We will also set a minimum size amount (200). These parameters could be adjusted based on specific experiment and variable region.

`02-Mothur-QC.sh`

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=mothur_QC
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate eDNA_COI

dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/Mothur_data"
proj_name="OSW_2023_invert" 

cd ${dir}

mothur "#screen.seqs(inputdir=., outputdir=., fasta=${proj_name}.paired.trim.contigs.fasta, group=${proj_name}.paired.contigs.groups, maxambig=0, maxlength=400, minlength=200)"

mothur "#summary.seqs(fasta=${proj_name}.paired.trim.contigs.good.fasta)"
```

Output files:    
- proj_name.trim.contigs.good.fasta   
- proj_name.trim.contigs.bad.accnos    
- proj_name.contigs.good.groups  

OSW example: Removed 3,652,435 sequences that fall outside of 200-400 bp length bin with no ambiguous calls so 6,777,603 remaining sequences. 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       213     213     0       3       1
2.5%-tile:      1       251     251     0       4       169441
25%-tile:       1       365     365     0       5       1694401
Median:         1       365     365     0       6       3388802
75%-tile:       1       365     365     0       6       5083203
97.5%-tile:     1       368     368     0       7       6608163
Maximum:        1       400     400     0       105     6777603
Mean:   1       361     361     0       5
# of Seqs:      6777603
```

Changed the name of the log files to:    
- mothur.screen.seqs.summary.logfile   
- mothur.screen.seqs.logfile  

## Step 6: Determining and counting unique sequences with Mothur

Identifying unique sequences and counting those. 

`03-Mothur-unique.sh`

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=mothur_uniq
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate eDNA_COI

dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/Mothur_data"
proj_name="OSW_2023_invert" 

cd ${dir}

mothur "#unique.seqs(fasta=${proj_name}.paired.trim.contigs.good.fasta)"

mothur "#count.seqs(name=${proj_name}.paired.trim.contigs.good.names, group=${proj_name}.paired.contigs.good.groups)"

mothur "#summary.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.fasta, count=${proj_name}.paired.trim.contigs.good.count_table)"
```

OSW example output: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       213     213     0       3       1
2.5%-tile:      1       251     251     0       4       169441
25%-tile:       1       365     365     0       5       1694401
Median:         1       365     365     0       6       3388802
75%-tile:       1       365     365     0       6       5083203
97.5%-tile:     1       368     368     0       7       6608163
Maximum:        1       400     400     0       105     6777603
Mean:   1       361     361     0       5
# of unique seqs:       3169531
total # of seqs:        6777603
```

There are 3,169,531 unique sequences in this OSW example dataset. 

## Step 7: Taxonomic Assignment with Mothur 

Using MetaZooGene Global database as reference. This database has sequences from both BOLD and NCBI GenBank. 

https://metazoogene.org/mzgdb/. Navigate to the [Worlds oceans](https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/index-o00.html) page and click MZGdb Data Access page. We have downloaded both the Global and North Atlantic versions.   

Download and update database:

```
cd /work/gmgi/databases/COI/MetaZooGene

wget https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/data-src/MZGfasta-coi__T4000000__o00__A.fasta.gz
wget https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/data-src/MZGfasta-coi__T4000000__o02__A.fasta.gz

## unzip into new name file 
## replace with version number 
gunzip -c MZGfasta-coi__T4000000__o00__A.fasta.gz > MZG_v2023-m07-15_Global_modeA.fasta
gunzip -c MZGfasta-coi__T4000000__o02__A.fasta.gz > MZG_v2023-m07-15_NorthAtlantic_modeA.fasta
```

Align the database prior to using as input in Mothur:

`MZG_align.sh`

```
#!/bin/bash
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=MZG_MAFFT
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### changing name format to MAFFT_xxx_aligned.fasta 

# Assign input and output file names from command-line arguments
INPUT_FILE=$1
OUTPUT_FILE=$2

# Run MAFFT with the provided input and output files
mafft --auto "$INPUT_FILE" > "$OUTPUT_FILE"
```

Use format `sbatch MZG_align.sh input output`    
Example: `sbatch MZG_align.sh MZG_v2023-m07-15_NorthAtlantic_modeA.fasta MAFFT_MZG_v2023-m07-15_NorthAtlantic_modeA_aligned.fasta`   


Running taxonomic assignment via Mothur: 

`04-Mothur-tax.sh`

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=mothur_tax
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate eDNA_COI

dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/Mothur_data"
proj_name="OSW_2023_invert" 
ref="/work/gmgi/databases/COI/MAFFT_MZG_v2023-m07-15_NorthAtlantic_modeA_aligned.fasta"

cd ${dir}

mothur "#align.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.fasta, reference=${ref}, flip=T, processors=48)"
mothur "#summary.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.align, count=${proj_name}.paired.trim.contigs.good.count_table)"
```

## Step 8: QC Taxonomic Alignment 

`05-Mothur-taxQC.sh`

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=mothur_taxQC
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate eDNA_COI

dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/Mothur_data"
proj_name="OSW_2023_invert" 

cd ${dir}

## Identifying those that don't meet the criteria 

mothur "#screen.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.align, count=${proj_name}.paired.trim.contigs.good.count_table, start=1968, end=11550, maxhomop=8)"
mothur "#summary.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.align, count=${proj_name}.paired.trim.contigs.good.count_table)"
mothur "#count.groups(count=${proj_name}.paired.trim.contigs.good.count_table)"

## Filtering those out
mothur "#filter.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.align, vertical=T, trump=.)"
mothur "#summary.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.good.filter.fasta, count=${proj_name}.paired.trim.contigs.good.count_table)"

## ID unique sequences again
mothur "#unique.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.good.filter.fasta, count=${proj_name}.paired.trim.contigs.good.good.count_table)"
```

## Step 9: Pre-clustering aka denoising

Now we need to further polish and cluster the data with pre.cluster. The purpose of this step is to remove noise due to sequencing error. The rational behind this step assumes that the most abundant sequences are the most trustworthy and likely do not have sequencing errors. Pre-clustering then looks at the relationship between abundant and rare sequences - rare sequences that are "close" (e.g., 1 nt difference) to highly abundant sequences are likely due to sequencing error. This step will pool sequences and look at the maximum differences between sequences within this group to form ASV groupings.

In this step, the number of sequences is not reduced, but they are grouped into amplicon sequence variants ASV's which reduces the error rate. 

Other programs that conduct this "denoising" are DADA2, UNOISE, and DEBLUR. However, these programs remove the rare sequences, which can distort the relative abundance of remaining sequences. DADA2 also removes all sigletons (sequences with single representation) which disproportionately affects the sequence relative abundance. Mothur avoids the removal of rare sequences for this reason.

`06-Mothur-denoise.sh`

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=mothur_taxQC
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate eDNA_COI

dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/Mothur_data"
proj_name="OSW_2023_invert" 

cd ${dir}

mothur "#pre.cluster(fasta=${proj_name}.paired.trim.contigs.good.unique.good.filter.unique.fasta, count=${proj_name}.paired.trim.contigs.good.unique.good.filter.count_table, diffs=1)"
mothur "#summary.seqs(fasta=${proj_name}.paired.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=${proj_name}.paired.trim.contigs.good.unique.good.filter.unique.precluster.count_table)"


```




Currently on:  
- creating MAFFT version of ref  
- re-run the mothur tax script 

