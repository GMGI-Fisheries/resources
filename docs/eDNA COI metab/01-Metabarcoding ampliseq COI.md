# Metabarcoding workflow for COI amplicon sequencing 

The COI region is commonly used for metabarcoding practices and consequently there are many primer options to choose from. The Fisheries team at GMGI has optimized the Leray Geller set (outlined in red box below). Citation: [Leray et al 2013](https://link.springer.com/article/10.1186/1742-9994-10-34).

We have primarily use this set for invertebrate targets and 12S for vertebrate communities. 

![](https://github.com/GMGI-Fisheries/resources/blob/master/img/eDNA_meta_COI_primerset.png?raw=true)

Workflow done on HPC. Scripts to run: 

1. 01-fastqc.sh   
2. 02-multiqc.sh  
3. 03-metadata.R
4. 04-ampliseq.sh
5. 05-taxonomicID.sh  

Taxonomic identification is completed using 1) NCBI at 90% percent identity + Least Common Ancestor (LCA) and 2) RDP Classifier via DADA2 with MetaZooGene. After filtering and annotation confirmation, the remaining unassigned sequences are 3) assigned to family-level at best from NCBI at 80% identity. 

This document does quality filtering (fastQC + multiQC), database updates, ASV generation (ampliseq DADA2), taxonomic identification with MetaZooGene (ampliseq DADA2), and lastly, taxonomic identification with NCBI (blast). The blast step is completed at 80% threshold and in future R scripts, this will be split into 90%+ and 80%+ at different steps.

Typical project folder structure (once analysis is run) for Fisheries team. Create `scripts`, `results` (for ampliseq output), `fastqc`, `blast`, and `metadata`. E.g., `mkdir results`. 

```
blast  fastqc  metadata  multiqc_raw_data  multiqc_raw.html  raw_data  results  scripts
```

## Basics

Create/edit a slurm script: `nano name.sh`  
View a file without editing: `less name.sh`  
Start an interactive node: `srun --pty bash`  

## Step 1: Confirm conda environment is available 

The conda environment is started within each slurm script, but to activate conda environment outside of the slurm script to update packages or check what is installed:

```
# Activate conda
source /work/gmgi/miniconda3/bin/activate

# Activate fisheries eDNA conda environment 
conda activate fisheries_eDNA

# List all available environments 
conda env list 

# List all packages installed in fisheries_eDNA
conda list 

# Update a package
conda update [package name]

# Update nextflow ampliseq workflow 
nextflow pull nf-core/ampliseq
``` 
 
## Step 2: Assess quality of raw data  

Background information on [FASTQC](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/05_qc_running_fastqc_interactively.html). This should take ~10 seconds per file for a typical full MiSeq run. If <2 seconds, then check output/error files.  

Create slurm script: `nano 01-fastqc.sh`. Copy below script into file and save (Ctrl+X; Y; Enter).

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

## SET PATHS (USER EDITS)
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
- Start slurm array (e.g., with 138 files) = `sbatch --array=0-137 01-fastqc.sh`. *Note that this starts with file 0*. 

Count the number of files by navigating to raw data folder and using `ls *.gz | wc`. The first column is the number of rows in the generated list (files). 

Notes:  

- This is going to output *many* error and output files. After job completes, navigate to `cd scripts/output/fastqc_output` and use `cat *output.* > ../fastqc_output.txt` to create one file with all the output and `cat *error.* > ../fastqc_error.txt` to create one file with all of the error message outputs. 
- Within the `out_dir` output folder, use `ls *html | wc` to count the number of html output files (1st/2nd column values). This should be equal to the --array range used and the number of raw data files. If not, the script missed some input files so address this before moving on.  


## Step 3: Visualize quality of raw data  

Background information on [MULTIQC](https://multiqc.info/docs/#:~:text=MultiQC%20is%20a%20reporting%20tool%20that%20parses%20results,experiments%20containing%20multiple%20samples%20and%20multiple%20analysis%20steps).

Create slurm script: `nano 02-multiqc.sh`. Copy below script into file and save (Ctrl+X; Y; Enter).

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
## fastqc_output = output from 01-fastqc.sh; fastqc program
fastqc_output="" 
multiqc_dir="" 

## RUN MULTIQC 
multiqc --interactive ${fastqc_output} -o ${multiqc_dir} --filename multiqc_raw.html
```

For typical MiSeq runs, this multiqc step is quick. Use an interactive node instead of sbatch. You'll generate output and error files if you create a .sh script and run with `bash` within an interactive node. Otherwise you could run the code line by line. 

```
## Start an interactive node
# Use srun to claim a node
srun --pty bash 

bash 02-multiqc.sh
```

## Step 4: Downloading and updating reference databases 

### Download and/or update NBCI blast nt database

NCBI is updated daily and therefore needs to be updated each time a project is analyzed. This is the not the most ideal method but we were struggling to get the `-remote` flag to work within slurm because I don't think NU slurm is connected to the internet? NU help desk was helping for awhile but we didn't get anywhere.

Within `/work/gmgi/databases/ncbi`, there is a `update_nt.sh` script with the following code. To run `sbatch update_nt.sh`. This won't take long as it will check for updates rather than re-downloading every time. 

`/work/gmgi/databases/ncbi/update_nt.sh` (this script is already in shared gmgi folder): 

```
#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=update_ncbi_nt
#SBATCH --mem=50G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate fisheries_eDNA

# Create output directory if it doesn't exist
cd /work/gmgi/databases/ncbi/nt

# Update BLAST nt database
update_blastdb.pl --decompress nt

# Print completion message
echo "BLAST nt database update completed"
```

This can take hours so plan ahead - start updating NCBI nt database when you know you are a couple days out from analyzing data. 

View the `update_ncbi_nt.out` file to confirm the echo printed at the end.

### Download Taxonkit 

Download taxonkit to your home directory. **You only need to do this once**, otherwise the program is set-up and you can skip this step. Taxonkit program is in `/work/gmgi/databases/taxonkit` but the program needs required files to be in your home directory.

```
# Move to home directory and download taxonomy information
cd ~
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extract the downloaded file
tar -zxvf taxdump.tar.gz

# Create the TaxonKit data directory 
mkdir -p $HOME/.taxonkit

# Copy the required files to the TaxonKit data directory
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit
```

### Download MetaZooGene

Check the following webpages for the most recent release's version number. If that version number is the same as the database we have on our server, then you don't need to update the database. If not, we need to download the most recent version.

[MetaZooGene North Altantic](https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/html-src/data__MZGdbALL__o02.html)  
[MetaZooGene Global](https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/html-src/data__MZGdbALL__o00.html) 

Instructions to update:

Open RStudio in NU Discovery Cluster OOD and run the following R script to update both the Global and Atlantic versions of MetaZooGene. User needs to update the file name in `base::write()` to reflect the correct version number. This script is very quick, plan for ~3-5 minutes total to run the script, change the version number, and confirm output. 

Path to R script: `/work/gmgi/databases/COI/MetaZooGene/MZG_to_DADA2.R` (this script is already in gmgi shared folder).

```
## GLOBAL DB 

library(magrittr)
input_file = "https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/data-src/MZGdata-coi__MZGdbALL__o00__A.csv.gz"


Global_MZG <- 
  # Download
  readr::read_csv(input_file, 
                              col_names = FALSE, 
                              trim_ws = TRUE, 
                              na = c("", "NA", "N/A"),
                              col_types = cols(.default = "c"),
                              show_col_types = FALSE) %>%
  # Select ID, Sequence and taxa columns
  dplyr::select(X34, X31, X33, X2) 

Edited_Global_MZG <- Global_MZG %>% 
  # Separate taxa column into separate columns
  tidyr::separate(X34, into = as.character(1:21), sep = ";") %>%
  
  # Merge wanted taxa columns and add ">"
  tidyr::unite(col = "Taxa", 1,4,5,7,8,9,10,11,12,16,18,X2, sep = ";" ) %>%
  dplyr::mutate(Taxa = paste(">", Taxa, sep = "")) %>%
  
  # Restructure into one single column vector
  dplyr::select(X33, Taxa, X31) %>% 
  tidyr::pivot_longer(cols = c("Taxa", "X31")) %>%
  dplyr::pull(value)

# Write vector to fasta 
### USER TO CHANGE VERSION NUMBER
Edited_Global_MZG %>% base::write("/work/gmgi/databases/COI/MetaZooGene/DADA2_MZG_v2023-m07-15_Global_modeA.fasta")
  
## ATLANTIC DB

input_file_ATL = "https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/data-src/MZGdata-coi__T4000000__o02__A.csv.gz"

ATL_MZG <- 
  # Download
  readr::read_csv(input_file_ATL, 
                  col_names = FALSE, 
                  trim_ws = TRUE, 
                  na = c("", "NA", "N/A"),
                  col_types = cols(.default = "c"),
                  show_col_types = FALSE) %>%
  # Select ID, Sequence and taxa columns
  dplyr::select(X34, X31, X33, X2) 

Edited_ATL_MZG <- ATL_MZG %>% 
  # Separate taxa column into separate columns
  tidyr::separate(X34, into = as.character(1:21), sep = ";") %>%
  
  # Merge wanted taxa columns and add ">"
  tidyr::unite(col = "Taxa", 1,4,5,7,8,9,10,11,12,16,18,X2, sep = ";" ) %>%
  dplyr::mutate(Taxa = paste(">", Taxa, sep = "")) %>%
  
  # Restructure into one single column vector
  dplyr::select(X33, Taxa, X31) %>% 
  tidyr::pivot_longer(cols = c("Taxa", "X31")) %>%
  dplyr::pull(value)

# Write vector to fasta
### USER TO CHANGE VERSION NUMBER
Edited_ATL_MZG %>% base::write("/work/gmgi/databases/COI/MetaZooGene/DADA2_MZG_v2023-m07-15_NorthAtlantic_modeA.fasta")
```

## Step 4: nf-core/ampliseq 

[Nf-core](https://nf-co.re/): A community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/).  
Nextflow: scalable and reproducible scientific workflows using software containers, used to build wrapper programs like the one we use here.  

[https://nf-co.re/ampliseq/2.11.0]: nfcore/ampliseq is a bioinformatics analysis pipeline used for amplicon sequencing, supporting denoising of any amplicon and supports a variety of taxonomic databases for taxonomic assignment including 16S, ITS, CO1 and 18S. 

![](https://raw.githubusercontent.com/nf-core/ampliseq/2.11.0//docs/images/ampliseq_workflow.png)

We use ampliseq for the following programs:  

- [Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primer. If primer trimming would be omitted, artifactual amplicon sequence variants might be computed by the denoising tool or sequences might be lost due to become labelled as PCR chimera.  
- [DADA2](https://www.nature.com/articles/nmeth.3869) performs fast and accurate sample inference from amplicon data with single-nucleotide resolution. It infers exact amplicon sequence variants (ASVs) from amplicon data with fewer false positives than many other methods while maintaining high sensitivity. This program also assigns taxonomy using a RDP classifier.    

### COI primer sequences (required)

Below is what we used for COI amplicon sequencing. This results in ~313 bp expected ASV. 

LG COI amplicon F: GGWACWGGWTGAACWGTWTAYCCYCC      
LG COI amplicon R: TAIACYTCIGGRTGICCRAARAAYCA       

Ampliseq will automatically calculate the reverse compliment and include this for us.

### Metadata sheet (optional) 

The metadata file has to follow the [QIIME2 specifications](https://docs.qiime2.org/2021.2/tutorials/metadata/). Below is a preview of the sample sheet used for this test. Keep the column headers the same for future use. The first column needs to be "ID" and can only contain numbers, letters, or "-". This is different than the sample sheet. NAs should be empty cells rather than "NA". 

### Create samplesheet sheet for ampliseq 

This file indicates the sample ID and the path to R1 and R2 files. Below is a preview of the sample sheet used in this test. File created on RStudio Interactive on Discovery Cluster using (`create_metadatasheets.R`).  

- sampleID (required): Unique sample IDs, must start with a letter, and can only contain letters, numbers or underscores (no hyphons!).  
- forwardReads (required): Paths to (forward) reads zipped FastQ files  
- reverseReads (optional): Paths to reverse reads zipped FastQ files, required if the data is paired-end  
- run (optional): If the data was produced by multiple sequencing runs, any string  

| sampleID | forwardReads              | reverseReads              | run |
|----------|---------------------------|---------------------------|-----|
| sample1  | ./data/S1_R1_001.fastq.gz | ./data/S1_R2_001.fastq.gz | A   |
| sample2  | ./data/S2_fw.fastq.gz     | ./data/S2_rv.fastq.gz     | A   |
| sample3  | ./S4x.fastq.gz            | ./S4y.fastq.gz            | B   |
| sample4  | ./a.fastq.gz              | ./b.fastq.gz              | B   |

*This is an R script, not slurm script. Open RStudio interactive on Discovery Cluster to run this script.*

Prior to running R script, use the `rawdata` file created for the fastqc slurm array from within the raw data folder to create a list of files. Below is an example from our Offshore Wind project but the specifics of the sampleID will be project dependent. This project had four sequencing runs with different file names. 

On OOD, open RStudio. Create a new R script file:

`03-metadata.R`

```
## Creating ampliseq metadata sheet for GOBLER
## Load libraries 
library(dplyr)
library(stringr)
library(strex) 

### Read in sample sheet 
sample_list <- read.delim2("/work/gmgi/ecosystem-diversity/Gobler/COI/raw_data/rawdata", header=F) %>% 
  dplyr::rename(forwardReads = V1) %>%
  mutate(sampleID = str_after_nth(forwardReads, "data/", 1),

         ### Confirm that _S is correct for your file names, sometimes this is _R
         sampleID = str_before_nth(sampleID, "_S", 1))

         ## add any other mutate functions as needed for sample names
         ## sampleID needs to match a column in your metadata file for later analyses 

# creating sample ID 
sample_list$sampleID <- gsub("-", "_", sample_list$sampleID)

# keeping only rows with R1
sample_list <- filter(sample_list, grepl("R1", forwardReads, ignore.case = TRUE))

# duplicating column 
sample_list$reverseReads <- sample_list$forwardReads

# replacing R1 with R2 in only one column 
sample_list$reverseReads <- gsub("R1", "R2", sample_list$reverseReads)

# rearranging columns 
sample_list <- sample_list[,c(2,1,3)]

sample_list %>% write.csv("/work/gmgi/ecosystem-diversity/Gobler/COI/metadata/samplesheet.csv", 
                          row.names=FALSE, quote = FALSE)
```

### Run nf-core/ampliseq (Cutadapt & DADA2)

Update ampliseq workflow if needed: `nextflow pull nf-core/ampliseq`. 

Create slurm script: `nano 04-ampliseq.sh`. Copy below script into file and save (Ctrl+X; Y; Enter).

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=ampliseq
#SBATCH --mem=70GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for project 
## 2. Adjust SBATCH options above (time, mem, ntasks, etc.) as desired  
## 3. Fill in F primer information based on primer type (no reverse compliment needed)
## 4. Adjust parameters as needed (below is Fisheries team default for COI)

# LOAD MODULES
### This script does not use the eDNA conda environment
module load singularity/3.10.3
module load nextflow/23.10.1

# SET PATHS (PROJECT SPECIFIC)
metadata="/work/gmgi/ecosystem-diversity/Gobler/COI/metadata" 
output_dir="/work/gmgi/ecosystem-diversity/Gobler/COI/results_Global"

# SET PATHS (EDIT VERSION # ON MZG IF NEEDED)
assignTaxonomy="/work/gmgi/databases/COI/MetaZooGene/DADA2_MZG_v2023-m07-15_Global_modeA.fasta"
taxlevels="Kingdom,Phylum,Subphylum,Superclass,Class,Subclass,Infraclass,Superorder,Order,Family,Genus,Species"

nextflow run nf-core/ampliseq -resume \
   -c /work/gmgi/Fisheries/scripts/COI_ampliseq.config \
   -profile singularity \
   --input ${metadata}/samplesheet.csv \
   --FW_primer "GGWACWGGWTGAACWGTWTAYCCYCC" \
   --RV_primer "TAIACYTCIGGRTGICCRAARAAYCA" \
   --outdir ${output_dir} \
   --trunclenf 220 \
   --trunclenr 220 \
   --trunc_qmin 25 \
   --max_ee 2 \
   --min_len_asv 300 \
   --max_len_asv 330 \
   --dada_ref_tax_custom ${assignTaxonomy} \
   --dada_assign_taxlevels ${taxlevels} \
   --skip_dada_addspecies \
   --sample_inference pseudo \
   --ignore_failed_trimming
```

To run:   
- `sbatch 04-ampliseq.sh` 

#### Files generated by ampliseq 

Pipeline summary reports:  

- `summary_report/`
- `summary_report.html`: pipeline summary report as standalone HTML file that can be viewed in your web browser.
- `*.svg*`: plots that were produced for (and are included in) the report.
- `versions.yml`: software versions used to produce this report.

Preprocessing:  

- FastQC: `fastqc/` and `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.  
- Cutadapt: `cutadapt/` and `cutadapt_summary.tsv`: summary of read numbers that pass cutadapt  
- MultiQC: `multiqc`, `multiqc_data/`, `multiqc_plots/` with `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser; 

ASV inferrence with DADA2:  

- `dada2/`, `dada2/args/`, `data2/log/` 
   - `ASV_seqs.fasta`: Fasta file with ASV sequences.
   - `ASV_table.tsv`: Counts for each ASV sequence.
   - `DADA2_stats.tsv`: Tracking read numbers through DADA2 processing steps, for each sample.
   - `DADA2_table.rds`: DADA2 ASV table as R object.
   - `DADA2_table.tsv`: DADA2 ASV table.  
- `dada2/QC/`
   - `*.err.convergence.txt`: Convergence values for DADA2's dada command, should reduce over several magnitudes and approaching 0.  
   - `*.err.pdf`: Estimated error rates for each possible transition. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. The estimated error rates (black line) should be a good fit to the observed rates (points), and the error rates should drop with increased quality.  
   - `*_qual_stats.pdf`: Overall read quality profiles: heat map of the frequency of each quality score at each base position. The mean quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position.  
   - `*_preprocessed_qual_stats.pdf`: Same as above, but after preprocessing.  

We add an ASV length filter that will output `asv_length_filter/` with:  

- `ASV_seqs.len.fasta`: Fasta file with filtered ASV sequences.  
- `ASV_table.len.tsv`: Counts for each filtered ASV sequence.  
- `ASV_len_orig.tsv`: ASV length distribution before filtering.  
- `ASV_len_filt.tsv`: ASV length distribution after filtering.  
- `stats.len.tsv`: Tracking read numbers through filtering, for each sample.  

## Step 5: Running additional NCBI taxonomic assignment ID script 

Ampliseq uses MetaZooGene Global to assign taxonomy but we also use NCBI Eukaryotic database sequences and blastn to find 90+% hits then use Least Common Ancestor approach. After filtering and confirming annotations, we then take NCBI hits at 80+% at family-level only.

Create new folder `blast` in project folder. 

### Blast NCBI at 80% identity 

Create slurm script: `nano 05-blast_NCBI.sh`. Copy below script into file and save (Ctrl+X; Y; Enter).

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=tax_ID
#SBATCH --mem=30GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for project; change db path if not 12S

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate fisheries_eDNA

# SET PATHS (PROJECT SPECIFIC)
## ASV fasta file path, excluding the file name (already included below)
ASV_fasta=""
out=""

# SET PATHS (NO EDITS)
ncbi="/work/gmgi/databases/ncbi/nt"
taxonkit="/work/gmgi/databases/taxonkit"

#### DATABASE QUERY ####
### NCBI database 
blastn -db ${ncbi}/"nt" \
   -query ${ASV_fasta}/ASV_seqs.len.fasta -num_threads 16 \
   -taxids 2759 -culling_limit 50 -max_target_seqs 50 -qcov_hsp_perc 95 -evalue 1e-30 \
   -perc_identity 80 \
   -out ${out}/BLASTResults_NCBI.txt \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

############################

#### TAXONOMIC CLASSIFICATION #### 
## creating list of staxids from all three files 
awk -F $'\t' '{ print $4}' ${out}/BLASTResults_NCBI.txt | sort -u > ${out}/NCBI_sp.txt

## annotating taxid with full taxonomic classification
cat ${out}/NCBI_sp.txt | ${taxonkit}/taxonkit reformat -I 1 -r "Unassigned" > ${out}/NCBI_taxassigned.txt
```

To run:  
- `sbatch 05-blast_NCBI.sh` 


### Reformat taxonkit output for LCA 

Create `05-blast_sort.R` and use RStudio interface on NU OOD to run this script to create a taxid list to use as input for LCA. 

```
## Sorting BLAST output for LCA 
library(tidyverse)

### Load data: user edits project path
NCBI_input="/work/gmgi/ecosystem-diversity/Gobler/COI/blast/BLASTResults_NCBI.txt"
TAXID_output="/work/gmgi/ecosystem-diversity/Gobler/COI/blast/taxids.txt"

blast <- read.table(NCBI_input, header=F,
                    col.names = c("ASV_ID", "sseqid", "sscinames", "staxid", "pident", "length", "mismatch",
                                  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
                    colClasses = c(rep("character", 3), "integer", rep("numeric", 9)))

### Print number of ASVs that had a BLAST output
length(unique(blast$ASV_ID))

### Filter to top pident; this will include ties so any entries with the top pident
blast_filtered <- blast %>% group_by(ASV_ID) %>% slice_max(pident, n=1) 

## nrow should = the number of unique ASVs 
staxid_list <- blast_filtered %>% 
  dplyr::select(ASV_ID, staxid) %>%
  group_by(ASV_ID) %>%
  
  ## combine tax ids into one row and separate with a space 
  summarise(tax_list = paste(staxid, collapse = " ")) 

### USER EDITS PATH
staxid_list %>% write_tsv(TAXID_output)
```

### Run LCA on BLAST output 

Running LCA on Blast output. At this point, the blast hits have been filtered down to the top hit(s). LCA is performed on those top hit(s) rather than everything above 80%. This way, we can take 90% in our first step then 80% later on. 

Create slurm script: `nano 06-LCA.sh`. Copy below script into file and save (Ctrl+X; Y; Enter).

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=LCA
#SBATCH --mem=10GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for project 

# SET PATHS (PROJECT SPECIFIC)
taxid_list="/work/gmgi/ecosystem-diversity/Gobler/COI/blast/taxids.txt"
out="/work/gmgi/ecosystem-diversity/Gobler/COI/blast"

# SET PATHS (NO EDITS)
taxonkit="/work/gmgi/databases/taxonkit"

###################

# Output file for LCA results
lca_results_file="${out}/lca_results.txt"
lca_reformat_file="${out}/lca_results_reformatted.txt"

##### LCA by TaxonKit #####
## Running LCA 
${taxonkit}/taxonkit lca ${taxid_list} -i 2 -o ${lca_results_file}

## Reformatting output 
${taxonkit}/taxonkit reformat ${lca_results_file} -I 3 -o ${lca_reformat_file}

```

### Download files for R analysis

In working folder (for Fisheries team, this is your folder on Box), create the following folders: `data`, `metadata`, and `scripts`. Within `data`, make a folder for `Taxonomic_assignments`. 

Download and place in `data`:  
- `results/summary_report/summary_report.html`    
- `results/overall_summary.tsv`    
- `results/pipeline_info/execution_report_*date*.html`    
- `results/asv_length_filter/ASV_seqs.len.fasta`    
- `results/ASV_table.len.tsv`  

Download and place in `data/Taxonomic_assignments`:  
- `results/phyloseq/dada2_phyloseq.rds`    
- `blast/NCBI_taxassigned.txt`    
- `blast/lca_results_reformatted.txt`    


