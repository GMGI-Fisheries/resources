# Metabarcoding workflow for COI amplicon sequencing 

The COI region is commonly used for metabarcoding practices and consequently there are many primer options to choose from. The Fisheries team at GMGI has optimized the Leray Geller set (outlined in red box below). Citation: [Leray et al 2013](https://link.springer.com/article/10.1186/1742-9994-10-34).

We primarily use this set for invertebrate targets and 12S for vertebrate communities. 

![](https://github.com/GMGI-Fisheries/resources/blob/master/img/eDNA_meta_COI_primerset.png?raw=true)

Workflow done on HPC. Scripts to run: 

1. 00-fastqc.sh   
2. 00-multiqc.sh  
3. 01a-metadata.R
4. 01b-ampliseq.sh
5. 02-taxonomicID.sh  

## Step 1: Confirm conda environment is available 

The conda environment is started within each slurm script, but to activate conda environment outside of the slurm script to update packages or check what is installed:

```
# Activate conda
source ~/../../work/gmgi/miniconda3/bin/activate

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
source ~/../../work/gmgi/miniconda3/bin/activate
conda activate fisheries_eDNA

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
source ~/../../work/gmgi/miniconda3/bin/activate
conda activate fisheries_eDNA

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

## Step 4: Downloading and updating reference databases 

### Download and/or update BOLD database 

Visit the [Figshare cite for v4](https://figshare.scilifelab.se/articles/dataset/COI_reference_sequences_from_BOLD_DB/20514192/4) and check for any latest versions. If a new version is available, download the COI references sequences from this webpage: bold_clustered.assignTaxonomy.fasta.gz and bold_clustered.addSpecies.fasta.gz. Via NU Cluster OOD, upload these files to `/work/gmgi/databases/COI/BOLD`. 

I downloaded `taxref_reformat_coidb.sh` from the ampliseq [github page](https://github.com/nf-core/ampliseq/blob/2.8.0/bin/taxref_reformat_coidb.sh). I then changed the first line to be only .gz files (`for f in *.gz; do` instead of `for f in $(ls); do`)

```
# Start on a computing node 
srun --pty bash 

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate fisheries_eDNA

cd /work/gmgi/databases/COI/BOLD 
## upload new .gz files via OOD 

# Edit names to reflect the version #
mv bold_clustered.addspecies.fasta.gz bold_v4_clustered.addspecies.fasta.gz                                                                                           
mv bold_clustered.assigntaxonomy.fasta.gz bold_v4_clustered.assigntaxonomy.fasta.gz 

bash taxref_reformat_coidb.sh
```

These resulting files `addSpecies.fna` and `assignTaxonomy.fna` will be fed into the ampliseq script below.

#### Download and/or update NBCI blast nt database

NCBI is updated daily and therefore needs to be updated each time a project is analyzed. This is the not the most ideal method but we were struggling to get the `-remote` flag to work within slurm because I don't think NU slurm is connected to the internet? NU help desk was helping for awhile but we didn't get anywhere.

Within `/work/gmgi/databases/ncbi`, there is a `update_nt.sh` script with the following code. To run `sbatch update_nt.sh`. This won't take long as it will check for updates rather than re-downloading every time. 

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

View the `update_ncbi_nt.out` file to confirm the echo printed at the end.

*Emma is still troubleshooting the -remote flag to also avoid storing the nt db within our /work/gmgi folder.* 

### Other options 

- [COInr](https://github.com/meglecz/mkCOInr): downloads NCBI and BOLD (Barcode of Life Database) databases to create one database for comparison. COInr is already downloaded in the conda environment and pulls NCBI And BOLD directly.  
- [MARES (MARine Eukaryote Species)](https://www.nature.com/articles/s41597-020-0549-9) and [github](https://github.com/wpearman1996/MARES_database_pipeline): This database combines sequences from both GenBank and BOLD to increase taxonomic coverage and confidence for marine eukaryotes. MARES Github [repo](https://github.com/wpearman1996/MARES_database_pipeline). Paper [link](https://doi.org/10.1038/s41597-020-0549-9).     
- [MIDORI](https://www.reference-midori.info/#:~:text=MIDORI2%20is%20a%20reference%20database%20of%20DNA%20and): A database specifically for COI sequences. MIDORI Reference pulls from GenBank.    

**COInr**

![](https://mkcoinr.readthedocs.io/en/latest/_images/COInr_flowchart_readme.png)

COInr database [instructions](https://mkcoinr.readthedocs.io/en/latest/content/tutorial.html). There are [options](https://mkcoinr.readthedocs.io/en/latest/content/tutorial.html#add-custom-sequences-to-a-database) to include custom sequences if needed.

The latest version of BOLD is 2015 so this 2022 set is the most updated. Use our own NCBI as well to catch recent entries. 

```
cd /work/gmgi/packages 
git clone https://github.com/meglecz/mkCOInr.git

cd /work/gmgi/databases/COI
wget https://zenodo.org/record/6555985/files/COInr_2022_05_06.tar.gz
tar -zxvf COInr_2022_05_06.tar.gz
rm COInr_2022_05_06.tar.gz
mv COInr_2022_05_06 COInr

## converting database information for blast 
perl /work/gmgi/packages/mkCOInr/scripts/format_db.pl -tsv COInr/COInr.tsv -outfmt blast -outdir /work/gmgi/databases/COI/COInr -out COInr_blast

## creating list of sseqID and taxIDs for R df step 
awk '{print $1 "\t" $2}' COInr.tsv > COInr_taxIDlist.tsv
```

**MIDORI**

Visit the [MIDORI website](https://www.reference-midori.info/download.php) to check for the most updated db. This folder is already formatted for blast searching so we don't need to create a blast formatted db. 

```
cd /work/gmgi/databases/COI/MIDORI

## download zip file from MIDORI website for CO1 sequences in BLAST format from nucleotide reference
wget https://www.reference-midori.info/download/Databases/GenBank261_2024-06-15/BLAST/uniq/fasta/MIDORI2_UNIQ_NUC_GB261_CO1_BLAST.fasta.zip
unzip MIDORI2_UNIQ_NUC_GB261_CO1_BLAST.fasta.zip 

## change notation if version is different 
makeblastdb -in MIDORI2_UNIQ_NUC_GB261_CO1_BLAST.fasta -dbtype nucl -out MIDORI2_UNIQ_NUC_GB261_CO1_BLAST.fasta
```

#### Download MetaZooGene

This needs to be done before each run to download most recent version

```
## GLOBAL DB FIRST

library(magrittr)

# input_file = "https://www.st.nmfs.noaa.gov/copepod/collaboration/metazoogene/atlas/data-src/MZGdata-coi__T4000000__o00__A.csv.gz"

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
Edited_ATL_MZG %>% base::write("/work/gmgi/databases/COI/MetaZooGene/DADA2_MZG_v2023-m07-15_NorthAtlantic_modeA.fasta")
```

## Step 4: nf-core/ampliseq 

[Nf-core](https://nf-co.re/): A community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/).  
Nextflow: scalable and reproducible scientific workflows using software containers, used to build wrapper programs like the one we use here.  

[https://nf-co.re/ampliseq/2.11.0]: nfcore/ampliseq is a bioinformatics analysis pipeline used for amplicon sequencing, supporting denoising of any amplicon and supports a variety of taxonomic databases for taxonomic assignment including 16S, ITS, CO1 and 18S. 

![](https://raw.githubusercontent.com/nf-core/ampliseq/2.11.0//docs/images/ampliseq_workflow.png)

We use ampliseq for the following programs:  

- [Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primer. If primer trimming would be omitted, artifactual amplicon sequence variants might be computed by the denoising tool or sequences might be lost due to become labelled as PCR chimera.  
- [DADA2](https://www.nature.com/articles/nmeth.3869) performs fast and accurate sample inference from amplicon data with single-nucleotide resolution. It infers exact amplicon sequence variants (ASVs) from amplicon data with fewer false positives than many other methods while maintaining high sensitivity.   

We skip the taxonomic assignment because we use 3-db approach described in the next section. 

Should we try BOLD through ampliseq?

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

`01a-metadata.R`

```
## Load libraries 

library(dplyr)
library(stringr)
library(strex) 

### Read in sample sheet 

sample_list <- read.delim2("/work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/rawdata", header=F) %>% 
  dplyr::rename(forwardReads = V1) %>%
  mutate(sampleID = str_after_nth(forwardReads, "data/", 1),
         sampleID = str_before_nth(sampleID, "_R", 1))

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

sample_list %>% write.csv("/work/gmgi/Fisheries/eDNA/offshore_wind2023/metadata/samplesheet.csv", 
                          row.names=FALSE, quote = FALSE)
```

### Run nf-core/ampliseq (Cutadapt & DADA2)

Update ampliseq workflow if needed: `nextflow pull nf-core/ampliseq`. 

Testing this on OSW work first. 

`01b-ampliseq.sh`:

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
module load singularity/3.10.3
module load nextflow/23.10.1

# SET PATHS 
metadata="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/metadata" 
output_dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/ampliseq_COIdb_results"
addSpecies="/work/gmgi/databases/COI/BOLD/addSpecies.fna"
assignTaxonomy="/work/gmgi/databases/COI/BOLD/assignTaxonomy.fna"

nextflow run nf-core/ampliseq -resume \
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
   --dada_ref_tax_custom_sp ${addSpecies} \
   --dada2_addspecies_allowmultiple TRUE \
   --sample_inference pseudo \
   --ignore_failed_trimming
```

Could add back in? 
   --max_len 200 \
   --skip_taxonomy 

From zach paper: removing the co-amplified putative nuclear mitochondrial 
pseudogenes (NUMTs) is highly recommended (Creedy et al., 2022; 
Porter & Hajibabaei, 2021; Song et al., 2008). - MetaWorks 
and VTAM implement a step of removing putative NUMTs OR multisample features matrix may be processed with metaMATE 
(Andújar et al., 2021) to remove putative NUMTs and other erroneous sequences (based on, e.g., length and relative read abundance)

DnoisE program - https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04115-6

python package called BOLDigger has been developed to help automate batch query submissions to the BOLD identification engine and can be used to identify COI, ITS, rbcL, and matK sequences (Buchner and Leese, 2020)

In addition to the potential of amino acid translation, the protein coding nature of COI leads to relatively stricter expectations of amplicon length

To run:   
- `sbatch 01b-ampliseq.sh` 

#### Testing MetaZooGene as database

Created folder `ampliseq_MZG_results`
cd to this folder and start script from there `sbatch ../scripts/01-ampliseq_MZG.sh` 

`01-ampliseq_MZG.sh`

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=ampliseq_MZG
#SBATCH --mem=70GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for project 
## 2. Adjust SBATCH options above (time, mem, ntasks, etc.) as desired  
## 3. Fill in F primer information based on primer type (no reverse compliment needed)
## 4. Adjust parameters as needed (below is Fisheries team default for COI)

# LOAD MODULES
module load singularity/3.10.3
module load nextflow/23.10.1

# SET PATHS 
metadata="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/metadata" 
output_dir="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/ampliseq_MZG_results"
assignTaxonomy="/work/gmgi/databases/COI/MetaZooGene/DADA2_MZG_v2023-m07-15_Global_modeA.fasta"
taxlevels="Kingdom, Phlyum, Subphylum, Superclass, Class, Subclass, Infraclass, Superorder, Order, Family, Genus, Species"

nextflow run nf-core/ampliseq -resume \
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

Took out: `--dada_ref_tax_custom_sp ${addSpecies} \` and `--dada2_addspecies_allowmultiple TRUE \` that correspond with `addSpecies="/work/gmgi/databases/COI/BOLD/addSpecies.fna"`

1-22: I tried this with DADA Global MZG fasta first that is DADA assign Taxonomy formatted. I could make Add Species compatible formatting in RStudio on server but we'll see what this output is first.     
1-23: Adding `dada_assign_taxlevels` with 12 entries to match MZG. 

Next try to change bootstrap to 0.8 from 0.5


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


## Step 5: Running additional taxonomic assignment ID script 

`02-taxonomicID.sh`: 

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

# SET PATHS 
ASV_fasta=""
out=""

COInr="/work/gmgi/databases/COI/COInr"
midori="/work/gmgi/databases/COI/MIDORI"
ncbi="/work/gmgi/databases/ncbi/nt"
taxonkit="/work/gmgi/databases/taxonkit"

#### DATABASE QUERY ####
### NCBI database 
blastn -db ${ncbi}/"nt" \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_NCBI.txt \
   -max_target_seqs 20 -perc_identity 99 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## MIDORI database 
blastn -db ${midori}/*.fasta" \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_midori.txt \
   -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid  pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

############################

#### TAXONOMIC CLASSIFICATION #### 
## creating list of staxids from all three files 
awk -F $'\t' '{ print $4}' ${out}/BLASTResults_NCBI.txt | sort -u > ${out}/NCBI_sp.txt

## annotating taxid with full taxonomic classification
cat ${out}/NCBI_sp.txt | ${taxonkit}/taxonkit reformat -I 1 -r "Unassigned" > ${out}/NCBI_taxassigned.txt
```

To run:  
- `sbatch 02-taxonomicID.sh` 