# COI metabarcoding : attempting Ryan Kelly + eDNA collaborative team's workflow 

https://github.com/MMARINeDNA/Pipeline/tree/main

Workflow:  
- FastQC + MultiQC  
- Cutadapt  
- DADA2  
- Blast to eukaryotic DNA database  
- LCA 

FastQC through DADA2 in ampliseq. Take this no taxon version and put through blast + LCA. This is similar to our workflow for 12S with DADA2 calling ASVs and then several database comparisons.

Make a eukaryotic blast db from our downloaded nt (if this works, put in update blast script)

```
# Claim a node 
srun --pty bash 

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate fisheries_eDNA
ncbi="/work/gmgi/databases/ncbi/nt"
nt_euk="/work/gmgi/databases/ncbi/nt_euk"

# Filter blast db to Eukaryotic only 
blastn -db ${ncbi}/"nt" -taxids 2759 -out ${nt_euk}/eukaryotic_nt.fasta

# Make custom blast db from this filtered .fasta
makeblastdb -in eukaryotic_nt.fasta -dbtype nucl -out eukaryotic_nt_db
```

### Blast to Eukaryotic DNA

Create a file `nano euk_blast.sh` and paste the below code in that file.

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
ASV_fasta="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/ampliseq_notax_results/asv_length_filter"
out="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/RKelly_testing/blast"

ncbi="/work/gmgi/databases/ncbi/nt"

#### DATABASE QUERY ####
### NCBI database 
blastn -db ${ncbi}/"nt" \
#   -num_threads 16 \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -num_threads 16 -taxids 2759 \
   -culling_limit 50 \
   -max_target_seqs 50 \
   -perc_identity 97 \
   -qcov_hsp_perc 95 \
   -evalue 1e-30 \
   -out ${out}/BLASTResults_NCBI.txt \
   -outfmt \"6 sscinames scomnames qseqid sseqid pident length mismatch gapopen qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs\"
```

Tried above script 2-5-2025  
- I got an error saying num_threads was already defined so I commented this out 

