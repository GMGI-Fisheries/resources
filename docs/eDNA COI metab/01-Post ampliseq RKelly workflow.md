# COI metabarcoding : attempting Ryan Kelly + eDNA collaborative team's workflow 

https://github.com/MMARINeDNA/Pipeline/tree/main

Workflow:  
- FastQC + MultiQC  
- Cutadapt  
- DADA2  
- Blast to eukaryotic DNA database  
- LCA 

FastQC through DADA2 in ampliseq. Take this no taxon version and put through blast + LCA. This is similar to our workflow for 12S with DADA2 calling ASVs and then several database comparisons.

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
### NCBI database 97% 
blastn -db ${ncbi}/"nt" \
   -query ${ASV_fasta}/ASV_seqs.len.fasta -num_threads 16 \
   -taxids 2759 -culling_limit 50 -max_target_seqs 50 -qcov_hsp_perc 95 -evalue 1e-30 \
   -perc_identity 97 \
   -out ${out}/BLASTResults_NCBI_97.txt \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

### NCBI database 95%
blastn -db ${ncbi}/"nt" \
   -query ${ASV_fasta}/ASV_seqs.len.fasta -num_threads 16 \
   -taxids 2759 -culling_limit 50 -max_target_seqs 50 -qcov_hsp_perc 95 -evalue 1e-30 \
   -perc_identity 95 \
   -out ${out}/BLASTResults_NCBI_95.txt \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

### NCBI database 93%
blastn -db ${ncbi}/"nt" \
   -query ${ASV_fasta}/ASV_seqs.len.fasta -num_threads 16 \
   -taxids 2759 -culling_limit 50 -max_target_seqs 50 -qcov_hsp_perc 95 -evalue 1e-30 \
   -perc_identity 93 \
   -out ${out}/BLASTResults_NCBI_93.txt \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

### NCBI database 91%
blastn -db ${ncbi}/"nt" \
   -query ${ASV_fasta}/ASV_seqs.len.fasta -num_threads 16 \
   -taxids 2759 -culling_limit 50 -max_target_seqs 50 -qcov_hsp_perc 95 -evalue 1e-30 \
   -perc_identity 91 \
   -out ${out}/BLASTResults_NCBI_91.txt \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'
```

Tried above script:  
- (2/5/2025) I got an error saying num_threads was already defined so I commented this out   
- (2/6/2025) this can't find the ASV len fasta? Took out the -num_threads line that was commented out. The path and file name are correct.  
- (2/6/2025) ASV len path issue fixed but now error that there are too many positional arguments. This line was `-outfmt \"6 sscinames scomnames qseqid sseqid pident length mismatch gapopen qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs\"` but I changed it to `-outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'` to reflect our other scripts.

This worked! 

I then added multiple sections to try out varoius percent IDs.    

97% pident: 542 ASVS (22.81% of reads) assigned    
95% pident: 592 ASVS (22.91% of reads) assigned   
93% pident: 662 ASVS (25.20% of reads) assigned  
91% pident: 776 ASVS (31.69% of reads) assigned  

### R Script for decision tree and taxon ids 

See NU R. Kelly folder. 

### Taxonkit LCA

https://bioinf.shenwei.me/taxonkit/usage/#lca

Use `nano LCA.sh` to create this script and `sbatch LCA.sh` to run it.

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

# SET PATHS
taxid_list="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/RKelly_testing/taxids.txt"
out="/work/gmgi/Fisheries/eDNA/offshore_wind/invertebrate/RKelly_testing"
taxonkit="/work/gmgi/databases/taxonkit"
###################

### No user edits beyond this point 

# Output file for LCA results
lca_results_file="${out}/lca_results.txt"
lca_reformat_file="${out}/lca_results_reformatted.txt"

##### LCA by TaxonKit #####
## Running LCA 
${taxonkit}/taxonkit lca ${taxid_list} -i 2 -o ${lca_results_file}

## Reformatting output 
${taxonkit}/taxonkit reformat ${lca_results_file} -I 3 -o ${lca_reformat_file}
```

This was crazy fast! 

Exporting back to R script to look at assignments. 
