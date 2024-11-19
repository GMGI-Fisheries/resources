# Creating a conda environment for eDNA projects 

Background information on [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html). 

GMGI Fisheries has a conda environment set-up with all the packages needed for this workflow. Code below was used to create this conda environment.

**DO NOT REPEAT** every time user is running this workflow.

```
# Activate conda
source ~/../../work/gmgi/miniconda3/bin/activate

# Creating conda 
conda create --name fisheries_eDNA

# Installing packages needed for this workflow 
conda install -c bioconda fastqc 
conda install multiqc 
conda install bioconda::nextflow 
conda install conda-forge::singularity
conda install bioconda::blast
conda install nextflow
conda install blast
conda install singularity
conda install -c bioconda vsearch -y
pip install nsdpy
conda install wget
conda install -c bioconda coidb
```

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
