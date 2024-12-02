# Creating a conda environment for eDNA projects 

### eDNA general 

Background information on [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html). 

GMGI Fisheries has a conda environment set-up with all the packages needed for this workflow. Code below was used to create this conda environment.

**DO NOT REPEAT** every time user is running this workflow.

```
# Activate conda
source /work/gmgi/miniconda3/bin/activate

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
conda install bioconda::mothur
```

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

### eDNA COI 

I had issues when installing mothur on the eDNA general conda environment. I got the following error: `mothur: error while loading shared libraries: libgsl.so.25: cannot open shared object file: No such file or directory`. I then installed gsl with `conda install -c conda-forge gsl`. This was already installed and is version 2.7 but I think I need 2.5.. Tried `conda install -c conda-forge gsl=2.5`. This downgraded mothur from 1.48.0-h9f4bb92_2 --> 1.36.1-0. gsl 2.7 will be superseded by higher priority channgel as 2.5. This resulted in a similar error so I updated gsl. Mothur wouldn't update.

Tried to start a different conda environment for mothur 

**DO NOT REPEAT** every time user is running this workflow.

```
# Activate conda
source /work/gmgi/miniconda3/bin/activate

# Creating conda 
conda create --name eDNA_COI mothur

# Installing different versions of packages
conda install -c conda-forge gsl=2.5
conda install -c bioconda mafft
```

The gsl version downgraded mothur (1.48 -> 1.44) and vsearch (2.15 -> 2.13) as well.

