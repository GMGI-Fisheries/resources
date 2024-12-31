# Introduction to Workflow Managers

Workflow Management Systems (WfMS), such as Snakemake, Galaxy, and Nextflow have been developed specifically to manage computational data-analysis workflows.

Benefits:

- **Run time management**: Management of program execution on the operating system and splitting tasks and data to run at the same time in a process called parallelisation.  
- **Software management**: Use of technology like containers, such as Docker or Singularity, that packages up code and all its dependencies so the application runs reliably from one computing environment to another.  
- **Portability & Interoperability**: Workflows written on one system can be run on another computing infrastructure e.g., local computer, compute cluster, or cloud infrastructure.  
- **Reproducibility**: The use of software management systems and a pipeline specification means that the workflow will produce the same results when re-run, including on different computing platforms.  
- **Reentrancy**: Continuous checkpoints allow workflows to resume from the last successfully executed steps.

### Nextflow 

Basic concepts of Nextflow: [8.1 Nextflow Basic Concepts](https://avantonder.github.io/Ghana_course/materials/08-workflow_managers/8.1-workflow_managers.html#nextflow-basic-concepts).

Nextflow workflows have three main parts:  

1. Processes - describe independent tasks to be run with defined input and output. A process script can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).      
2. Channels - used to manipulate the flow of data from one process to the next.   
3. Workflows - set of instructions for programs involved in the pipeline. A single pipeline can contain many subworkflows. 

![](https://avantonder.github.io/Ghana_course/materials/fig/channel-process_fqc.png)

### nf-core 

[https://nf-co.re/](https://nf-co.re/)

An organization that creates and maintains a curated set of open-source analysis pipelines built using Nextflow. For example, this is a workflow built for rnaseq:

![](https://raw.githubusercontent.com/nf-core/rnaseq/3.8.1/docs/images/nf-core-rnaseq_metro_map_grey.png)

### Snakemake 

[https://snakemake.github.io/](https://snakemake.github.io/)

Similar to Nextflow, the Snakemake workflow management system is a tool for creating reproducible and scalable data analyses. The main difference is that workflows are described via a human readable, Python based language. They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment.

Resources:  

- University of Ghana [documentation](https://avantonder.github.io/Ghana_course/materials/08-workflow_managers/8.1-workflow_managers.html)    
- [Intro to Snakemake](https://ngs-docs.github.io/2023-snakemake-book-draft/intro.html)