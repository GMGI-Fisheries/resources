# Introduction to High Performance Computing Clusters 

High Performance Computing (HPC) cluster/server:  

- The ability to carry out large scale computations to solve complex problems, that either need to process a lot of data, or to have a lot of computing power at their disposal. Basically, any computing system that doesn't fit on a desk can be described as HPC.  
- HPC systems are actually networks of processors. The key principle of HPC lies in the possibility to run massively parallel code to benefit from a large acceleration in runtime.   
- Homogeneous machines only have CPUs while the hybrids have both GPUs and CPUs.   

HPC systems can look like: 

![](https://miro.medium.com/max/3072/1*Ush2RzWxWiAKLXvK_98ZEw.jpeg)

Or smaller, lab-size: 

![](https://telecoms.adaptit.tech/wp-content/uploads/2020/10/hpc-illustration_V2-02-scaled-1-2048x1449.jpg)

### System organization

A HPC consists of several computers connected in a network. These are called **nodes**:  

- The login nodes are the machines that we connect to and from where we interact with the HPC. These should not be used to run resource-intensive tasks   
- The compute nodes are the high-performance machines on which the actual heavy computations run. Jobs are submitted to the compute nodes through a job scheduler.   

The **job scheduler** is used to submit scripts to be run on the compute nodes:  

- The role of this software is to manage large numbers of jobs being submitted and prioritise them according to their resource needs    
- E.g., slurm  

The filesystem on a HPC is often split between a small (backed) **home directory**, and a large and high-performance (non-backed) **scratch space**:    

- `/home`: The user's home is used for things like configuration files and local software instalation     
- `/scratch`: The scratch space is used for the data and analysis scripts. Usually scrubbed after ~30 days depending on system. Once deleted, cannot be recovered!      
- `/work`: Large, backed working and/or storage space for data    

Regardless of size, HPC clusters will have the following organization. A user upon login, will be on a 'login' node. To run tasks, users can either enter a computing node and run processes themselves or submit a script to a HPC scheduler (e.g., slurm). The job scheduler balances workload distribution and ensures nodes are not overloaded or underutilized. 

![](https://hpc.gwu.edu/files/2021/01/ClusterExample-965x1024.png)


### Resources types 

**CPU (central processing units)** is the "brain" of the computer, performing a wide range of operations and calculations. CPUs can have several "cores", which means they can run tasks in parallel, increasing the throughput of calculations per second. A typical personal computer may have a CPU with 4-8 cores. A single compute node on the HPC may have 32-48 cores (and often these are faster than the CPU on our computers).

**RAM (random access memory)** is a quick access storage where data is temporarily held while being processed by the CPU. A typical personal computer may have 8-32Gb of RAM. A single compute nodes on a HPC may often have >100Gb RAM. 

**GPUs (graphical processing units)** are similar to CPUs, but are more specialised in the type of operations they can do. While less flexible than CPUs, each GPU can do thousands of calculations in parallel. This makes them extremely well suited for graphical tasks, but also more generally for matrix computations and so are often used in machine learning applications. 


### Helpful links  

[https://avantonder.github.io/Ghana_course/materials/12-intro_hpc/12.1-intro_hpc.html](https://avantonder.github.io/Ghana_course/materials/12-intro_hpc/12.1-intro_hpc.html)  

