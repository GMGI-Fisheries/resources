# Introduction to High Performance Computing Clusters 

High Performance Computing (HPC):  

- The ability to carry out large scale computations to solve complex problems, that either need to process a lot of data, or to have a lot of computing power at their disposal. Basically, any computing system that doesn't fit on a desk can be described as HPC.  
- HPC systems are actually networks of processors. The key principle of HPC lies in the possibility to run massively parallel code to benefit from a large acceleration in runtime.   
- Homogeneous machines only have CPUs while the hybrids have both GPUs and CPUs.   

HPC systems can look like: 

![](https://miro.medium.com/max/3072/1*Ush2RzWxWiAKLXvK_98ZEw.jpeg)

Or smaller, lab-size: 

![](https://telecoms.adaptit.tech/wp-content/uploads/2020/10/hpc-illustration_V2-02-scaled-1-2048x1449.jpg)

### System organization

Regardless of size, HPC clusters will have the following organization. A user upon login, will be on a 'login' node. To run tasks, users can either enter a computing node and run processes themselves or submit a script to a HPC scheduler (e.g., slurm). The job scheduler balances workload distribution and ensures nodes are not overloaded or underutilized. 

![](https://hbctraining.github.io/Intro-to-shell-flipped/img/compute_cluster.png)




Resources:  

- [https://avantonder.github.io/Ghana_course/materials/12-intro_hpc/12.1-intro_hpc.html](https://avantonder.github.io/Ghana_course/materials/12-intro_hpc/12.1-intro_hpc.html)  
- 

