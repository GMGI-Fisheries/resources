# GMGI in-house Computing Resources 

Gloucester Marine Genomics Institute (GMGI) has 2 in-house servers that are used for bioinformatic analyses. 

![](https://github.com/GMGI-Fisheries/resources/blob/master/img/GMGI_computing_inhouse.png?raw=true)

**Ubuntu Linux operating system aka Humarus** 

Humarus is primarily used for large-scale jobs (e.g., genome assemblies) and thus not the primary working area for Fisheries. 

**Red Hat Enterprise Linux (RHEL) aka Gadus** 

RHEL/Gadus is the primary working area, storage space, and is data is backed up daily to the Synology RackStation in-house. 

## Logging in

Use ssh with username and the correct IP address that can be found on Lab Archives. Follow instructions for entering password. New users will need to get set-up with Jen while onboarding. 

```
ssh username@123.456.7.8
```

## Server Structure 

Once logged in, users are directed to their home directory (`~/`) by default. This space has limited storage and is not intended for regular work. The Fisheries team primarily uses the NU Discovery Cluster for active projects and GMGI's in-house resources for long-term storage and data archiving. Consequently, team members typically use their home directory only for data transfers. 

General server structure: 

Do not edit any folder other than `data`. Only the RHEL main contact is responsible for downloading modules or setting up users. 

```
[estrand@gadus ~]$ cd ../../
[estrand@gadus /]$ ls
bin  boot  data  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
```

Subdirectories within `data`:  
- `prj`: Each lab has their own folder (e.g., `prj/Fisheries`) that is a working area for data and bioinformatic analyses.    
- `resources`: Shared resources like common databases, modules, and scripts live here.  
- `usr` and `var` are for the RHEL main contact only. 

```
[estrand@gadus /]$ cd data/
[estrand@gadus data]$ ls
prj  resources  usr  var
```

Fisheries folders (`prj/Fisheries`):  

We organize these folders by type of analyses or project. I.e., all eDNA projects should be nested within `edna`. 

```
[estrand@gadus Fisheries]$ ls
202402_negatives  edna  epiage  JonahCrabGenome  JonahCrabPopGen  lobster  SandLanceData
```

## Programs and modules 

We run programs as 'modules' that are downloaded by the RHEL's main contact (Jen). If you need a program, send Jen a slack or email with the program name, download link, and if an R package, specify if it is a Bioconductor package or regular CRAN repository package. Once a program is downloaded as a module, this is available for all users. Global installation of programs and R packages helps keep the server uncluttered and not waste space with multiple installations. *Do not install your own copies.* 

Common commands:   
- To find already installed programs: `module avail`      
- To get information about a module: `module help [module/version]` or `module whatis [module/version]`. "help" will provide what the module is, package information including version and install date, and a link to the documentation/github. "whatis" will provide a short, one line description of the program.    
- To load a module: `module load [module/version]` (e.g., `module load bamUtil/v1.0.15`). Loading a module will put all the necessary executables and dependencies for that program in your path so you can call the commands from any location (i.e. your working directory).    

Replace "[module/version]" with the information for your module of interest, as it shows up in "module avail" list.

## Running bioinformatic scripts 

The GMGI RHEL does not currently have a job scheduler program so each user needs to be extremely careful with how much memory and resources their scripts take up. 

Common commands:  
- Check all jobs that are running: `top` and to exit that screen, click Q  
- Check only our user: `top -u username` and to exit that screen, click Q  




