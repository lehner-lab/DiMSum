<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Installation Instructions

## System requirements

DiMSum is expected to work on all Unix-like operating systems.

## Installing DiMSum dependencies

**REQUIRED:** Before [installing the DiMSum package](#installing-dimsum), please ensure that the following required software is installed:

* **[R](https://www.r-project.org/) >=v3.5.2**
* **[Pandoc](https://pandoc.org/installing.html) >=v1.17.2**

Pandoc comes bundled with [RStudio](https://rstudio.com/products/rstudio/download/) and the *pandoc* binary can be found in the RStudio *bin/pandoc* directory.

**OPTIONAL:** Additionally, if raw FASTQ files will be processed (with DiMSum *WRAP*), the following software needs to be installed:

* **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.3**
* **[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) v2.4**
* **[USEARCH 32-bit](https://drive5.com/usearch/download.html) v10.0**
* **[Starcode](https://github.com/gui11aume/starcode) v1.3**

**NOTE:** Please ensure that these external binaries are available on the command-line prompt under their lower-case names i.e. *fastqc*, *cutadapt*, *usearch*, *pandoc* and *starcode*. You can rename a binary from its default by creating a symbolic link if necessary, for example:
```
ln -s usearch10.0.240_i86linux32 usearch 
```
Also please ensure that the *$PATH* vairable is set so that these external binaries are available from the command-line prompt. You can add a directory (containing an external binary or symblic link) to your path by adding the following line at the bottom of your *~/.bashrc* file:
```
export PATH=EXTERNAL_BINARY_DIRECTORY:$PATH
```

## Installing DiMSum

Before installing DiMSum, please ensure that the required [software dependencies](#installing-dimsum-dependencies) are available.

Firstly, clone the DiMSum repository:
```
git clone https://github.com/lehner-lab/DiMSum.git
```
Then, from the same location run R and enter:
```
if(!require(devtools)) install.packages("devtools")
devtools::install_deps('DiMSum')
devtools::install('DiMSum')
```
Finally, add the cloned DiMSum repository base directory to your path. You can do this by adding the following line at the bottom of your *~/.bashrc* file:
```
export PATH=CLONED_DIMSUM_REPOSITORY:$PATH
```
Display the DiMSum usage information by typing the following on the command-line:
```
DiMSum -h
```