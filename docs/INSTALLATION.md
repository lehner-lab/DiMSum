**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Installation Instructions

## System requirements

DiMSum is expected to work on all Unix-like operating systems.

## Installing DiMSum dependencies

**REQUIRED:** Before [installing the DiMSum package](#installing-dimsum), please ensure that the following required software is installed:

* **[_R_](https://www.r-project.org/) >=v3.6**
* **[_Pandoc_](https://pandoc.org/installing.html) >=v1.17.2**

Pandoc comes bundled with [RStudio](https://rstudio.com/products/rstudio/download/) and the *pandoc* binary can be found in the RStudio *bin/pandoc* directory.

**OPTIONAL:** Additionally, if raw FASTQ files will be processed (with DiMSum *WRAP*), the following software needs to be installed:

* **[_FastQC_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) >=v0.11.3**
* **[_Cutadapt_](https://cutadapt.readthedocs.io/en/stable/) v2.4**
* **[_VSEARCH_](https://github.com/torognes/vsearch) v2.14.2**
* **[_Starcode_](https://github.com/gui11aume/starcode) v1.3**

**TIP:** The easiest way to install these dependencies is with the [Conda](https://docs.conda.io/en/latest/) package management system:
```
conda install -c conda-forge pandoc=1.17.2
conda install -c bioconda fastqc=0.11.9
conda install -c bioconda cutadapt=2.4
conda install -c bioconda vsearch=2.14.2
conda install -c bioconda starcode=1.3
```

**NOTE:** Please ensure that the *$PATH* vairable is set so that these external binaries are available from the command-line prompt. You can add a directory (containing an external binary or symblic link) to your path by adding the following line at the bottom of your *~/.bashrc* file:
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

## Demo DiMSum

In order to test that you have a working installation of the DiMSum pipeline and all necessary software dependencies, run the [Demo](DEMO.md).
