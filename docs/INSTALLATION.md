**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Installation Instructions

## System requirements

DiMSum is expected to work on all Unix-like operating systems.

## Installing DiMSum using Conda (recommended)

The easiest way to install DiMSum is using the [bioconda package](https://anaconda.org/bioconda/r-dimsum).

Firstly, install the conda package manager (if you don't already have it).

On MacOS, run:
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh
```
On Linux, run:
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

After installing conda you will need to add the bioconda channel as well as the other channels bioconda depends on. Open a new console window/tab and run the following:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Next, optionally, create a dedicated environment for DiMSum and it's dependencies. This is recommended if you already have _R_ and/or _Python_ installations that you need to maintain separately.
```
conda create --name dimsum
conda activate dimsum
```

Finally, install the [DiMSum bioconda package](https://anaconda.org/bioconda/r-dimsum):
```
conda install -c bioconda r-dimsum
```

To check that you have a working installation of DiMSum, run the [Demo](DEMO.md)

## Installing DiMSum dependencies

Installing DiMSum dependencies manually is not recommended. The easiest way to install DiMSum (and its dependencies) is by using the [DiMSum bioconda package](https://anaconda.org/bioconda/r-dimsum). See [Installing DiMSum using Conda](installing-dimsum-using-conda-recommended).

**REQUIRED:** Before [installing DiMSum from GitHub](#installing-dimsum-from-github), please ensure that the following required software is installed:

* **[_R_](https://www.r-project.org/) >=v3.6**
* **[_Pandoc_](https://pandoc.org/installing.html) >=v1.17.2**

Pandoc comes bundled with [RStudio](https://rstudio.com/products/rstudio/download/) and the *pandoc* binary can be found in the RStudio *bin/pandoc* directory.

**OPTIONAL:** Additionally, if raw FASTQ files will be processed (with DiMSum *WRAP*), the following software needs to be installed:

* **[_FastQC_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) >=v0.11.3**
* **[_Cutadapt_](https://cutadapt.readthedocs.io/en/stable/) v2.4**
* **[_VSEARCH_](https://github.com/torognes/vsearch) v2.14.2**
* **[_Starcode_](https://github.com/gui11aume/starcode) v1.3**

**NOTE:** Please ensure that the *$PATH* vairable is set so that these external binaries are available from the command-line prompt. You can add a directory (containing an external binary or symblic link) to your path by adding the following line at the bottom of your *~/.bashrc* or *~/.bash_profile* file:
```
export PATH=EXTERNAL_BINARY_DIRECTORY:$PATH
```

## Installing DiMSum from GitHub

Installing DiMSum from GiHub is not recommended. The easiest way to install DiMSum (and its dependencies) is by using the [DiMSum bioconda package](https://anaconda.org/bioconda/r-dimsum). See [Installing DiMSum using Conda](installing-dimsum-using-conda-recommended).

Before installing DiMSum from GitHub, please ensure that the required [software dependencies](#installing-dimsum-dependencies) are available.

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
Finally, add the cloned DiMSum repository base directory to your path. You can do this by adding the following line at the bottom of your *~/.bashrc* or *~/.bash_profile* file:
```
export PATH=CLONED_DIMSUM_REPOSITORY:$PATH
```

## Demo DiMSum

In order to test that you have a working installation of the DiMSum pipeline and all necessary software dependencies, run the [Demo](DEMO.md).
