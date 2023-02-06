**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Installation Instructions

## System requirements

DiMSum is expected to work on all Unix-like operating systems.

## Installing DiMSum and its dependencies using Conda (recommended)

The easiest way to install DiMSum is by using the [bioconda package](http://bioconda.github.io/recipes/r-dimsum/README.html).

Firstly, install the [Conda](https://docs.conda.io/) package/environment management system (if you don't already have it).

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

**IMPORTANT:** If in doubt, respond with "yes" to the following question during installation: "Do you wish the installer to initialize Miniconda3 by running conda init?". In this case Conda will modify your shell scripts (*~/.bashrc* or *~/.bash_profile*) to initialize Miniconda3 on startup. Ensure that any future modifications to your *$PATH* variable in your shell scripts occur **before** this code to initialize Miniconda3.

After installing Conda you will need to add the bioconda channel as well as the other channels bioconda depends on. Start a new console session (e.g. by closing the current window and opening a new one) and run the following:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Next, optionally, create a dedicated environment for DiMSum and it's dependencies. This is recommended if you already have _R_ and/or _Python_ installations that you would like to maintain in a separate environment.
```
conda create --name dimsum
conda activate dimsum
```
**TIP:** See [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for more information about managing conda environments.

Finally, install the [DiMSum bioconda package](http://bioconda.github.io/recipes/r-dimsum/README.html):
```
conda install -c bioconda r-dimsum
```

To check that you have a working installation of DiMSum, run the [Demo](DEMO.md)

## Installing DiMSum dependencies using Conda

Alternatively, once Conda is installed (or if you already have it) you can install DiMSum dependencies alone by creating a Conda environment from the [dimsum.yaml](../dimsum.yaml) file.

Firstly, clone the DiMSum repository:
```
git clone https://github.com/lehner-lab/DiMSum.git
```
Then use it to create the environment for DiMSum dependencies and activate it:
```
conda env create -f DiMSum/dimsum.yaml
conda activate dimsum
```

Once DiMSum dependencies have been installed successfully you can then [install DiMSum from source](#installing-dimsum-from-source-github) (you can skip the first step to clone the DiMSum repository).

## Installing DiMSum dependencies manually

Installing DiMSum dependencies manually is not recommended. The easiest way to install DiMSum (and its dependencies) is by using the [DiMSum bioconda package](http://bioconda.github.io/recipes/r-dimsum/README.html). See [Installing DiMSum using Conda](#installing-dimsum-using-conda-recommended).

**REQUIRED:** Before [installing DiMSum from source](#installing-dimsum-from-source-github), please ensure that the following required software is installed:

* **[_R_](https://www.r-project.org/) >=v3.6**
* **[_Pandoc_](https://pandoc.org/installing.html) >=v1.17.2**

Pandoc comes bundled with [RStudio](https://rstudio.com/products/rstudio/download/) and the *pandoc* binary can be found in the RStudio *bin/pandoc* directory.

**OPTIONAL:** Additionally, if raw FASTQ files will be processed (with DiMSum *WRAP*), the following software needs to be installed:

* **[_FastQC_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) >=v0.11.3**
* **[_Cutadapt_](https://cutadapt.readthedocs.io/en/stable/) v2.4**
* **[_VSEARCH_](https://github.com/torognes/vsearch) >=v2.17**
* **[_Starcode_](https://github.com/gui11aume/starcode) v1.3**

**NOTE:** Please ensure that the *$PATH* vairable is set so that these external binaries are available from the command-line prompt. You can add a directory (containing an external binary or symblic link) to your path by adding the following line at the bottom of your *~/.bashrc* or *~/.bash_profile* file:
```
export PATH=EXTERNAL_BINARY_DIRECTORY:$PATH
```

## Installing DiMSum from source (GitHub)

Installing DiMSum from source is not recommended. The easiest way to install DiMSum (and its dependencies) is by using the [DiMSum bioconda package](http://bioconda.github.io/recipes/r-dimsum/README.html). See [Installing DiMSum using Conda](#installing-dimsum-using-conda-recommended).

Before installing DiMSum from source, please ensure that the required [software dependencies](#installing-dimsum-dependencies-using-conda) are available.

Firstly, clone the DiMSum repository:
```
git clone https://github.com/lehner-lab/DiMSum.git
```
Then, from the same location run R and enter:
```
if(!require(devtools)) install.packages("devtools")
#devtools::install_deps('DiMSum') #Uncomment this line to install DiMSum dependencies
devtools::install('DiMSum')
```
Finally, add the cloned DiMSum repository base directory to your path. You can do this by adding the following line at the bottom of your *~/.bashrc* or *~/.bash_profile* file:
```
export PATH=CLONED_DIMSUM_REPOSITORY:$PATH
```

## Demo DiMSum

In order to test that you have a working installation of the DiMSum pipeline and all necessary software dependencies, run the [Demo](DEMO.md).
