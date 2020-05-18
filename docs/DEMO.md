**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Demo

## DiMSum _STEAM_ Demo

In order to run the DiMSum _STEAM_ demo (variant counts to fitness and error), no optional software is required (see [Installation Instructions](docs/INSTALLATION.md) for details).

Simply type the following on the command-line:
```
DiMSum --runDemo T
```
A table of variant counts included in the DiMSum R package is processed and results are saved to a folder named 'DiMSum_Project' in the current working directory. 

## Full DiMSum Demo (_WRAP+STEAM_)

**REQUIRED:** In order to run the full DiMSum demo (FASTQ files to variant fitness and error), _FastQC_, _Cutadapt_, _USEARCH_ and _Starcode_ need to be installed (see [Installation Instructions](docs/INSTALLATION.md) for details). The full demo also requires demo FASTQ files, which can be downloaded from [here](https://www.dropbox.com/s/633skyevl49i0ts/FASTQ.zip?dl=0) and saved locally.

Then, unzip the FASTQ folder and run the full demo from the command-line as follows:
```
unzip FASTQ.zip
DiMSum --runDemo T -i FASTQ
```
Results are saved to a folder named 'DiMSum_Project' in the current working directory.

