# g-LDSC ```V1.0.0```
g-LDSC is a tool for estimating heritability and functional enrichment from GWAS summary statistics. g-LDSC is written under R-4.1.0. 
# Tutorial
g-LDSC is written under R-4.1.0. In this tutorial, we would give a demonstration of how g-LDSC is run under a Linux-based system.
## Installation
Install gldsc (R package) via devtools 
```
#install.packages('devtools')
devtools::install_github("xzw20046/gldsc")
```
## Input file
3 files are required to run g-LDSC:
- GWAS summary statistics
- Pre-calculated (partitioned) LD score matrix
- g-LDSC function file
## GWAS summary statistics
For GWAS summmary statistics, the required formate is shown as follows:
```
SNP A1 A2 N Z
rs1000000 G A 361194 1.5397
rs10000010 T C 361194 -0.850433
rs1000002 C T 361194 0.672368
```
To convert your GWAS result in such format, you could use ```munge_sumstats.py``` from [ldsc](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability).
## Pre-calculated LD score matrix
This file contains the information of LD matrix and annotation information, to get this file, see tutorial of **Calculate LD score matrix** in the futher section.
## g-LDSC function file
The R file ```functions.R``` that contain all g-LDSC functions.
# Usage
## Calculate LD score matrix
To calculate LD score matrix you could use the command shown as follow:
```
Rscript gldsc.run.R \
LDpath=ldblk_ukbb_eur \
annopath=/your file path/baseline \
mafpath=/your file path/1000G_frq \
function=mlfun.R \
out=/your out path/ \
cores=4
```
- ```LDpath``` This flag tells ```g-LDSC``` which LD matrix files to use in calculating LD score matrix. Under this folder all LD matrixs file should be in ```.hdf5``` format. LD matrix of 1000G and UKBB could be download [here](https://github.com/getian107/PRScsx). 
- ```annopath``` & ```mafpath``` This two flag give the location of the annotation and MAF of SNPs. The input format here remain the same as ```.annot``` and ```.M_5_50```in  ```ldsc```. Detail see [here](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability).
- ```out``` This flag tells ```g-LDSC``` where to print the the output.
- ```cores``` This flag tells ```g-LDSC``` how many cores you would like to use for parallel computing.

## Output of LD score matrix calculation
This function will output a file called ```LDSM.pannel.Rdata```. The size of this output file is approximately 10GB.

## Estimate heritability and functional enrichment
You could use the command shown as follow:
```
Rscript gldsc.run.R \
panel=LDSM.pannel.Rdata \
gwas=BMI.sumstats \
function=mlfun.R \
out=/your out path/BMI \
cores=4 
```
## Out put of heritability and functional enrichment estimation
The output of this process will return a data frame with rows represent the result of each functional annotation and columns shown as follow:
```
Taus  Partition_H2  Partition_H2_SD  Enrichment  Enrichment_SD  intercept  intercept_SD  e.stat  P  tau_SD   
```
- ```Taus``` annotation contributor
- ```Taus_SD``` standard error of annotation contributor
- ```Partition_H2``` partitioned SNP-heritability
- ```Partition_H2_SD``` standard error of partitioned SNP-heritability
- ```Enrichment``` fold of enrichment
- ```Enrichment_SD``` jackknife standard error of fold of enrichment
- ```intercept``` confounding bias
- ```intercept_SD``` standard error of confounding bias
- ```e.stat``` t-statistics of function enrichment
- ```P``` P-value of function enrichment (t-test)

# More information
Author: Zewei Xiong (the University of Hong Kong): xzw20046@163.com
