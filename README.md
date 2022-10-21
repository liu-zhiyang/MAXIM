About MAXIM
===========

MAXIM is a model-based analysis and pipeline of dCas9 Capture-3C-Seq data of multiplexed version. It uses multiplescale Bayesian models for the significance calling of chromatin interactions. It oprovides a versatile and flexible pipeline to analyze the dCas9 Capture-3C-Seq data V2.0 from raw sequencing reads to chromatin loops. MAXIM integrates all steps required for the data analysis, and it supports the multiplexed version that uses batchs of multiple targets (sgRNAs) in one experiments.

Original code can be found in [ChenYong-RU/MAXIM](https://github.com/ChenYong-RU/MAXIM) or the `original` branch of this repo, this repo is a fixed version of MAXIM which also containg a detailed tutorial.

## 1. Requirements

* python
* numpy
* scipy
* statsmodels
* matplotlib
* seaborn
* pysam
* bowtie2
* samtools


## 2. Installation

```bash
  > conda create -n MAXIM
  > conda activate MAXIM
  > conda install numpy scipy statsmodels matplotlib seaborn pysam bowtie2 samtools
  > git clone https://github.com/liu-zhiyang/MAXIM.git
  > cd MAXIM
  > python setup.py install
```


## 3. Usage

There are two scripts. 

1. *runMAXIM.py* is designed to run whole pipeline including raw reads mapping, reads trimming and remapping, fixing pairs, plotting PETs distribution, calling inter/intra-chromosomal interaction and calculating significance.
2. *runMAXIM2.py* is designed to run partial pipeline skipping raw reads mapping, reads trimming and remapping and fixing pairs, which allows you to anlysis interactions from other baits.

### 3.1 *runMAXIM.py*

```bash
usage: runMAXIM.py -x hg38 -1 sample_R1.fastq.gz [sample_R1.fastq.gz ...] -2 sample_R2.fastq.gz [sample_R2.fastq.gz ...] --prefix prefix
               [--bait chr11:5305934] [--extendsize 100000] [--readlen 36]
               [--seed 1024] [--smooth-window 100] [--peakstart 5305834] [--peakend 5306034] 
               [--nperm 10000] [-w "."] [p 10]
```


* Required parameters:

| parameter                              | description                                                    |
| -------------------------------------- | -------------------------------------------------------------- |
| -x hg38, --genome hg38                 | Bowtie2 built genome.                                          |
| -1 S1_R1.fastq.gz [S2_R1.fastq.gz ...] | Read 1 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed. |
| -2 S1_R2.fastq.gz [S2_R2.fastq.gz ...] | Read 2 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed. |
| --prefix                               | Prefix of result files.                                        |



* Optional parameters:

| parameter            | description                                               |
| -------------------- | --------------------------------------------------------- |
| --bait chr11:5305934 | Bait genomic locus. [Default="chr11:5305934"]             |
| --extendsize 100000  | Length to be extended from bait regions. [Defaut=100000]. |
| --readlen 36         | Read length. [Default=36]                                 |
| --seed 1024          | Seed to generate random values. [Default=1024].           |
| -smooth-window 101   | Smooth window for peak size inference. [Default=101].     |
| --nperm 10000        | Number of permutatons. [Default=10000].                   |
| --peakstart          | Start position of bait peak.[Default=5305834]             |
| --peakend            | End position of bait peak.[Default=5306034]               |
| -w "."               | Working directory. [Default="."].                         |
| -p 10                | Number of processes. [Default=10].                        |



### 3.2 *runMAXIM2.py*

```bash
usage: runMAXIM2.py -f pairs.gz -m MAPPINGDIR --prefix prefix [--bait chr11:5305934] [--extendsize 100000]
                  [--readlen 36] [--seed 1024] [--smooth-window 100] [--peakstart 5305834] [--peakend 5306034]
                  [--nperm 10000] [-w "."] [-p 10]
```


* Required parameters:

| parameter                    | description                                           |
| ---------------------------- | ----------------------------------------------------- |
| -f / --tbffile pairs.gz      | tbffile of fixed mate pairs from runMAXIM.py          |
| -m / --mappingdir MAPPINGDIR | Directory containing mapping results from runMAXIM.py |
| --prefix                     | Prefix of result files.                               |


* Optional parameters:

| parameter            | description                                               |
| -------------------- | --------------------------------------------------------- |
| --bait chr11:5305934 | Bait genomic locus. [Default="chr11:5305934"]             |
| --extendsize 100000  | Length to be extended from bait regions. [Defaut=100000]. |
| --readlen 36         | Read length. [Default=36]                                 |
| --seed 1024          | Seed to generate random values. [Default=1024].           |
| -smooth-window 101   | Smooth window for peak size inference. [Default=101].     |
| --nperm 10000        | Number of permutatons. [Default=10000].                   |
| --peakstart          | Start position of bait peak.[Default=5305834]             |
| --peakend            | End position of bait peak.[Default=5306034]               |
| -w "."               | Working directory. [Default="."].                         |
| -p 10                | Number of processes. [Default=10].                        |


**Note:** 

Although `--bait`, `--peakstart` and `--peakend` are optional parameters, I find that they are all needed in the pipeline. So set them with proper value.


## 4. Tutorial