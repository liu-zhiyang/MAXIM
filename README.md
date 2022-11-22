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



## 4. Tutorial

```
# data preparation, example data were from GSE139109 titled CAPTURE-3C-seq_K562_sgHS1-5_combined
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/016/SRR10312416/SRR10312416_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/016/SRR10312416/SRR10312416_2.fastq.gz

# run MAXIM for the first round
runMAXIM.py -x path/to/bowtie2_index/GRCh38 -1 SRR10312416_1.fastq.gz -2 SRR10312416_2.fastq.gz --prefix K562_sgLCR -w HS1 --bait "chr11:5274942-5276575" --nperm 10000 -p 16 --peakstart 5274942 --peakend 5276575
# after MAXIM successfully end, you will get 3 folder under your work path "HS1".
# 010ReadMapping contains data reqiured by futher running of MAXIM
# 020Plotting
# 030Modeling contains interaction data in longrang format for WashU Epigenome Browser virtualizing

# run MAIM for another bait, this step can be repeated for more baits
runMAXIM2.py -t HS1/010ReadMapping/K562_sgLCR.pairs.gz -m HS1/010ReadMapping/ --prefix K562_sgLCR -w HS2 --bait "chr11:5280070-5281875" --nperm 10000 -p 16 --peakstart 5280070 --peakend 5281875
# note: you should set another uniq work for a new bait, otherwise the output might be overwrited

# prepare your longrange data for virtualizing
cat `find . -name *wu.longrange` | sort -k1,1 -k2,2n | bgzip -c > mergeBaits.longrange.gz
tabix -p bed mergeBaits.longrange.gz

# virtualizing interactions
# please open http://epigenomegateway.wustl.edu/browser/ , select proper genome, click tracks --> local tracks --> choose 'longrange -- long range interaction data in longrange format' --> choose track file "mergeBaits.longrange.gz"&"mergeBaits.longrange.gz.tbi"
```
