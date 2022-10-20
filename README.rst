About MAXIM 
=============================
MAXIM is a model-based analysis and pipeline of dCas9 Capture-3C-Seq data of multiplexed version. It uses multiplescale Bayesian models for the significance calling of chromatin interactions. It oprovides a versatile and flexible pipeline to analyze the dCas9 Capture-3C-Seq data V2.0 from raw sequencing reads to chromatin loops. MAXIM integrates all steps required for the data analysis, and it supports the multiplexed version that uses batchs of multiple targets (sgRNAs) in one experiments.

Original code can be found in https://github.com/ChenYong-RU/MAXIM, this repo is a fixed version of MAXIM.

1. Prerequisition
-------------------
- Python 3.7 packages (automatically installed)

  - numpy >= 1.13.1
  - scipy >= 0.19.1
  - statsmodels >=0.8.0
  - pandas >= 0.15.2
  - matplotlib >= 2.0.2
  - seaborn >= 0.7.1
  - pysam >= 0.11.2.2

- Other tools

  - bowtie2 >= 2.2.2
  - samtools >= 1.5
  
2. Installation
----------------

::

  > git clone https://github.com/YONGCHENUTD/MAXIM.git
  > cd MAXIM
  > python setup.py install --user

3. Usage of MAXIM
----------------

There are two python script to run MAXIM:
  1. runMAXIM.py
  2. runMAXIM2.py

- Usage of runMAXIM.py

::

  usage: runMAXIM.py -x hg38 -1 sample_R1.fastq.gz [sample_R1.fastq.gz ...] -2
                 sample_R2.fastq.gz [sample_R2.fastq.gz ...] --prefix prefix
                 [--bait chr11:5305934] [--extendsize 100000] [--readlen 36]
                 [--seed 1024] [--smooth-window 100] [--nperm 10000] [-w "."]
                 [-p 10]

- Required parameters:

+--------------------------------------+--------------------------------------------------------------+
|-x hg38, --genome hg38                |Bowtie2 built genome.                                         |
+--------------------------------------+--------------------------------------------------------------+
|-1 S1_R1.fastq.gz [S2_R1.fastq.gz ...]|Read 1 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.|
+--------------------------------------+--------------------------------------------------------------+
|-1 S1_R2.fastq.gz [S2_R2.fastq.gz ...]|Read 2 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.|
+--------------------------------------+--------------------------------------------------------------+
|--prefix                              |Prefix of result files.                                       |
+--------------------------------------+--------------------------------------------------------------+


- Optional parameters:

+--------------------------------------+--------------------------------------------------------------+
|--bait chr11:5305934                  |Bait genomic locus. [Default="chr11:5305934"]                 |
+--------------------------------------+--------------------------------------------------------------+
|--extendsize 100000                   |Length to be extended from bait regions. [Defaut=100000].     |
+--------------------------------------+--------------------------------------------------------------+
|--readlen 36                          |Read length. [Default=36]                                     |
+--------------------------------------+--------------------------------------------------------------+
|--seed 1024                           |Seed to generate random values. [Default=1024].               |
+--------------------------------------+--------------------------------------------------------------+
|-smooth-window 101                    |Smooth window for peak size inference. [Default=101].         |
+--------------------------------------+--------------------------------------------------------------+
|--nperm 10000                         |Number of permutatons. [Default=10000].                       |
+--------------------------------------+--------------------------------------------------------------+
|-w "."                                |Working directory. [Default="."].                             |
+--------------------------------------+--------------------------------------------------------------+
|-p 10                                 |Number of processes. [Default=10].                            |
+--------------------------------------+--------------------------------------------------------------+

- Usage of runMAXIM2.py

::

  usage: runMAXIM2.py -f pairs.gz -m MAPPINGDIR --prefix prefix [--bait chr11:5305934] [--extendsize 100000]
                    [--readlen 36] [--seed 1024] [--smooth-window 100] [--peakstart 5305834] [--peakend 5306034]
                    [--nperm 10000] [-w "."] [-p 10]

- Required parameters:

+--------------------------------------+--------------------------------------------------------------+
|-f / --tbffile  pairs.gz              |tbffile of fixed mate pairs from runMAXIM.py                  |
+--------------------------------------+--------------------------------------------------------------+
|-m / --mappingdir MAPPINGDIR          |Directory containing mapping results from runMAXIM.py         |
+--------------------------------------+--------------------------------------------------------------+
|--prefix                              |Prefix of result files.                                       |
+--------------------------------------+--------------------------------------------------------------+


- Optional parameters:

+--------------------------------------+--------------------------------------------------------------+
|--bait chr11:5305934                  |Bait genomic locus. [Default="chr11:5305934"]                 |
+--------------------------------------+--------------------------------------------------------------+
|--extendsize 100000                   |Length to be extended from bait regions. [Defaut=100000].     |
+--------------------------------------+--------------------------------------------------------------+
|--readlen 36                          |Read length. [Default=36]                                     |
+--------------------------------------+--------------------------------------------------------------+
|--seed 1024                           |Seed to generate random values. [Default=1024].               |
+--------------------------------------+--------------------------------------------------------------+
|-smooth-window 101                    |Smooth window for peak size inference. [Default=101].         |
+--------------------------------------+--------------------------------------------------------------+
|--nperm 10000                         |Number of permutatons. [Default=10000].                       |
+--------------------------------------+--------------------------------------------------------------+
|-w "."                                |Working directory. [Default="."].                             |
+--------------------------------------+--------------------------------------------------------------+
|-p 10                                 |Number of processes. [Default=10].                            |
+--------------------------------------+--------------------------------------------------------------+



4. Demo examples:
-----------------

::

  # hg38 genome
  > wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  > samtools faidx hg38.fa.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 \
  chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 \
  chr22 chrX chrY >hg38_basic.fa

  # Build Bowtie2 genome
  > bowtie2-build hg38_basic.fa hg38_basic
  
  # run runmaxim.py on fastq data
  > runMAXIM.py -x hg38_basic -1 K562_dCas9_sg3-HS1_3C-Capture_R1.fastq.gz \
  -2 K562_dCas9_sg3-HS1_3C-Capture_R2.fastq.gz --prefix sg3-HS1 -w sg3-HS1_raw \
  --bait 'chr11:5226276' 


5. Interpretation of the results
----------------------------------
