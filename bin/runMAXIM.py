#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 20 Oct 2022
# Errors were fixed and the pipeline can run successfully

##########################################################################
#                      Module/Scripts Description
# 
# Copyright (c) 2020 Rowan Unversity
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @version: 2.0.0
# @design: Yong Chen <chenyong@rowan.edu>
# @implementation: Yunfei Wang <yfwang0405@gmail.com>
# @corresponding author:  Yong Chen <chenyong@rowan.edu>
##########################################################################

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import pandas
import argparse


import maxim

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='MAXIM is a model-based analysis and pipeline of dCas9 Capture-3C-Seq data of multiplexed version. It uses multiplescale Bayesian models for the significance calling of chromatin interactions. It oprovides a versatile and flexible pipeline to analyze the dCas9 Capture-3C-Seq data V2.0 from raw sequencing reads to chromatin loops. MAXIM integrates all steps required for the data analysis, and it supports the multiplexed version that uses batchs of multiple targets (sgRNAs) in one experiments.',add_help=False,epilog='dependency numpy, scipy, pandas, pysam, statsmodels')

    pr = p.add_argument_group('Required')
    pr.add_argument("-x","--genome",dest="genome",type=str,metavar="hg38", required=True, help="Bowtie2 built genome.")
    pr.add_argument("-1",dest="fq1",type=str,metavar='sample_R1.fastq.gz',nargs="+",required=True,help="Read 1 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.")
    pr.add_argument("-2",dest="fq2",type=str,metavar='sample_R2.fastq.gz',nargs="+",required=True,help="Read 2 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.")
    pr.add_argument("--prefix",dest="prefix",type=str,metavar='prefix',required=True,help="Prefix of result files.")

    po = p.add_argument_group('Optional')
    po.add_argument("--bait",dest="bait",type=str,metavar="chr11:5305934",default="chr11:5305934",help="Bait genomic locus. [Default=\"chr11:5305934\"]")
    po.add_argument("--extendsize",dest="extendsize",type=int,metavar="100000",default=100000,help="Length to be extended from bait regions. [Defaut=100000]")
    po.add_argument("--readlen",dest="readlen",type=int,metavar="36",default=36,help="Read length. [Default=36]")
    po.add_argument("--seed",dest="seed",type=int,metavar="1024",default=1024,help="Seed to generate random values. [Default=1024].")
    po.add_argument("--smooth-window",dest="smooth_window",type=int,metavar="100",default=100,help="Smooth window for peak size inference. [Default=100].")
    po.add_argument("--peakstart", dest="peakstart", type=int, metavar=5305834, default=0, help="User defined peak start. Used together with '--peakend'.")
    po.add_argument("--peakend", dest="peakend", type=int, metavar=5306034, default=0, help="User defined peak end. Used together with '--peakend'.")
    po.add_argument("--nperm",dest="nperm",type=int,metavar="10000",default=10000,help="Number of permutatons. [Default=10000].")
    po.add_argument("-w",dest="wdir",type=str,metavar='"."',default=".",help="Working directory. [Default=\".\"].")
    po.add_argument("-p",dest='proc',type=int,metavar='10',default=10,help="Number of processes. [Default=10]")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    args = argParser()
    # check parameters
    peaksize = None
    if args.peakstart!=-1 and args.peakend!=-1:
        if args.peakstart==0 and args.peakend==0:
            peaksize = None
        else:
            peaksize = args.peakend - args.peakstart
        maxim.Utils.touchtime("Use user defined peak size.")
    elif args.peakstart!=-1 or args.peakend!=-1:
        maxim.Utils.touchtime("ERROR: both '--peakstart' and '--peakend' should be provided.")

    # Mapping reads to genome
    if "-" in args.bait:
        start, end = args.bait.split(":")[1].split('-')
        args.bait = "{}:{}".format(args.bait.split(":")[0],round((int(start)+int(end))/2))
    mappingdir = args.wdir+"/010ReadMapping"
    fq1, fq2 = ",".join(args.fq1), ",".join(args.fq2)
    mappingdir = maxim.Utils.touchdir(mappingdir)

    # 1st round of mapping
    maxim.Utils.touchtime("FIRST ROUND OF MAPPING ...")
    maxim.Utils.touchtime("MAPPING READ 1 ...")
    maxim.Algorithms.bowtie2_SE(args.genome,fq1,args.prefix+"_R1",proc=args.proc,wdir=mappingdir,min_qual=30)
    maxim.Utils.touchtime("MAPPING READ 2 ...")
    maxim.Algorithms.bowtie2_SE(args.genome,fq2,args.prefix+"_R2",proc=args.proc,wdir=mappingdir,min_qual=30)
    maxim.Utils.touchtime()

    # Split the reads by GATC sites and take the larger one
    maxim.Utils.touchtime("Split read by GATC sites ...")
    maxim.Algorithms.ParseGATCSites("{0}/{1}_R1_un.fastq.gz".format(mappingdir,args.prefix),"{0}/{1}_R1_split.fastq.gz".format(mappingdir,args.prefix))
    maxim.Algorithms.ParseGATCSites("{0}/{1}_R2_un.fastq.gz".format(mappingdir,args.prefix),"{0}/{1}_R2_split.fastq.gz".format(mappingdir,args.prefix))
    maxim.Utils.touchtime()

    # 2nd round of mapping
    maxim.Utils.touchtime("SECOND ROUND OF MAPPING ...")
    maxim.Utils.touchtime("MAPPING READ 1 ...")
    maxim.Algorithms.bowtie2_SE(args.genome,"{0}/{1}_R1_split.fastq.gz".format(mappingdir,args.prefix),args.prefix+"_R1_remap",min_qual=30,proc=args.proc,wdir=mappingdir)
    maxim.Utils.touchtime("MAPPING READ 2 ...")
    maxim.Algorithms.bowtie2_SE(args.genome,"{0}/{1}_R2_split.fastq.gz".format(mappingdir,args.prefix),args.prefix+"_R2_remap",min_qual=30,proc=args.proc,wdir=mappingdir)
    maxim.Utils.touchtime()

    # Fix mate pairs
    maxim.Utils.touchtime("Merge bam files and fix mate pairs ...")
    bams = [mappingdir+args.prefix+f for f in ["_R1.bam", "_R2.bam", "_R1_remap.bam", "_R2_remap.bam"]]
    tbffile = maxim.Algorithms.FixMatePairs(bams,mappingdir+args.prefix,args.proc)
    maxim.Utils.touchtime()

    # Infer peak characteristics from the bait region
    plotdir = args.wdir+"/020Plotting"
    plotdir = maxim.Utils.touchdir(plotdir)
    maxim.Utils.touchtime("Draw bait figures ...")
    tbf = maxim.TabixFile(tbffile,peaksize)
    tbf.setChromSizes(bams[0])    
    tbf.BaitStatsPlot(args.bait,
                      plotdir+args.prefix+"_stats.pdf",
                      extendsize=args.extendsize,
                      readlen=args.readlen,
                      smooth_window=args.smooth_window)
    maxim.Utils.touchtime()

    # Calculate intra- and inter-chrom interactions
    modeldir = args.wdir+"/030Model"
    modeldir = maxim.Utils.touchdir(modeldir)
    maxim.Utils.touchtime("Permutation on intra-chromosomal interactions ...")
    ns, ps = tbf.GetIntraChromLinks(nperm=args.nperm)
    #for n,p in zip(ns,ps):
    #    print n,p
    maxim.Utils.touchtime("Permutation on inter-chromosomal interactions ...")
    n, p = tbf.GetInterChromLinks(nperm=args.nperm)
    
    #print n, p
    maxim.Utils.touchtime()

    # Calculate p values for intra- and inter-chrom interactions.
    maxim.Utils.touchtime("Calculate p values ...")
    tbf.InferBaitPval(modeldir+args.prefix)
    maxim.Utils.touchtime()

    # Ending
    maxim.Utils.touchtime("RunMAXIM finished successfully. Thank you for using MAXIM!")



