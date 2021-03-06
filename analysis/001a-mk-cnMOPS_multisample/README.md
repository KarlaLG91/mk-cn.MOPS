# 001a-mk-cnMOPS_multisample - Pipeline for CNV detection in targeted exome NGS data for multiple samples simultaneously.

### Abbreviations

-**HTS**: High-Throughput Sequencing

-**cn.MOPS**: Copy Number estimation by a Mixture Of PoissonS

-**BAM**: Binary Alignment Map

-**BED**: Browser Extensible Data

-**CNVs**: Copy Number Variations

## About 001a-mk-cn.MOPS_multisample

**Objective:**

This module detects CNVs from sorted BAM files using the cn.MOPS algorithm. The module uses read depth approach to estimate CNVs in multiple samples simultaneously.

## Module description

In this module cn.MOPS converted each input BAM file into read count matrices for each sample and set segments using a BED file. Then the read counts at each segment were considered to build a model across all samples. The algorith uses de linear relation between read counts and copy number to determine CNVs. 

As a result of this module a .csv file is created recording the levels of copy number clases (CN) of all samples. It gives the _Genomic location (chr, start, end)_ and then four metadata columns. These are _SampleName_, _Median_, _Mean_ and _CN_. 

 ````
NOTE:  At least 6 samples are recommended for proper parameter estimation. [1]
````

## Pipeline configuration

**Data formats:**

* Input data

 -Sorted by genome coordinate BAM files 
 
 -BED file for target baits.
 
 * Output data
 
Files in .csv format including the levels of copy number found for each region for all samples.

### Software dependencies:
 
 
 * [mk](https://9fans.github.io/plan9port/man/man1/mk.html "A successor for make.") 
 
 * [cn.MOPS](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html "Copy Number estimation by a Mixture Of PoissonS.") 
 
 
### Configuration file

This pipeline includes a config.mk file (located at 001a-mk-cnMOPS_multisample/config.mk, were you can adjust the following parameters:

````
# Path to targeted regions bed file
 TARGET_BED: Must point to an appropriate BED file.
 
 ````
 
### Module parameters

Used by cn.MOPS: used by bin/cnmops.R

````
args = commandArgs(trailingOnly = TRUE)
library(cn.mops)   ->  Load package
BAMFiles <- list.files(path= args[1], pattern=".bam$", full.names = TRUE)   ->   Get input data from BAM files.
BEDFile <- read.table(args[2], sep="\t", as.is=TRUE)   ->  Get segments from BED file.
genomic_ranges <- GRanges(BEDFile[,1],IRanges(BEDFile[,2],BEDFile[,3]))  ->  The initial segments in which the reads are counted should be chosen as the regions of the baits, targets or exons established by the BED file.
Read_counts <- getSegmentReadCountsFromBAM(BAMFiles,GR=genomic_ranges) ->  The read count matrix is generated. It requires the genomic coordinates of the predefined segments as GRanges object.
resCNMOPS <- exomecn.mops(Read_counts)  ->  run cn.MOPS algorithm.
IntegerCN <- calcIntegerCopyNumbers(resCNMOPS)   ->  Calculate integer copy number.
CNVs <- as.data.frame(cnvs(IntegerCN))  ->  Extract cnvs results.
write.table(CNVs,file= args[3], quote=FALSE, sep= "\t", row.names=FALSE  ->  Export cnvs results as tab delimited file.

````

## 001a-mk-cn.MOPS_multisample directory structure


````
001a-mk-cnMOPS_multisample	##Submodule main directory.
├── bin				##Executables directory.
│   ├── cnmops.R		##Script to run cn.MOPS.
│   └── create_targets	##Script to print every directory required by this module.
├── config.mk			##Configuration file for this module.
├── data -> ../../test/test_data/	##Symbolic link to data for processing
├── mkfile			##File in mk format, specifying the rules for building every result requested by bin/create_targets.
├── README.md		##This document. General workflow description.
└── results			##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.

````


## References

\[1\][cn.MOPS algorithm](https://academic.oup.com/nar/article/40/9/e69/1136601) 

#### Author info

Developed by Karla Lozano (klg1219sh@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/) 2018.