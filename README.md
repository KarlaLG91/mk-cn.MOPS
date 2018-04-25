# mk-cn.MOPS - Pipeline for CNV detection in targeted NGS data (Illumina)

### Abbreviations

-**HTS**: Hight-Throughput Sequencing

-**cn.MOPS**: Copy Number estimation by a Mixture Of PoissonS

-**BAM**: Binary Alignment Map

-**BED**: Browser Extensible Data

-**CNVs**: Copy Number Variations

## About mk-alignment

**Objective:**

This module detects CNVs from sorted BAM files using the cn.MOPS algorithm. The module uses read depth approach to estimate CNVs in multiple samples simultaneously. 

## Module description

[The cn.MOPS algorithm](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html) works by converting files in BAM format [1] into read count matrices or genomic ranges objects, which are then the input objects. Said genomic ranges are stablished using a BED format file [2].

cn.mops models the depths of coverage across samples at each genomic position. Thus, preventing read count biases across the chromosomes. It applies a mixture of Poisson model that is able to differentiate between variations caused by copy number and variations caused by noise. Since a mixture of Poisson model is created for each genomic position, there is no interference in the results due to read count biases across the entire chromosome. Furthermore, a Bayesian approach is used to determine integer copy numbers based on read variations across samples based on its mixture components. The model set a constant copy number of 2 for all samples and uses the lineal relation between copy numbers and read counts to elucidate any trend that moves further away from the posterior probability expected. The Poisson distribution accounts for background noise and the model assumes that read counts in a segment are similarly distributed across all samples. A CNV is called when there is a deviation across samples  in a number of consecutive segments along a chromosome. [3]

This tool converted each input BAM file into read count matrices for each sample and the targeted regions needed to be specified using a bait file. Consequently, the program partitioned the genome into segments and consider the read counts at a segment to build a model across all samples. 

As a result of this module a .csv file is created recording the levels of copy number clases (CN). It gives the _Genomic location (chr, start, end)_ and then four metadata columns. These are _SampleName_, _Median_, _Mean_ and _CN_. 
"_CN_" gives the estimated integer copy number of the CNV. The copy number classes default is CN0, CN1, CN2, CN3, .., CN8. CN2 is the normal copy number for diploid samples. CN1 is a heterozygous deletion and CN0 is a homozygous deletion. CN3 through CN8 are amplifications.
"_Median_" and "_Mean_" give the median or mean individual high informative/non-informative call (I/NI call) for this copy number segment. The individual I/NI call is something like the expected log foldchange. Log foldchanges are often used in context with copy number detection.



## Pipeline configuration

### Data formats:

* Input data

 -Sorted by genome coordinate BAM files 
 
 -BED file for target baits.
 
 * Output data
 
Files in .csv format including the levels of copy number found for each region.

 ````
NOTE:  At least 6 samples are recommended for proper parameter estimation. [3]
````

### Software dependencies:
 
 
 * [mk](https://9fans.github.io/plan9port/man/man1/mk.html "A successor for make.") 
 
 * [cn.MOPS](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html "Copy Number estimation by a Mixture Of PoissonS.") 
 
 
### Configuration file

This pipeline includes a config.mk file (located at mk-cn.MOPS/config.mk, were you can adjust the following parameters:

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


## mk-cn.MOPS directory structure


````
mk-cn.MOPS			##Pipeline main directory.
├── bin				##Executables directory.
│   ├── cnmops.R		##Script to run cn.MOPS.
│   └── create_targets	##Script to print every directory required by this module.
├── config.mk			##Configuration file for this module.
├── data				##Directory for sample data
│   └── test_cnmops	##Directory containing a minimum of 6 samples for cn.MOPS analisis. Must follow the following naming condition, "*_cnmops".
│       ├── 01-130529-TM_S1.GATK.realigned.recal.bam	##Example of BAM file.
│       └── 01-130529-TM_S1.GATK.realigned.recal.bam.bai	##Example of index file of the corresponding BAM file.
├── mkfile			##File in mk format, specifying the rules for building every result requested by bin/create_targets.
├── notes			##Notes about proper execution of modules.
│   └── mk-cnMOPS.md	##Notes for module execution.
├── README.md		##This document. General workflow description.
├── reference			##Directory for reference files used by the module
│       └── 039970_D_BED_20120404_pad50.bed	##Exampe of BED file.
└── results			##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.

````


## References

\[1\][BAM format](https://genome.sph.umich.edu/wiki/BAM) 

\[2\][BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) 

\[3\][cn.MOPS algorithm](https://academic.oup.com/nar/article/40/9/e69/1136601) 



#### Author info
Developed by Karla Lozano (klg1219sh@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/) 2018.