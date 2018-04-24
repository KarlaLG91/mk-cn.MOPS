# mk-cn.MOPS - Pipeline for CNV detection in targeted NGS data (Illumina)

### Abbreviations

-**HTS**: Hight-Throughput Sequencing

-**cn.MOPS**: Copy Number estimation by a Mixture Of PoissonS

-**BAM**: Binary Alignment Map

-**CNVs**: Copy Number Variations

## About mk-alignment

**Objective:**

This module detects CNVs from sorted BAM files using the cn.MOPS algorithms. The module used read depth approach to estimate CNVs in multiple sample simultaneously. 

## Module description

[The cn.MOPS algorithm](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html) works by converting files in BAM format [1] into read count matrices or genomic ranges objects, which are then the input objects.

cn.mops models the depths of coverage across samples at each genomic position. Thus, preventing read count biases across the chromosomes. It applies a mixture of Poisson model that is able to differentiate between variations caused by copy number and variations caused by noise. Since a mixture of Poisson model is created for each genomic position, there is no interference in the results due to read count biases across the entire chromosome. Furthermore, a Bayesian approach is used to determine integer copy numbers based on read variations across samples based on its mixture components. The model set a constant copy number of 2 for all samples and uses the lineal relation between copy numbers and read counts to elucidate any trend that moves further away from the posterior probability expected. The Poisson distribution accounts for background noise and the model assumes that read counts in a segment are similarly distributed across all samples. A CNV is called when there is a deviation across samples  in a number of consecutive segments along a chromosome. [2]

This tool converted each input BAM file into read count matrices for each sample and the targeted regions needed to be specified using a bait file. Consequently, the program partitioned the genome into segments and consider the read counts at a segment to build a model across all samples. 

As a result of this module a .csv file is created recording the levels of copy number clases (cn).

## Pipeline configuration

### Data formats:

* Input data

 Sorted by genome coordinate BAM files 
 
 * Output data
 
 Files in .csv format including the levels of copy number found for each region.

 ````
NOTE:  At least 6 samples are recommended for proper parameter estimation. [2]
````

### Software dependencies:
 
 
 * [mk](https://9fans.github.io/plan9port/man/man1/mk.html "A successor for make.") 
 
 * [cn.MOPS](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html "Copy Number estimation by a Mixture Of PoissonS.") 
 
 
### Configuration file

This pipeline includes a config.mk file (located at mk-cn.MOPS/config.mk, were you can adjust the following parameters:

````

# Path to targeted regions bed file
 TARGET_BED: Must point to an appropriate .bed file.
 
 ````
 
 
### Module parameters

Used by cn.MOPS:

````
args = commandArgs(trailingOnly = TRUE)
library(cn.mops)   ->  Load package
BAMFiles <- list.files(path= args[1], pattern=".bam$", full.names = TRUE)   ->   Get input data from BAM files.
segments <- read.table(args[2], sep="\t", as.is=TRUE)   ->  Get segments from BED file.
gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))  ->  The initial segments in which the reads are counted should be chosen as the regions of the baits, targets or exons established by the BED file.
X <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)  ->  The read count matrix is generated. It requires the genomic coordinates of the predefined segments as GRanges object.
resCNMOPS <- exomecn.mops(X)  ->  run cn.MOPS algorithm.
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)   ->  Calculate integer copy number.
CNVs <- as.data.frame(cnvs(resCNMOPS))  ->  Extract cnvs results.
write.csv(CNVs,file= args[3])  ->  Export cnvs results as .csv.

````

## mk-cn.MOPS directory structure


````


````


## References

\[1\][BAM format](https://genome.sph.umich.edu/wiki/BAM) 
\[2\][cn.MOPS algorithm](https://academic.oup.com/nar/article/40/9/e69/1136601) 




#### Author info
Developed by Karla Lozano (klg1219sh@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/) 2018.