# mk-cn.MOPS - Pipeline for CNV detection in targeted exome NGS data.

### Abbreviations

-**HTS**: High-Throughput Sequencing

-**cn.MOPS**: Copy Number estimation by a Mixture Of PoissonS

-**BAM**: Binary Alignment Map

-**BED**: Browser Extensible Data

-**CNVs**: Copy Number Variations

## About mk-cn.MOPS

**Objective:**

This module detects CNVs from sorted BAM files using the cn.MOPS algorithm. The module uses read depth approach to estimate CNVs in multiple samples simultaneously (001a-mk-cnMOPS_multisample), as well as individual samples (001b-mk-cnMOPS_singlesample).

## Module description

[The cn.MOPS algorithm](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html) works by converting files in BAM format [1] into read count matrices or genomic ranges objects, which are then the input objects. Said genomic ranges are stablished using a BED format file [2].

cn.mops models the depths of coverage across samples at each genomic position. Thus, preventing read count biases across the chromosomes. It applies a mixture of Poisson model that is able to differentiate between variations caused by copy number and variations caused by noise. Since a mixture of Poisson model is created for each genomic position, there is no interference in the results due to read count biases across the entire chromosome. Furthermore, a Bayesian approach is used to determine integer copy numbers based on read variations across samples based on its mixture components. The model set a constant copy number of 2 for all samples and uses the lineal relation between copy numbers and read counts to elucidate any trend that moves further away from the posterior probability expected. The Poisson distribution accounts for background noise and the model assumes that read counts in a segment are similarly distributed across all samples. A CNV is called when there is a deviation across samples  in a number of consecutive segments along a chromosome. [3]

````
IMPORTANT NOTE:

cn.MOPS is reported to not be suitable for single sample analysis. The reason for this is that, the algorithm differentiates variations along samples whether they are copy numbers or noise. The quality of this discrimination increases with the number of samples. [3]
However, the package offers a single sample analysis command `singlecn.mops`. [4]
During the development of this module said command was tested. 
The results obtained turned out to be dissimilar when compared to the results obtained through multi-sample analysis.

````
cn.MOPS installation instructions can be found at: 
[https://bioconductor.org/packages/release/bioc/html/cn.mops.html](https://bioconductor.org/packages/release/bioc/html/cn.mops.html) 

cn.MOPS publication can be found at: [Klambauer G, Schwarzbauer K, Mayr A, Mitterecker A, Clevert D, Bodenhofer U and Hochreiter S (2012). “cn.MOPS: Mixture of Poissons for Discovering Copy Number Variations in Next Generation Sequencing Data with a Low False Discovery Rate.” Nucleic Acids Research, 40, pp. e69. doi: 10.1093/nar/gks003]

## Submodules descriptions

**- 001a-mk-cnMOPS_multisample**

In this module cn.MOPS converted each input BAM file into read count matrices for each sample and specified the targeted regions needed using a BED file. Consequently, the program partitioned the genome into segments and consider the read counts at each segment to build a model across all samples. The algorith uses de linear relation between read counts and copy number to determine CNVs. 

As a result of this module a .csv file is created recording the levels of copy number clases (CN) of all samples. It gives the _Genomic location (chr, start, end)_ and then four metadata columns. These are _SampleName_, _Median_, _Mean_ and _CN_. 
"_CN_" gives the estimated integer copy number of the CNV. The copy number classes default is CN0, CN1, CN2, CN3, .., CN8. CN2 is the normal copy number for diploid samples. CN1 is a heterozygous deletion and CN0 is a homozygous deletion. CN3 through CN8 are amplifications.
"_Median_" and "_Mean_" give the median or mean individual high informative/non-informative call (I/NI call) for this copy number segment. The individual I/NI call is something like the expected log foldchange. Log foldchanges are often used in context with copy number detection.

**- 001b-mk-cnMOPS_singlesample**

In this module cn.MOPS converted each input BAM file into read count matrices for each sample and specified the targeted regions needed using a BED file. Consequently, the program partitioned the genome into segments and consider the read counts at each segment to build a model. The algorith uses de linear relation between read counts and copy number to determine CNVs. 

As a result of this module a .csv file is created recording the levels of copy number clases (CN). It gives the _Genomic location (chr, start, end)_ and then four metadata columns. These are _SampleName_, _Median_, _Mean_ and _CN_. 
"_CN_" gives the estimated integer copy number of the CNV. The copy number classes default is CN0, CN1, CN2, CN3, .., CN8. CN2 is the normal copy number for diploid samples. CN1 is a heterozygous deletion and CN0 is a homozygous deletion. CN3 through CN8 are amplifications.
"_Median_" and "_Mean_" give the median or mean individual high informative/non-informative call (I/NI call) for this copy number segment. The individual I/NI call is something like the expected log foldchange. Log foldchanges are often used in context with copy number detection.

## Pipeline configuration

#### 001a-mk-cnMOPS_multisample

**Data formats:**

* Input data

 -Sorted by genome coordinate BAM files 
 
 -BED file for target baits.
 
 * Output data
 
Files in .csv format including the levels of copy number found for each region for all samples.

 ````
NOTE:  At least 6 samples are recommended for proper parameter estimation. [3]
````

#### 001b-mk-cnMOPS_singlesample

**Data formats:**

* Input data

 -Sorted by genome coordinate BAM files 
 
 -BED file for target baits.
 
 * Output data
 
Files in .csv format including the levels of copy number found for each region for every sample.


### Software dependencies:
 
 
 * [mk](https://9fans.github.io/plan9port/man/man1/mk.html "A successor for make.") 
 
 * [cn.MOPS](https://bioconductor.riken.jp/packages/3.0/bioc/html/cn.mops.html "Copy Number estimation by a Mixture Of PoissonS.") 
 

### Configuration file

This pipeline includes config.mk files (located at every submodule under analysis/SUBMODULE/), where you can adjust several parameters.

Every config.mk file is independent. Please keep it in mind.
 
### Reference files

Under the test_reference/ directory, we provide some files needed for test-runs.

* Genome target regions: BED files where the first column is a chromosome (e.g. "1"), the second and third columns are start and end position of a region, respectively. BED files must end with the suffix ".bed".

## mk-cn.MOPS directory structure


````
mk-cn.MOPS			##Pipeline main directory.
├── analysis			## Directory for submodule organization.
│   ├── 001a-mk-cnMOPS_multisample		##Submodule for CNV detection in multiple samples.
│   └── 001b-mk-cnMOPS_singlesample	####Submodule for CNV detection in single samples.
├── notes			##Notes about proper execution of modules.
│   ├── mk-cnMOPS_multisample.md		##Notes for module execution.
│   └── mk-cnMOPS_singlesample.md		##Notes for module execution.
├── README.md		##This document. General workflow description.
└── test				##Directory for data required for running tests in the module.
    ├── test_data		##Directory for data files used by the module.
    │   └── test_sample	##Directory with sample BAM files for test run.
    │       ├── 06-130529-DH_S1.GATK.realigned.recal.bam	##Example of BAM file.
    │       └── 06-130529-DH_S1.GATK.realigned.recal.bam.bai	##Example of index file of the corresponding BAM file.
    └── test_reference	##Directory for reference files used by the module.
              └── 039970_D_BED_20120404_pad50.bed	##Exampe of BED file.

````


## References

\[1\][BAM format](https://genome.sph.umich.edu/wiki/BAM) 

\[2\][BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) 

\[3\][cn.MOPS algorithm](https://academic.oup.com/nar/article/40/9/e69/1136601) 

\[4\][cn.MOPS manual](https://www.bioconductor.org/packages/3.7/bioc/manuals/cn.mops/man/cn.mops.pdf) 

#### Author info
Developed by Karla Lozano (klg1219sh@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/) 2018.