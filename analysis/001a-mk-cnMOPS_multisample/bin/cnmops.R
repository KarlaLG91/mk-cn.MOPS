## DESCRIPTION:
## R script to run cn.MOPS algorithm in order to determine CNVs in multiple samples.
## The script uses BAM files to obtain the read counts, and a BED file to establish segment regions.
## The script also generates a tab delimited file with the copy number clases found on each sample.
##
## USAGE:
## This R script is called by the mkfile of this module with the following command:
##	`Rscript --vanilla bin/cnmops.R <Path to BAM files> <BED file> <Output file>`; 
##					where <Path to BAM files> is the path to any directory found in data/, 
##					<BED file> is the variable $TARGET_BED found in the `config.mk` file in this module, 
##					and <Output file> is any line printed by the `/bin/create_targets` script in this module.
## 
## AUTHOR:
##      Karla Lozano (klg1219sh@gmail.com), for Winter Genomics (http://www.wintergenomics.com/) - 2018

# Reading arguments from shell
#
args = commandArgs(trailingOnly = TRUE)
# Loading the package cn.mops
#
library(cn.mops)
# Loading BAM files
#
BAMFiles <- list.files(path= args[1], pattern=".bam$", full.names = TRUE)
# Loading BED file
#
BEDFile <- read.table(args[2], sep="\t", as.is=TRUE)
# Setting genomic ranges stablished by the BED file
#
genomic_ranges <- GRanges(BEDFile[,1],IRanges(BEDFile[,2],BEDFile[,3]))
# Obtaining read counts form BAM file
#
Read_counts <- getSegmentReadCountsFromBAM(BAMFiles,GR=genomic_ranges)
#Running cn.MOPS algorithm on read counts
#
resCNMOPS <- exomecn.mops(Read_counts)
#Calculating integer copy numbers
#
IntegerCN <- calcIntegerCopyNumbers(resCNMOPS)
CNVs <- as.data.frame(cnvs(IntegerCN))
#Generating tab delimited file
#
write.table(CNVs,file= args[3], quote=FALSE, sep= "\t", row.names=FALSE)
