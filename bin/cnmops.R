#Reading arguments from shell
args = commandArgs(trailingOnly = TRUE)
#FOR DEBUGGING ONLY
#Reading the entry values dor the paths to input file and output file.
#args[1] is the path to the input BAM files
#args[1] <- "/data/test_cnmops/"
#args[2] is the path to the input bed file
#args[2] <- "/reference/039970_D_BED_20120404_pad50.bed"
#args[3] is the path to the output file with the .csv file.
#args[3] <- "results/test_cnmops/cnv_results.csv"

library(cn.mops)
BAMFiles <- list.files(path= args[1], pattern=".bam$", full.names = TRUE)
BEDFile <- read.table(args[2], sep="\t", as.is=TRUE)
genomic_ranges <- GRanges(BEDFile[,1],IRanges(BEDFile[,2],BEDFile[,3]))
Read_counts <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)
resCNMOPS <- exomecn.mops(Read_counts)
IntegerCN <- calcIntegerCopyNumbers(resCNMOPS)
CNVs <- as.data.frame(cnvs(IntegerCN))
write.table(CNVs,file= args[3], quote=FALSE, sep= "\t", row.names=FALSE)
