#Reading arguments from shell
args = commandArgs(trailingOnly = TRUE)
#FOR DEBUGGING ONLY
#Reading the entry values dor the paths to input file and output file.
#args[1] is the path to the input BAM files
#args[1] <- "/home/consultorwg/mk-cn.MOPS/data/test/"
#args[2] is the path to the input bed file
#args[2] <- "/home/consultorwg/mk-cn.MOPS/data/039970_D_BED_20120404_pad50.bed"
#args[3] is the path to the output file with the .csv file.
#args[3] <- "results/test.csv"

library(cn.mops)
BAMFiles <- list.files(path= args[1], pattern=".bam$", full.names = TRUE)
segments <- read.table(args[2], sep="\t", as.is=TRUE)
gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
X <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)
resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs,file= args[3])
