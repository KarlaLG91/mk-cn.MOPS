## DESCRIPTION:
## mk module for the detection of Copy Number Variations in multiple samples simultaneously, using BAM type of data and a BED file to establish targeted regions. 
## The script uses Copy Number estimation by a Mixture Of PoissonS (cn.MOPS) to detect copy number clases (CN) from a BAM file, defining segments using a BED file.
## 
## USAGE:
## Alternative 1: Single target execution.
##	`mk <SPECIFIC_TARGET>`; where SPECIFIC_TARGET is any line printed by the `/bin/create_targets` script in this module.
##
## Alternative 2: Multiple target tandem execution.
##	`bin/targets | xargs mk` 
##
## AUTHOR:
##	Karla Lozano (klg1219sh@gmail.com), for Winter Genomics (http://www.wintergenomics.com/) - 2018

MKSHELL=/bin/bash

# Load configurations from file
< config.mk

## An example of what happens with the variables set by config.mk
##FOR DEV PURPOSES
config_variables:V: $TARGET_BED
	echo "[DEBUGGING] The BED file is in $prereq"

#Identify copy number clases in .csv file format
#
results/%_cnmops/cnv_results.csv:: data/%_cnmops 
	set -x
	mkdir -p $(dirname $target)
	Rscript --vanilla bin/cnmops.R $prereq $TARGET_BED $target.build \
	&& mv $target.build $target
