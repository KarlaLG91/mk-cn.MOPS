MKSHELL=/bin/bash

# To load the configurations of a config.mk we use the following line
# This line searches the config.mk in the same directory where the mkfile is.
< config.mk

## An example of what happens with the variables set by config.mk
config_variables:V: $TARGET_BED
	echo "[DEBUGGING] The BED file is in $prereq"

#Generate copy numbers in .csv file format
results/%_cnmops/cnv_results.csv:: data/%_cnmops/ 
	set -x
	mkdir -p $(dirname $target)
	Rscript --vanilla bin/cnmops.R $prereq $TARGET_BED $target.build \
	&& mv $target.build $target

all:V:
	bin/create_targets | xargs mk
