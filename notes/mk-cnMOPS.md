#	Module mk-cn.MOPS notes	#

## Module execution:

Given proper configuration (see the README.md for this module), this module can be executed from the directory where mkfile is located, by any of the following commands:

1) `mk <SPECIFIC_TARGET>`; where SPECIFIC_TARGET is any line printed by the `/bin/create_targets` script in the module  .

2) `bin/create_targets | xargs mk`; every possible target printed by bin/create_targets will be generated **in tandem**.

## Expected output example:

````
* Requested output:
results/
└── test_cnmops		
    └── cnv_results.csv	##Tab delimited file with copy number results.

````

## Module observations:

* The module uses read depth approach to estimate CNVs in multiple samples simultaneously using the cn.MOPS algorithm.

* A .csv file is generated recording the levels of copy number clases (CN) for each sample. It gives the _Genomic location (chr, start, end)_ and then four metadata columns. These are _SampleName_, _Median_, _Mean_ and _CN_. 