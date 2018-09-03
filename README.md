# hiv
These scripts contain examples for basic exome analysis. 
Files and scripts used in this project have a shortened name: "sc", which will not be obvious to if you have not used this set of data before.

The first exmple is a two-step method:
1. sc_genotype.sh
2. sc_filter.sh

The input is a vcf of 16 samples that are genotyped.

First, sc_genotype.sh shows a basic genotype method, similar to what was done by someone else for the input file. 
It then carries on demonstrating how basic filtering is done. 
Common SNPs and indels are tagged, based on relatively reliable databases. Very common variants are removed.

To get to a smaller list of candidate disease-causing (or in this case, possibly protective) rare variants, some filtering methods are used.
This is demonstrated in the second script; sc_filter.sh.

More sophisticated options are possbible when we have enough controls or non-case samples from the same sequencing runs and library preps.
This will be shown in the next commit.
