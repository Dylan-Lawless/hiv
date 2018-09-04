# hiv
These scripts contain examples for basic exome analysis. 
Files and scripts used in this project have a shortened name: "sc", which will not be obvious to if you have not used this set of data before.

The first exmple is a two-step method:
1. sc_genotype.sh
2. sc_filter.sh

The input is a vcf of 16 samples that are genotyped.
For an example of how to get to this point, see "automatic_exomes.sh".

First, sc_genotype.sh shows a basic genotype method, similar to what was done by someone else for the input file. 
It then carries on demonstrating how basic filtering is done. 
Common SNPs and indels are tagged, based on relatively reliable databases. Very common variants are removed.

To get to a smaller list of candidate disease-causing (or in this case, possibly protective) rare variants, some filtering methods are used.
This is demonstrated in the second script; sc_filter.sh.

More sophisticated options are possbible when we have enough controls or non-case samples from the same sequencing runs and library preps.
This will be shown in the next commit.

The R script (gene_of_interest_allele_freq.R) contains code used to make several figures. 
These show the frequencies at which variants were found in the candidate gene. They also demonstrate the rate of variants gene-wide using the GnomAD database. 
Lastly, the variants of interest, and most likely to be importnat for us but combining all of the datasets and showing the key samples, possible index swapping, and population frequency.

In one or two places, the variant IDs are replaced with place holders (variant1,2,3 etc.) and will not work on the original source data, until the column headers are swapped to match.
