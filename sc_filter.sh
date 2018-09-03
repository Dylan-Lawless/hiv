#!/bin/bash

########################################
### This will cover:                ####
### All genes			            ####
### getFuncVar, shared, biallelic   ####
### routine cleanup                 ####
### output an xls to for quick ref  ####
########################################

########################################
####    No Immune list used         ####
####    All exome                   ####
########################################
####    Uncomment to use my list    ####
####    of all known immune genes   ####
####    to norrow the search        ####
####    when a clear immuno         ####
####    driven cause is suspected   ####
########################################

    # location format "X:1-2000", of -b for a bed file or list file with 1 per line.
    # perl /~/vcfhacks-v0.2.0/filterVcfOnLocation.pl \
    # -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vep.vcf \
    # -b ~/txt \
    # -o ~/hiv/hiv_sc/filter/SCcases.filter.vcf

########################################
####    Get variants                ####
########################################

# Always check sample IDs from input source
# in case you move ahead with the wrong source
# ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vep.vcf
# 16 samples: 
# Sample IDs not shown here
    
# get any functional variants
    perl /~/vcfhacks-v0.2.0/getFunctionalVariants.pl \
    -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vep.vcf \
    -s all \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.vcf && \
    
# get functional only when >/=2 have variants in the shared gene
    perl /~/vcfhacks-v0.2.0/getFunctionalVariants.pl \
    -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vep.vcf \
    -f -n 2 \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.SharedGenes.vcf && \
    
# find findBiallelic
# only variants that are common in ALL -s are considered. Use -n to specify the number. Try -n 1.
    perl /~/vcfhacks-v0.2.0/findBiallelic.pl \
    -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vep.vcf \
    -s all \
    -n 1 \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.findBiallelic.all.vcf && \

############################################
####    Rank, annontate and simplify    ####
############################################

############################################
####    Set 1 functional variants       ####
############################################

    perl /~/vcfhacks-v0.2.0/rankOnCaddScore.pl \
    -c /data/shared/cadd/v1.3/*.gz \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.vcf  \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.cadd_ranked.vcf --progress && \
    
    perl /~/vcfhacks-v0.2.0/geneAnnotator.pl \
    -d /~/vcfhacks-v0.2.0/data/geneAnnotatorDb \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.cadd_ranked.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.cadd_ranked.gene_anno && \
    
    perl /~/vcfhacks-v0.2.0/annovcfToSimple.pl \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.cadd_ranked.gene_anno \
    --vep --gene_anno --canonical_only \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.cadd_ranked.gene_anno.simple.xlsx && \
    
    perl /~/vcfhacks-v0.2.0/annovcfToSimple.pl \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.cadd_ranked.gene_anno \
    --vep --gene_anno --canonical_only -u --contains_variant \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.all.cadd_ranked.gene_anno.summarise.simple.xlsx && \

############################################
####    Set 2 functional var shared     ####
############################################

    perl /~/vcfhacks-v0.2.0/rankOnCaddScore.pl \
    -c /data/shared/cadd/v1.3/*.gz \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.SharedGenes.vcf  \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.SharedGenes.cadd_ranked.vcf --progress && \
    
    perl /~/vcfhacks-v0.2.0/geneAnnotator.pl \
    -d /~/vcfhacks-v0.2.0/data/geneAnnotatorDb \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.SharedGenes.cadd_ranked.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.SharedGenes.cadd_ranked.gene_anno && \
    
    perl /~/vcfhacks-v0.2.0/annovcfToSimple.pl \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.SharedGenes.cadd_ranked.gene_anno \
    --vep --gene_anno \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.getFunctionalVariantsVep.SharedGenes.cadd_ranked.gene_anno.simple.xlsx && \

############################################
####    Set 3 biallelic unranked        ####
############################################

    perl /~/vcfhacks-v0.2.0/geneAnnotator.pl \
    -d /~/vcfhacks-v0.2.0/data/geneAnnotatorDb \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.findBiallelic.all.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.findBiallelic.all.gene_anno && \
    
    perl /~/vcfhacks-v0.2.0/annovcfToSimple.pl \
    -i ~/hiv/hiv_sc/filter/SCcases.filter.findBiallelic.all.gene_anno \
    --vep --gene_anno \
    -o ~/hiv/hiv_sc/filter/SCcases.filter.findBiallelic.all.gene_anno.simple.xlsx && \

exit
