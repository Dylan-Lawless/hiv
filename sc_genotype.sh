#!/bin/bash

########################################
####    This will cover:            ####
####    Example joint genotype      ####
####    using databases to assess   ####
####    indel snp, dbsnp, evs       ####
####    annotate with vep           ####
########################################

########################################
####    genotyped "SC" samples      ####
########################################

# ~/hiv/hiv_vnp/data/vcfs/SC.cases.vcf
# VNP.cohort.vcf has 17 samples:  VNP001-VNP0017

########################################
####    Typical Joint genotype      ####
####    done with GATK              ####
####    for example                 ####
########################################
# java -Xmx12g -jar ~/GATK/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
# -T GenotypeGVCFs \
# -R database/human_g1k_v37/human_g1k_v37_decoy.fa \
# -D ~/ref/b37/dbSnp146.b37.vcf.gz -stand_call_conf 30 \
# -stand_emit_conf 10 \
# -V ~/gvcf/sample...sort.dedup.indelrealn.recal.HC.g.vcf \
# -V ~/gvcf/sample...sort.dedup.indelrealn.recal.HC.g.vcf \
# -V ~/gvcf/sample...sort.dedup.indelrealn.recal.HC.g.vcf \
# -V ~/gvcf/sample...sort.dedup.indelrealn.recal.HC.g.vcf \
#  ...
#  ...
# -V ~/gvcf/sample...sort.dedup.indelrealn.recal.HC.g.vcf \
# -o ~/hiv/hiv_sc/filter/SCcases.genotype.vcf -nda --showFullBamList -nt 16 && \

############################################
####    Start here with genotyped       ####
####    samples                         ####
############################################

# I normally prefer varaint recal, but this is the first time using 
# these samples. I don't know about the quality etc until I test the original 
# fastq files. Cohort is not large and I don't have info of the library prep / location
# for each sample yet.
# Best to do the basics first; hard filtering and look at output. Expect ~200-300 variants from 16 samples.

# Hard filter
# first, SNPs selected
    java -Xmx12g -jar ~/GATK/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R database/human_g1k_v37/human_g1k_v37_decoy.fa \
    -selectType SNP \
    --variant ~/hiv/hiv_sc/SC.cases.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.raw-snps.vcf && \

# then, INDELs selected
    java -Xmx12g -jar ~/GATK/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R database/human_g1k_v37/human_g1k_v37_decoy.fa \
    --variant ~/hiv/hiv_sc/SC.cases.vcf \
    -selectType INDEL -selectType MNP \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.raw-indels.vcf && \

# filter both
    java -Xmx8g -jar ~/GATK/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R database/human_g1k_v37/human_g1k_v37_decoy.fa \
    -V ~/hiv/hiv_sc/filter/SCcases.genotype.raw-snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "snp_hard_filter" \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.raw-snps.filtered.snvs.vcf && \
    
    java -Xmx8g -jar ~/GATK/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R database/human_g1k_v37/human_g1k_v37_decoy.fa \
    -V ~/hiv/hiv_sc/filter/SCcases.genotype.raw-indels.vcf \
    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filterName "indel_hard_filter" \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.raw-indels.filtered.indels.vcf && \

# combine the outputs
    java -Xmx8g -jar ~/GATK/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
    -T CombineVariants -R database/human_g1k_v37/human_g1k_v37_decoy.fa \
    --variant ~/hiv/hiv_sc/filter/SCcases.genotype.raw-snps.filtered.snvs.vcf \
    --variant ~/hiv/hiv_sc/filter/SCcases.genotype.raw-indels.filtered.indels.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.vcf \
    --genotypemergeoption UNSORTED && \

# filter variants in EdbSNP >/= 1% and not listed as pothogenic by ClinVar
    perl ~/vcfhacks-v0.2.0/annotateSnps.pl \
    -d ~/ref/b37/dbSnp146.b37.vcf.gz ~/ref/b37/clinvar_20160531.vcf.gz -f 1 -pathogenic \
    -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.vcf -t 16 && \

# filter variants in EVS greater >/= 1%
    perl ~/vcfhacks-v0.2.0/filterOnEvsMaf.pl -d ~/ref/evs/ -f 1 --progress \
    -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vcf -t 16 && \

# annotate with vep
    perl ~/variant_effect_predictor/variant_effect_predictor.pl \
    --offline --vcf --everything \
    --dir_cache ~/variant_effect_predictor/vep_cache \
    --dir_plugins ~/variant_effect_predictor/vep_cache/Plugins \
    --plugin Condel,~/variant_effect_predictor/vep_cache/Plugins/config/Condel/config/ \
    --plugin ExAC,~/ref/ExAC/ExAC.r0.3.sites.vep.vcf.gz \
    --plugin SpliceConsensus \
    --fasta ~/variant_effect_predictor/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vcf \
    -o ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vep.vcf \
    --fork 12 && \

exit

# double check
# getSamplesname
# perl ~/vcfhacks-v0.2.0/getSampleNames.pl \
# -i ~/hiv/hiv_sc/filter/SCcases.genotype.fltd-combinedvars.1pcdbsnp.1pcEVS.vep.vcf

