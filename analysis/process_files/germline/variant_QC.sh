#!/bin/bash
date

vcf=/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.vcf.gz
vcfanno ExAC_config.toml /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.vcf.gz > /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.annotated.vcf

# git clone --recursive https://github.com/vcflib/vcflib.git
# cd vcflib
# make

# /Users/khuang/bin/vcflib/bin/vcffilter /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.annotated.vcf \
# -f "!(AF_Adj > 0.01)"

/Users/khuang/bin/vcflib/bin/vcffilter /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.annotated.vcf \
-f "ExAC_AF_Adj < 0.01" > /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.annotated.rare.vcf


# get genotype file
gsutil cp gs://dinglab/isb-cgc/tcga/genotyping/BRCA/c6f8b1ba-1acd-4ffa-8b02-0d5e6fac0887--TCGA-E2-A10B-10A-01D-A10L-01.vcf* .
~/bin/vcftools_0.1.13/bin/vcftools --gzvcf /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.vcf.gz \
--gzdiff /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/c6f8b1ba-1acd-4ffa-8b02-0d5e6fac0887--TCGA-E2-A10B-10A-01D-A10L-01.vcf.gz \
--diff-site --not-chr X --not-chr MT --out inWES_v_inGENO
# Comparing sites in VCF files...
# Found 14826 sites common to both files.
# Found 206109 sites only in main file.
# Found 486488 sites only in second file.
# Found 13 non-matching overlapping sites.
# After filtering, kept 220948 out of a possible 226043 Sites
# Run Time = 4.00 seconds

# concordant sites
awkt '$4=="B" && $5==$6 && $7==$8' inWES_v_inGENO.diff.sites_in_files > concordant_sites
wc -l concordant_sites
# # discordant
# awkt '$4=="B" && ($5!=$6 || $7!=$8)' inWES_v_inGENO.diff.sites_in_files

~/bin/vcftools_0.1.13/bin/vcftools --gzvcf /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/DVboost/PCA_brca_MC3_example/combine/prefilter.snp_indel.vcf.gz \
--positions concordant_sites --recode --recode-INFO-all --out prefilter.snp_indel_concordant

RVboost.R input_vcf concordant_sites

date