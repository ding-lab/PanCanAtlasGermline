# describe workflow to conduct variant QC, concordance calculation and annotation
# by Kuan-lin Huang and Jay Mashl 2017 May @ WashU

variant_QC_annotation.sh which runs script: 
filter_VCF_AF_AD.py

and other dependencies: 
vcfanno
ExAC_config.toml [vcfanno configuration file]
ExAC_nonTCGA frequency file (ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz and ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz.tbi)
vcftools
bed file: gs://dinglab/isb-cgc/tcga/reference_files/all_CDS_ncRNA_ENCODE_multicell_ROI.bed


Note: # this is rough ballpark number check the exact number of this sample H_LS-E2-A10B-10A-01D-A10M-09 in the most recent ppt
Variant QC and filter: 
1. Annotate AF, filter with AF and AD [3 min]
226K -> 102K variants
2. Extract variants in ROI (exon + encode all cell regulatory region), index [28 min]
102K -> 21K variants 


Concordance [may be added in to the last variant QC and filter step before annotation]
run_calc_vcf_concordance.sh which runs:
calc_vcf_concordance.py

kuan-merge-genotype-bigmem which has the genotype VCF already trimmed down to region of interests:
${cancer}.normal.merge.allCDS.vcf.gz
${cancer}.normal.merge.ENCODE.vcf.gz