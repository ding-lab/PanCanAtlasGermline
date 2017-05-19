#!/bin/bash
# by Kuan-lin Huang 2017 May @ WashU
echo "start time"
date
# requires: vcftools & vcfanno (https://github.com/brentp/vcfanno)

### AF and AD frequency filter ###
vcfFile=$1
bedFile=../../../../TCGA_data/reference_files/all_CDS_ncRNA_ENCODE_multicell_ROI.bed
vcfannoConfigFile=ExAC_config.toml
AF_thres=0.01
AD_thres=3

# annotate with ExAC frequency
annotated_VCF=${vcfFile/vcf.gz/annotated.vcf}
echo "using vcfanno to annotate" ${vcfFile} "into" ${annotated_VCF}
vcfanno ${vcfannoConfigFile} ${vcfFile} > ${annotated_VCF}

filtered_VCF=${annotated_VCF/vcf/ExAC_AF.${AF_thres}.AD.${AD_thres}.vcf}
echo "filtering" ${annotated_VCF} "into" ${filtered_VCF}
python filter_VCF_AF_AD.py ${annotated_VCF} $AF_thres $AD_thres
date

### extract ROI (region of interest) ###
extracted_VCF=${filtered_VCF/vcf/ROI.vcf.gz}
echo "extracting" ${filtered_VCF} "based on" bedFile "into" ${extracted_VCF}
vcftools --vcf $filtered_VCF \
--bed $bedFile \
--keep-INFO-all --recode -c  | bgzip -c  > ${extracted_VCF}

# tabix 
echo "Indexing extracted VCF"
tabix -p vcf ${extracted_VCF}
# remove intermediate VCF files
rm -f $annotated_VCF
rm -f $filtered_VCF
date

# possible option: calculate concordance with genotype file here

# annotate the resulting VCF using VEP

# run CharGer

# move the resulting VCF and files back to storage

echo "end time"
date
