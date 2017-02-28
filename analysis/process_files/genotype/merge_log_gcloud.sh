#!/bin/bash

# 2 vCPU; 60GiB <- not enough memory; maybe because vcf-merge open up lots of vcf to do re-indexing?
# 8vCPU; 52GiB mem, 60GiB

# set unlimit file number higher so it works for breast cancer
ulimit -n 2500

nohup bash merge_genotype_by_cancer.sh cancer_type.txt > merge_genotype_by_cancer.log &

# merge
vcf-merge $(ls -1 *.normal.merge.vcf.gz | perl -pe 's/\n/ /g') > all.normal.merge.vcf
# job ended because disk run out of space (60G, zip into partial vcf)
bgzip -c all.normal.merge.vcf > all.normal.merge_partial.vcf.gz
tabix -p vcf all.normal.merge_partial.vcf.gz

# get the remaining SNPs
nohup vcftools --gzvcf BRCA.normal.merge.vcf.gz --gzdiff all.normal.merge_partial.vcf.gz --diff-site --out inBRCA_v_inAllPartial &
awk -F '\t' '$4=="1"{print $1"\t"$2}' inBRCA_v_inAllPartial.diff.sites_in_files > leftover.sites_in_files.positions.txt
# get each cancer types leftover vcf
for file in *.normal.merge.vcf.gz; do
	echo $file
	echo leftover.$file
	vcftools --gzvcf $file --positions leftover.sites_in_files.positions.txt --recode --stdout | bgzip -c > leftover.$file
done
# tabix
for file in leftover*.normal.merge.vcf.gz; do
	tabix -p vcf $file
done

# merge the remaining sites
nohup vcf-merge $(ls -1 leftover*.normal.merge.vcf.gz | perl -pe 's/\n/ /g') > leftover.all.normal.merge.vcf &
bgzip -c leftover.all.normal.merge.vcf > leftover.all.normal.merge.vcf.gz
tabix -p vcf leftover.all.normal.merge.vcf.gz

# the last line is broken; clean it up
zcat all.normal.merge_partial.vcf.gz | head -n -1 | bgzip -c > all.normal.merge_partial_cleaned.vcf.gz
tabix -p vcf all.normal.merge_partial_cleaned.vcf.gz

# merge both and upload
vcf-concat all.normal.merge_partial_cleaned.vcf.gz leftover.all.normal.merge.vcf.gz | bgzip -c > all.normal.merge.vcf.gz &
tabix -p vcf all.normal.merge.vcf.gz
# $ zcat all.normal.merge.vcf.gz | wc -l
# 522763

gsutil cp all.normal.merge.vcf.gz* gs://dinglab/isb-cgc/tcga/genotyping/merge