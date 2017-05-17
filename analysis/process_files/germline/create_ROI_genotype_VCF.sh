#!/bin/bash
gcloud compute copy-files --zone us-central1-f ../../../../TCGA_data/reference_files/all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_ForMusic huangkuanlin@kuan-merge-genotype-bigmem:~/
gcloud compute copy-files --zone us-central1-f ../../../../TCGA_data/reference_files/ROI_MultiCell_perid.txt huangkuanlin@kuan-merge-genotype-bigmem:~/

while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "Processing cancer type: $line"
    cancer="${line%\\n}"
    
    # extract exonic region
    vcftools --gzvcf ${cancer}.normal.merge.vcf.gz \
    --bed all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_ForMusic \
    --keep-INFO-all --recode -c | bgzip -c  > ${cancer}.normal.merge.allCDS.vcf.gz &

    # extract encode region 
    vcftools --gzvcf ${cancer}.normal.merge.vcf.gz \
    --bed ROI_MultiCell_perid.txt \
    --keep-INFO-all --recode -c | bgzip -c  > ${cancer}.normal.merge.ENCODE.vcf.gz

    #index
    tabix -p vcf ${cancer}.normal.merge.allCDS.vcf.gz &
    tabix -p vcf ${cancer}.normal.merge.ENCODE.vcf.gz
    
done < cancer_type.txt