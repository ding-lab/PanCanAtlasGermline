#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "Processing cancer type: $line"
    cancer="${line%\\n}"
    #get file
    gsutil ls gs://dinglab/isb-cgc/tcga/genotyping/tarballs/genotyping.${cancer}.tar
    gsutil cp gs://dinglab/isb-cgc/tcga/genotyping/tarballs/genotyping.${cancer}.tar .
    tar -xvf genotyping.${cancer}.tar 
    
    #merge
    vcf-merge $(ls -1 ${cancer}/*.vcf.gz | grep TCGA-..-....-1.* |  perl -pe 's/\n/ /g') | bgzip -c > ${cancer}.normal.merge.vcf.gz
    #index
    tabix -p vcf ${cancer}.normal.merge.vcf.gz
    #upload
    gsutil cp ${cancer}.normal.merge.vcf.gz* gs://dinglab/isb-cgc/tcga/genotyping/merge
    # delete files
    rm -rf ${cancer}/*
    rm -rf genotyping.${cancer}.tar
done < "$1"

# set unlimit file number higher so it works for breast cancer
ulimit -n 2500