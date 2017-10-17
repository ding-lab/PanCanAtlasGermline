#!/bin/bash

batch=${1##*pca_table_split/}
echo "Start time"
date
echo "Processing batch "$batch
echo ""

while IFS='' read -r line || [[ -n "$line" ]]; do
sample="${line%\\n}"

        sampleNamePre=${sample##*production/}
        sampleNamePre2=${sampleNamePre%/combine/prefilter.snp_indel.vcf.gz}
        sampleName=${sampleNamePre2#*/}
        echo "Sample "$sample
        echo "Copying vcf for "$sampleName
        gsutil cp ${sample} ${batch}/${sampleName}.prefilter.snp_indel.vcf.gz
        python filter_VCF_AD.py ${batch}/${sampleName}.prefilter.snp_indel.vcf.gz 5 | bgzip -c > ${batch}/${sampleName}.AD.5.vcf.gz
        tabix -p vcf ${batch}/${sampleName}.AD.5.vcf.gz
        # copy the filtered VCF back to storage
        gsutil cp ${batch}/${sampleName}.AD.5.vcf.gz* gs://dinglab/isb-cgc/tcga/germline/release1.0/individualVCF/

done < "$1"

    #merge
    ~/bin/bcftools-1.5/bcftools merge --output-type z --output ${batch}.merge.vcf.gz $(ls -1 ${batch}/*.vcf.gz | perl -pe 's/\n/ /g')
    #index
    tabix -p vcf ${batch}.merge.vcf.gz
    #upload
    #gsutil cp ${batch}.merge.vcf.gz* gs://dinglab/isb-cgc/tcga/germline/production/merge
    # delete files
    #rm -rf ${batch}/*.vcf.gz
    #rm -rf ${batch}/*.vcf.gz.tbi
    echo "End time"
    date

