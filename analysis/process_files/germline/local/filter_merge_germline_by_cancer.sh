 #!/bin/bash
cancer=$1

    # set limit to more files
    ulimit -n 2500
    echo "Processing cancer type: $line"
    echo "Start time"
    date
    cancer="${line%\\n}"
        mkdir $cancer
        samples=$(gsutil ls gs://dinglab/isb-cgc/tcga/germline/production/${cancer}/)
    for sample in $samples; do 
        sampleNamePre=${sample##*$cancer/}
        sampleName=${sampleNamePre%/}
        echo "Sample "$sample
        echo "Copying vcf for "$sampleName
        gsutil cp ${sample}combine/*gz ${sampleName}.prefilter.snp_indel.vcf.gz
        gsutil cp ${sample}combine/*gz.tbi ${sampleName}.prefilter.snp_indel.vcf.gz.tbi
        gsutil cp ${sample}combine/prefilter.snp_indel.vcf.gz ${cancer}/${sampleName}.prefilter.snp_indel.annotated.ExAC_AF.0.01.AD.3.ROI.vcf.gz
        python filter_VCF_AD.py $sample 5 | bgzip -c > ${cancer}/${sampleName}.AD.5.vcf.gz
        tabix -p vcf ${cancer}/${sampleName}.AD.5.vcf.gz
        # copy the filtered VCF back to storage
        gsutil cp ${cancer}/${sampleName}.AD.5.vcf.gz* ${sample}combine/
    done
    
    # here we may need to limit to the best BAM

    #merge
    bcftools merge --output-type z --output ${cancer}.merge.vcf.gz $(ls -1 ${cancer}/*.vcf.gz | perl -pe 's/\n/ /g')
    #index
    tabix -p vcf ${cancer}.merge.vcf.gz
    #upload
    gsutil cp ${cancer}.merge.vcf.gz* gs://dinglab/isb-cgc/tcga/germline/production/merge
    # delete files
    rm -rf ${cancer}/*.vcf.gz
    rm -rf ${cancer}/*.vcf.gz.tbi
    echo "End time"
    date