#!/bin/bash
date

cancer=$1
outFileExon=${cancer}.exon.QCstat.tsv
outFileEncode=${cancer}.encodeROI.QCstat.tsv
touch $outFileExon
touch $outFileEncode

samples=$(gsutil ls gs://dinglab/isb-cgc/tcga/germline/production/${cancer}/)
for sample in $samples; do 
	sampleNamePre=${sample##*$cancer/}
	sampleName=${sampleNamePre%/}
    	echo "Copying vcf for "$sampleName

    	gsutil cp ${sample}combine/*gz ${sampleName}.prefilter.snp_indel.vcf.gz
	gsutil cp ${sample}combine/*gz.tbi ${sampleName}.prefilter.snp_indel.vcf.gz.tbi

	#echo "running" ${cancer}.normal.merge.vcf.gz ${sampleName}.prefilter.snp_indel.vcf.gz $outFile
        python calc_vcf_concordance.py ${sampleName}.prefilter.snp_indel.vcf.gz ${cancer}.normal.merge.allCDS.vcf.gz $outFileExon
        python calc_vcf_concordance.py ${sampleName}.prefilter.snp_indel.vcf.gz ${cancer}.normal.merge.ENCODE.vcf.gz $outFileEncode
	rm -f ${sampleName}.prefilter.snp_indel.vcf.gz*

done
date
