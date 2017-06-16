#!/bin/bash

# by Jay Mashl, June 2017
# PanCanAtlas Germline

#for cancerType in  ACC BLCA  BRCA CESC CHOL  COAD DLBC  ESCA  GBM HNSC  KICH   KIRC KIRP  LGG LIHC LUAD   LUSC MESO OV  PAAD PCPG  PRAD  READ SARC  SKCM  STAD TGCT  THCA THYM UCEC  UCS UVM  ; do
cancerType=$1

analysisList=analysisID_lists/$cancerType.ids

for analysisId in $(cat $analysisList) ; do

    # variable input
    inputPath=gs://dinglab/isb-cgc/tcga/germline/production/${cancerType}/${analysisId}/combine
    VCFIN=prefilter.snp_indel.vcf.gz

    outputPath=gs://dinglab/isb-cgc/tcga/germline/production/${cancerType}/${analysisId}/combine
    VCFOUT=${VCFIN/%vcf.gz/annotated.ExAC_AF.0.01.AD.3.ROI.vcf.gz}

    logsPath=gs://dinglab/isb-cgc/tcga/germline/production/${cancerType}/${analysisId}/combine/annotate_logs

    # fixed input
    EXAC=gs://dinglab/jay/annotation/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
    BEDFILE=gs://dinglab/isb-cgc/tcga/reference_files/all_CDS_ncRNA_ENCODE_multicell_ROI.bed

    ~/dsub/dsub/dsub \
	--project isb-cgc-06-0004 \
	--zones "us-central1-*" \
	--logging $logsPath \
	--input VCFIN=${inputPath}/${VCFIN}  VCFINIDX=${inputPath}/${VCFIN}.tbi  EXAC=${EXAC} EXACIDX=${EXAC}.tbi  BEDFILE=${BEDFILE} \
	--output  VCFOUT=${outputPath}/${VCFOUT}  VCFOUTIDX=${outputPath}/${VCFOUT}.tbi   \
	--command 'cd $(dirname ${VCFIN}) && mv $EXAC $EXACIDX $BEDFILE $(dirname ${VCFIN})  &&  /usr/local/bin/variant_QC_annotation.sh  $(basename ${VCFIN}) &&  mv $(basename ${VCFOUT}) $(basename ${VCFOUTIDX}) $(dirname ${VCFOUT})' \
	--disk-size 20 \
	--min-ram 4 \
	--min-cores 1 \
	--name annotate \
	--image gcr.io/isb-cgc-06-0004/dinglab_pca_analysis:0.1 \
	--scopes  https://www.googleapis.com/auth/compute https://www.googleapis.com/auth/devstorage.full_control https://www.googleapis.com/auth/genomics https://www.googleapis.com/auth/logging.write https://www.googleapis.com/auth/monitoring.write

done

#done
