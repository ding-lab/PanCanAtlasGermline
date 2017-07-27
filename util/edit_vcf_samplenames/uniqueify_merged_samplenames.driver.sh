#!/bin/bash

# Jay Mashl, July 2017

for ct in ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PCA PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM  ; do
    echo '------'
    echo $ct
    echo '------'
    gsutil -m cp gs://dinglab/isb-cgc/tcga/germline/production/merge/$ct.merge.vcf.gz  gs://dinglab/isb-cgc/tcga/germline/production/merge/$ct.merge.vcf.gz.tbi .
    gunzip -dc $ct.merge.vcf.gz | ./replace_vcf_header_sample_with_source.pl   > $ct.merge.newheader.txt
    tabix -r $ct.merge.newheader.txt   $ct.merge.vcf.gz   >  $ct.merge.fixedHeader.vcf.gz
    tabix -p vcf $ct.merge.fixedHeader.vcf.gz

    rm -f $ct.merge.vcf.gz $ct.merge.vcf.gz.tbi

done

    
    
