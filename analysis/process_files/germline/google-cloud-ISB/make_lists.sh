#!/bin/bash

# by Jay Mashl, June 2017
# PanCanAtlas Germline

# Get analysisId
LISTS_DIR=analysisID_lists
mkdir -p $LISTS_DIR
for cancerType in  ACC BLCA  BRCA CESC CHOL  COAD DLBC  ESCA  GBM HNSC  KICH   KIRC KIRP  LGG LIHC LUAD   LUSC MESO OV  PAAD PCPG  PRAD  READ SARC  SKCM  STAD TGCT  THCA THYM UCEC  UCS UVM  ; do
    echo $cancerType
    listFile=$LISTS_DIR/$cancerType.ids
    if [ ! -e $listFile ] ; then
	    gsutil ls -d gs://dinglab/isb-cgc/tcga/germline/production/$cancerType/* | sed -e 's/\// /g' | awk '{print $NF}' > $listFile
    fi
done
