#!/bin/bash
VCFs=$(gsutil ls gs://dinglab/isb-cgc/tcga/germline/production/merge_v2/3.annotated/whitelisted/by_tcgaBarcode/*vcf.gz)
for file in $VCFs; do 
	gsutil cp $file* .
	vcfName=${file##*/}
	echo "processing ${vcfName}"
	python find_shared_var_relatives.py TCGA_z1_z2_relatives_strict.tsv $vcfName &
	
	NPROC=$(($NPROC+1))
	if [ "$NPROC" -ge 8 ]; then
		wait
		NPROC=0
        fi
done 