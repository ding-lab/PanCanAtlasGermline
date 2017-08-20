# not in ExAC
VCFs=$(gsutil ls gs://dinglab/isb-cgc/tcga/germline/production/merge_v2/not-in-exac/3.annotated/whitelisted/by_tcgaBarcode/PCA.merge.tcgaBarcode.chr*.anno.whitelist.vcf.gz)
for file in $VCFs; do 
	bash find_segregating_var.sh $file &

	NPROC=$(($NPROC+1))
	if [ "$NPROC" -ge 8 ]; then
		wait
		NPROC=0
        fi
done

# in ExAC: gs://dinglab/isb-cgc/tcga/germline/production/merge.ExAConly/all/annotated/by_tcgaBarcode/PCA.merge.ExAConly.tcgaBarcode.anno.whitelist.vcf.gz*

file=gs://dinglab/isb-cgc/tcga/germline/production/merge.ExAConly/all/annotated/by_tcgaBarcode/PCA.merge.ExAConly.tcgaBarcode.anno.whitelist.vcf.gz
bash find_segregating_var.sh $file
