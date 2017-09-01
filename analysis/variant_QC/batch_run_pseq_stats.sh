# not in ExAC
#VCFs=$(gsutil ls gs://dinglab/isb-cgc/tcga/germline/production/merge_v2/not-in-exac/3.annotated/whitelisted/by_tcgaBarcode/PCA.merge.tcgaBarcode.chr*.anno.whitelist.vcf.gz)
#for file in $VCFs; do
for i in {17..22}; do 
	file=gs://dinglab/isb-cgc/tcga/germline/production/merge_v2/not-in-exac/3.annotated/whitelisted/by_tcgaBarcode/PCA.merge.tcgaBarcode.chr${i}.anno.whitelist.vcf.gz
	bash run_pseq_stats.sh $file &

	NPROC=$(($NPROC+1))
	if [ "$NPROC" -ge 8 ]; then
		wait
		NPROC=0
        fi
done

# in ExAC: gs://dinglab/isb-cgc/tcga/germline/production/merge.ExAConly/all/annotated/by_tcgaBarcode/PCA.merge.ExAConly.tcgaBarcode.anno.whitelist.vcf.gz*

#file=gs://dinglab/isb-cgc/tcga/germline/production/merge.ExAConly/all/annotated/by_tcgaBarcode/PCA.merge.ExAConly.tcgaBarcode.anno.whitelist.vcf.gz
#bash run_pseq_stats.sh $file
