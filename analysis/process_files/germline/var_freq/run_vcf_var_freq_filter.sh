# run vcf frequency filter to rule out (1) alleles with greater than 5% AF in the PCA cohort and without any variant in the final sample list
# update AN/AC fields to the final cohort of 9401 samples

# not in ExAC: gs://dinglab/isb-cgc/tcga/germline/production/merge_v2/not-in-exac/3.annotated/whitelisted/by_tcgaBarcode/PCA.merge.tcgaBarcode.chr*.anno.whitelist.vcf.gz


	file=$1
	gsutil cp $file* .
	vcfName=${file##*/}
	outVCF=${vcfName/.vcf.gz/.cohortAF0.05.vcf.gz}
	echo "Filtering ${vcfName} into $outVCF"
	# frequency check, recalculate AC, AN, and AF based on the cohort
	perl vcf_var_freq_filter.pl --vcf $vcfName | bgzip -c > $outVCF
	tabix -p $outVCF
	ls -klh ${outVCF}*

	gsutil cp ${outVCF}* gs://dinglab/isb-cgc/tcga/germline/production/merge_v2/not-in-exac/3.annotated/whitelisted/by_tcgaBarcode/cohort_AF_filtered
	rm -f $vcfName
	rm -f ${outVCF}*