# not in ExAC
VCFs=$(gsutil ls gs://dinglab/isb-cgc/tcga/germline/production/merge_v2/not-in-exac/3.annotated/whitelisted/by_tcgaBarcode/PCA.merge.tcgaBarcode.chr*.anno.whitelist.vcf.gz)
for file in $VCFs; do 
	gsutil cp $file* .
	vcfName=${file##*/}
	date
	echo "Examining vcf-stats for ${vcfName}"
	../plinkseq/build/execs/pseq $vcfName v-stats > ${vcfName}.pseq.vstats.tsv
	../plinkseq/build/execs/pseq $vcfName i-stats > ${vcfName}.pseq.istats.tsv
	#vcf-stats $vcfName > ${vcfName}.vcfstats.json
	
	rm -f $vcfName
	rm -f ${vcfName}.tbi
done

# in ExAC: gs://dinglab/isb-cgc/tcga/germline/production/merge.ExAConly/all/annotated/by_tcgaBarcode/PCA.merge.ExAConly.tcgaBarcode.anno.whitelist.vcf.gz*

file=gs://dinglab/isb-cgc/tcga/germline/production/merge.ExAConly/all/annotated/by_tcgaBarcode/PCA.merge.ExAConly.tcgaBarcode.anno.whitelist.vcf.gz
gsutil cp $file* .
   	vcfName=${file##*/}
        date
	echo "Examining vcf-stats for ${vcfName}"
        ../plinkseq/build/execs/pseq $vcfName v-stats > ${vcfName}.pseq.vstats.tsv
        ../plinkseq/build/execs/pseq $vcfName i-stats > ${vcfName}.pseq.istats.tsv
	#vcf-stats $vcfName > ${vcfName}.vcfstats.json
        
        rm -f $vcfName
        rm -f ${vcfName}.tbi
