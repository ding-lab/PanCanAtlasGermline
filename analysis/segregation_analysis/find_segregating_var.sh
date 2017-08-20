# run segregation python script
	file=$1
	gsutil cp $file* .
	vcfName=${file##*/}
	echo "Filtering ${vcfName} into $outVCF"
	# frequency check, recalculate AC, AN, and AF based on the cohort
	python find_shared_var_relatives.py TCGA_z1_z2_relatives.tsv $vcfName
	ls -klh ${outVCF}*
	
	rm -f $vcfName