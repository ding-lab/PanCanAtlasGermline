#!/bin/bash
# Kuan, Oct 2017 adapted from Jay's script	

for file in PCA_*.merge.vcf.gz; do
    echo '------'
    echo "Reheadering "$file
    chunk=${file%.merge.vcf.gz}
    echo '------'
    gunzip -dc $file | perl replace_vcf_header_sample_with_source_TCGA.pl > $chunk.merge.newheader.txt
    tabix -r $chunk.merge.newheader.txt $file >  $chunk.merge.TCGAbarcode.vcf.gz
    tabix -p vcf $chunk.merge.TCGAbarcode.vcf.gz

    rm -f $chunk.merge.newheader.txt #$file $file.tbi 
done