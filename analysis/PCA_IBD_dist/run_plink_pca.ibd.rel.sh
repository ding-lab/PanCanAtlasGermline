# use plink 1.9 to:
# calculate PCA/IBS/relatedness, before and after pruning
plink --vcf all.normal.merge.vcf.gz --pca --out all.normal.merge.vcf.pca
plink --vcf all.normal.merge.vcf.gz --pca --maf 0.15 --out all.normal.merge.vcf.MAF0.15.pca &

plink --vcf all.normal.merge.vcf.gz --distance square --indep 50 5 2 --maf 0.05 --out all.normal.merge.indep_50_5_2.vcf.dist

plink --vcf all.normal.merge.vcf.gz --make-rel square --indep 50 5 2 --maf 0.05 --out all.normal.merge.indep_50_5_2.vcf.rel

rm -f *temporary*

# make tar balls
tar -cvzf plink_dist_rel.tar *indep*
tar -cvzf plink_pca.tar *pca*

# transfer to local 
gcloud compute copy-files --zone us-central1-f huangkuanlin@kuan-merge-genotype-bigmem:~/plink*tar plink_out
# transfer to storage
gsutil cp plink*tar gs://dinglab/isb-cgc/tcga/analysis_files