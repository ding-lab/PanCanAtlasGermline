# do cancer type by cancer type to avoid file size
$ cat cancer_type.txt 
ACC
BLCA
BRCA
CESC
CHOL
COAD
DLBC
ESCA
GBM
HNSC
KICH
KIRC
KIRP
LAML # there are no LAML samples so this was throwing errors earlier
LGG
LIHC
LUAD
LUSC
MESO
OV
PAAD
PCPG
PRAD
READ
SARC
SKCM
STAD
TGCT
THCA
THYM
UCEC
UCS
UVM

$ cat merge_germline_by_cancer.sh 
#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "Processing cancer type: $line"
    echo "Start time"
    date
    cancer="${line%\\n}"
        mkdir $cancer
        samples=$(gsutil ls gs://dinglab/isb-cgc/tcga/germline/production/${cancer}/)
    for sample in $samples; do 
        sampleNamePre=${sample##*$cancer/}
        sampleName=${sampleNamePre%/}
        echo "Copying vcf for "$sampleName
        #gsutil cp ${sample}combine/*gz ${sampleName}.prefilter.snp_indel.vcf.gz
        #gsutil cp ${sample}combine/*gz.tbi ${sampleName}.prefilter.snp_indel.vcf.gz.tbi
        gsutil cp ${sample}combine/prefilter.snp_indel.annotated.ExAC_AF.0.01.AD.3.ROI.vcf.gz ${cancer}/${sampleName}.prefilter.snp_indel.annotated.ExAC_AF.0.01.AD.3.ROI.vcf.gz
        gsutil cp ${sample}combine/prefilter.snp_indel.annotated.ExAC_AF.0.01.AD.3.ROI.vcf.gz.tbi ${cancer}/${sampleName}.prefilter.snp_indel.annotated.ExAC_AF.0.01.AD.3.ROI.vcf.gz.tbi
    done
    
    #merge
    vcf-merge $(ls -1 ${cancer}/*.vcf.gz | perl -pe 's/\n/ /g') > ${cancer}.merge.vcf
    bgzip -c ${cancer}.merge.vcf > ${cancer}.merge.vcf.gz
    #index
    tabix -p vcf ${cancer}.merge.vcf.gz
    #upload
    gsutil cp ${cancer}.merge.vcf.gz* gs://dinglab/isb-cgc/tcga/germline/production/merge
    # delete files
    rm -rf ${cancer}.merge.vcf
    rm -rf ${cancer}/*ROI.vcf.gz
    rm -rf ${cancer}/*ROI.vcf.gz.tbi
    echo "End time"
    date
done < "$1"

# set unlimit file number higher so it works for breast cancer

ulimit -n 2500

nohup bash merge_germline_by_cancer.sh cancer_type.txt > merge_germline_by_cancer.log &
nohup bash merge_germline_by_cancer2.sh cancer_type2.txt > merge_germline_by_cancer2.log &
nohup bash merge_germline_by_cancer2.sh cancer_type3.txt > merge_germline_by_cancer3.log &
nohup bash merge_germline_by_cancer2.sh cancer_type4.txt > merge_germline_by_cancer4.log &
nohup bash merge_germline_by_cancer2.sh CESC.txt > merge_germline_by_cancerCESC.log &

# there is no LAML? if so delete
rm -rf LAML*

### merge verything
# copy file to big VM
echo "Start time"
date
gsutil cp gs://dinglab/isb-cgc/tcga/germline/production/merge/* .
#merge
echo "Merging"
#vcf-merge $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') > PCA.merge.vcf
bcftools merge --output-type z --output PCA.merge.vcf.gz $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') 
#echo "Zipping"
#bgzip -c PCA.merge.vcf > PCA.merge.vcf.gz
#index
echo "Indexing"
tabix -p vcf PCA.merge.vcf.gz
#upload
gsutil cp PCA.merge.vcf.gz* gs://dinglab/isb-cgc/tcga/germline/production/merge
echo "End time"
date

# # note: bcftools turned out much more efficient and can directly generate zipped file

# # note: other tools 
# # copy local files to VM
# gcloud compute copy-files --zone us-central1-f /Users/khuang/Downloads/GenomeAnalysisTK-3.7.tar.bz2 huangkuanlin@kuan-merge-germline-bigmem:~/
# gcloud compute copy-files --zone us-central1-f /Users/khuang/Downloads/picard-2.9.0.zip  huangkuanlin@kuan-merge-germline-bigmem:~/

# # set up
# wget ftp://genome.wustl.edu/pub/reference/GRCh37-lite/GRCh37-lite.fa.gz
# gunzip GRCh37-lite.fa.gz 
# samtools faidx GRCh37-lite.fa