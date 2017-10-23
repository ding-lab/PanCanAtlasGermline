#!/bin/bash
# Kuan, Oct 2017
# 1. filter and merge
bash make.bsub.commands.sh
# calls filter_merge_germline_by_sample.sh on each chunk of the data

# 2. reheader to TCGA barcode
bsubl -oo update_vcfHeader_to_TCGA.log 'bash update_vcfHeader_to_TCGA.sh'
# calls replace_vcf_header_sample_with_source_TCGA.pl to find TCGA barcode corresponding to each input file and update in VCF header using tabix

# 3. Merge everything
ls -1 *merge.TCGAbarcode.vcf.gz | perl -pe 's/\n/ /g'
bsub -q bigmem -oo merge_all.log '~/bin/bcftools-1.5/bcftools merge --output-type z --output PCA.r1.TCGAbarcode.merge.vcf.gz PCA_aa.merge.TCGAbarcode.vcf.gz PCA_ab.merge.TCGAbarcode.vcf.gz PCA_ac.merge.TCGAbarcode.vcf.gz PCA_ad.merge.TCGAbarcode.vcf.gz PCA_ae.merge.TCGAbarcode.vcf.gz PCA_af.merge.TCGAbarcode.vcf.gz PCA_ag.merge.TCGAbarcode.vcf.gz PCA_ah.merge.TCGAbarcode.vcf.gz PCA_ai.merge.TCGAbarcode.vcf.gz PCA_aj.merge.TCGAbarcode.vcf.gz PCA_ak.merge.TCGAbarcode.vcf.gz PCA_al.merge.TCGAbarcode.vcf.gz PCA_am.merge.TCGAbarcode.vcf.gz PCA_an.merge.TCGAbarcode.vcf.gz PCA_ao.merge.TCGAbarcode.vcf.gz PCA_ap.merge.TCGAbarcode.vcf.gz PCA_aq.merge.TCGAbarcode.vcf.gz PCA_ar.merge.TCGAbarcode.vcf.gz PCA_as.merge.TCGAbarcode.vcf.gz PCA_at.merge.TCGAbarcode.vcf.gz PCA_au.merge.TCGAbarcode.vcf.gz PCA_av.merge.TCGAbarcode.vcf.gz PCA_aw.merge.TCGAbarcode.vcf.gz PCA_ax.merge.TCGAbarcode.vcf.gz PCA_ay.merge.TCGAbarcode.vcf.gz PCA_az.merge.TCGAbarcode.vcf.gz PCA_ba.merge.TCGAbarcode.vcf.gz PCA_bb.merge.TCGAbarcode.vcf.gz PCA_bc.merge.TCGAbarcode.vcf.gz PCA_bd.merge.TCGAbarcode.vcf.gz PCA_be.merge.TCGAbarcode.vcf.gz PCA_bf.merge.TCGAbarcode.vcf.gz PCA_bg.merge.TCGAbarcode.vcf.gz PCA_bh.merge.TCGAbarcode.vcf.gz PCA_bi.merge.TCGAbarcode.vcf.gz PCA_bj.merge.TCGAbarcode.vcf.gz'
bsubl -oo tabix_all.log 'tabix -p vcf PCA.r1.TCGAbarcode.merge.vcf.gz'

# # alternatively split and merge first this may go faster
# for file in PCA_*merge.TCGAbarcode.vcf.gz; do 
# 	echo "bsubl 'tabix -B -h "$file" ../../../../TCGA_data/reference_files/all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_ForMusic | bgzip -c > "${file%.vcf.gz}".exon.vcf.gz'"
# done
# ls -1 *.merge.TCGAbarcode.exon.vcf.gz | perl -pe 's/\n/ /g'
# bsub -q bigmem -oo merge_all_exon.log '~/bin/bcftools-1.5/bcftools merge --output-type z --output PCA.r1.TCGAbarcode.merge.exon.vcf.gz PCA_aa.merge.TCGAbarcode.exon.vcf.gz PCA_ab.merge.TCGAbarcode.exon.vcf.gz PCA_ac.merge.TCGAbarcode.exon.vcf.gz PCA_ad.merge.TCGAbarcode.exon.vcf.gz PCA_ae.merge.TCGAbarcode.exon.vcf.gz PCA_af.merge.TCGAbarcode.exon.vcf.gz PCA_ag.merge.TCGAbarcode.exon.vcf.gz PCA_ah.merge.TCGAbarcode.exon.vcf.gz PCA_ai.merge.TCGAbarcode.exon.vcf.gz PCA_aj.merge.TCGAbarcode.exon.vcf.gz PCA_ak.merge.TCGAbarcode.exon.vcf.gz PCA_al.merge.TCGAbarcode.exon.vcf.gz PCA_am.merge.TCGAbarcode.exon.vcf.gz PCA_an.merge.TCGAbarcode.exon.vcf.gz PCA_ao.merge.TCGAbarcode.exon.vcf.gz PCA_ap.merge.TCGAbarcode.exon.vcf.gz PCA_aq.merge.TCGAbarcode.exon.vcf.gz PCA_ar.merge.TCGAbarcode.exon.vcf.gz PCA_as.merge.TCGAbarcode.exon.vcf.gz PCA_at.merge.TCGAbarcode.exon.vcf.gz PCA_au.merge.TCGAbarcode.exon.vcf.gz PCA_av.merge.TCGAbarcode.exon.vcf.gz PCA_aw.merge.TCGAbarcode.exon.vcf.gz PCA_ax.merge.TCGAbarcode.exon.vcf.gz PCA_ay.merge.TCGAbarcode.exon.vcf.gz PCA_az.merge.TCGAbarcode.exon.vcf.gz PCA_ba.merge.TCGAbarcode.exon.vcf.gz PCA_bb.merge.TCGAbarcode.exon.vcf.gz PCA_bc.merge.TCGAbarcode.exon.vcf.gz PCA_bd.merge.TCGAbarcode.exon.vcf.gz PCA_be.merge.TCGAbarcode.exon.vcf.gz PCA_bf.merge.TCGAbarcode.exon.vcf.gz PCA_bg.merge.TCGAbarcode.exon.vcf.gz PCA_bh.merge.TCGAbarcode.exon.vcf.gz PCA_bi.merge.TCGAbarcode.exon.vcf.gz PCA_bj.merge.TCGAbarcode.exon.vcf.gz'
# bsubl -oo tabix_all_exon.log 'tabix -p vcf PCA.r1.TCGAbarcode.merge.exon.vcf.gz'

# 4. Extract exon bed
#tabix -R region.txt my.vcf.gz
bsubl -oo extract_bed.log 'tabix -B -h PCA.r1.TCGAbarcode.merge.vcf.gz ../../../../TCGA_data/reference_files/all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_ForMusic | bgzip -c > PCA.r1.TCGAbarcode.merge.exon.vcf.gz'
# zcat PCA.r1.TCGAbarcode.merge.vcf.gz | head -n 10000 | awk '/^#/' > PCA.r1.TCGAbarcode.head.txt
# tabix -r PCA.r1.TCGAbarcode.head.txt PCA.r1.TCGAbarcode.merge.exon.nohead.vcf.gz > PCA.r1.TCGAbarcode.merge.exon.vcf.gz
bsubl 'tabix -p vcf PCA.r1.TCGAbarcode.merge.exon.vcf.gz'

# 5. VEP annotation [for both exon and everything]
# split the files by chromosome first
for i in {1..22} X Y
do
	echo "bsubl 'tabix -h PCA.r1.TCGAbarcode.merge.vcf.gz "$i" | bgzip -c > allVCF/PCA.r1.TCGAbarcode.merge.chr${i}.vcf.gz'"
	echo "bsubl 'tabix -h PCA.r1.TCGAbarcode.merge.exon.vcf.gz "$i" | bgzip -c > VCF/PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf.gz'"
done

# run VEP: ended up may only be able to do this on exon; careful of disk size limitation
bash run_VEP.v85.sh > run.cmd.sh
bash run.cmd.sh

# 6. Run CharGer (clinical classification)
# get the VEP annotated files to variant only first
nohup bash shorten.vcf.sh &

bash run_charger_on_vep_VCF.sh > run.charger.cmd.sh
bash run.charger.cmd.sh


