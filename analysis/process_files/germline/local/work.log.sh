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
# normalize the VCF so we can map the CharGer var back to VCF easily
for i in {1..22} X Y; do
	echo "bsub -q short '~/bin/bcftools-1.5/bcftools norm -m - PCA.r1.TCGAbarcode.merge.exon.chr"$i".vcf.gz | bgzip -c > PCA.r1.TCGAbarcode.merge.exon.chr"${i}".norm.vcf.gz'"
done

# run VEP: ended up may only be able to do this on exon; careful of disk size limitation
bash run_VEP.v85.sh > run.cmd.sh
bash run.cmd.sh

# 6. Run CharGer (clinical classification)
# get the VEP annotated files to variant only first
#nohup bash shorten.vcf.sh &

bash run_charger_on_vep_VCF.sh > run.charger.cmd.sh
bash run.charger.cmd.sh

# 7. Map CharGer results back to VCF (samples and genotype)
# normalize the VCF so we can map the CharGer var back to VCF easily
for i in {1..22} X Y; do
	echo "bsub -q short '~/bin/bcftools-1.5/bcftools norm -m - anno.PCA.r1.TCGAbarcode.merge.exon.chr"$i".vcf | bgzip -c > anno.PCA.r1.TCGAbarcode.merge.exon.chr"${i}".norm.vcf.gz'"
done

bash post_CharGer.sh > postCharGerCmd.sh
bash postCharGerCmd.sh

# check CharGer liftover result
bash check.charger.list.sh 

# # clean up
# for i in {1..22} X Y; do 
# 	python remove_wrong_liftover.py ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf.samples.tsv \
# 	> ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf.samples.cleaned.tsv
# done
charged.PCA.r1.TCGAbarcode.merge.exon.chr7.vcf.samples.tsv
cat charged.PCA.r1.TCGAbarcode.merge.exon.chr*.vcf.samples.tsv | grep -v HUGO_Symbol > charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.tsv
cat charged.PCA.r1.TCGAbarcode.merge.exon.chrY.vcf.samples.tsv charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.tsv > tmp # chrY is essentially header
mv tmp charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.tsv

# 8. Post processing: non-pass variants don't have the ExAC frequency filter points from VEP annotation
# cut -f26-34 charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.tsv > charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.tsv.pos.vcf
# grep "##" ../AnnotatedVCFs/anno.PCA.r1.TCGAbarcode.merge.exon.chr1.varOnly.vcf > vcf.header.txt
# cat vcf.header.txt charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.tsv.pos.vcf > tmp
# mv -f tmp charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.tsv.pos.vcf

# may already have the AF but ignored previously 
python expand_csq.py new_run/ChargedSample/vcf.header.txt new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.tsv > new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.tsv
# filter on frequency and PM2 
python recalc_AF_PM2.py new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.tsv \
 > new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.tsv

awkt 'NR==1 || ($90 < 0.05 && $16 != "BA1")' new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.tsv > \
 new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.lowAF.tsv

cut -f 1-25,27- new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.lowAF.tsv > new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.lowAF.sele.tsv
wc -l new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL*
    #  49124 new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.lowAF.sele.tsv
    #  49124 new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.lowAF.tsv
    # 104668 new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.tsv
    # 122756 new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.tsv
    # 122756 new_run/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.tsv

### write script to clean each AF fields and use that

### run charger on exac
bsubl -oo extract_bed.log 'tabix -B -h ExAC.r1.sites.vep.biallelic.combine.vcf.gz /gscmnt/gc3020/dinglab/medseq/Germline/projects/PanCanAtlasGermline/TCGA_data/reference_files/all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_ForMusic | bgzip -c > ExAC.r1.sites.vep.biallelic.combine.exon.vcf.gz'

# by chromosome
bsubl 'tabix -p vcf ExAC.r1.sites.vep.biallelic.combine.exon.vcf.gz'
for i in {1..22} X Y
do
	echo "bsubl 'tabix -h ExAC.r1.sites.vep.biallelic.combine.exon.vcf.gz "$i" | bgzip -c > byChromosome/ExAC.r1.sites.vep.biallelic.combine.exon.chr${i}.vcf.gz'"
done