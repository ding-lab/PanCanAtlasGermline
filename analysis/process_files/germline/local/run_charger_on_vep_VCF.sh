#!/bin/bash

#setEnvironment=". /gscuser/ascott/python2_7_13_env"

inputDir="/gscmnt/gc3020/dinglab/medseq/Germline/projects/PanCanAtlasGermline/analysis/process_files/germline/local/AnnotatedVCFs"

#clinvar="/gscmnt/gc2706/dinglab/medseq/ClinVar/MacArthurLab/DataSnapshots/201701_ead5de/clinvar_alleles.tsv.gz"
clinvar="/gscmnt/gc2706/dinglab/medseq/ClinVar/MacArthurLab/clinvar/output/b37/single/clinvar_alleles.single.b37.tsv.gz"
mmGenes="/gscmnt/gc2737/ding/Analysis/Germline/CharGer/20160301_Rahman_KJ_KH_gene_table_CharGer.txt"
mmVariants="/gscmnt/gc2737/ding/Analysis/VariantLists/emptyRemoved_20160428_pathogenic_variants_HGVSg_VEP.vcf"
hotspot="/gscmnt/gc2737/ding/Analysis/Germline/CharGer/MC3.noHypers.mericUnspecified.d10.r20.v114.clusters"
rareThreshold="0.0005"

results="Charged_VEP/"
if [ ! -d ${results} ]; then
        mkdir ${results}
fi

#queue="bigmem"
queue="long"
queue="ding-lab"
group="/khuang"


for i in {1..22} X Y
do
	sample="$inputDir/anno.PCA.r1.TCGAbarcode.merge.exon.chr${i}.varOnly.vcf"
	vcf="PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf"
	output="${results}charged.${vcf}.tsv"
	command="charger -f ${sample} -o ${output} -O -D -g ${mmGenes} -z ${mmVariants} -H ${hotspot} -l --mac-clinvar-tsv ${clinvar} --rare-threshold ${rareThreshold} > ${results}charger.${vcf}.out"
	log="${results}charger.${vcf}.log"
	echo "bsub -R\"select[type==LINUX64 && mem>80000] rusage[mem=80000]\" -M 60000000 -g ${group} -q ${queue} -oo ${log} \"${command}\""
done
