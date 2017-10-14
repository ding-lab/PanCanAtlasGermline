#!/bin/bash

inputDir="/gscmnt/gc3020/dinglab/medseq/Germline/projects/pan8000_germline_clinical/variant_files/201604_PCGP_variants/"
mmGenes="/gscmnt/gc2737/ding/Analysis/Germline/CharGer/20160301_Rahman_KJ_KH_gene_table_CharGer.txt"
mmVariants="/gscmnt/gc2737/ding/Analysis/VariantLists/emptyRemoved_20160428_pathogenic_variants_HGVSg_VEP.vcf"
hotspot="/gscmnt/gc2737/ding/Analysis/Germline/CharGer/MC3.noHypers.mericUnspecified.d10.r20.v114.clusters"
clinvar="/gscmnt/gc2706/dinglab/medseq/ClinVar/MacArthurLab/clinvar/output/b37/single/clinvar_alleles.single.b37.tsv.gz"
rareThreshold="0.0005"
results="Charged/"
if [ ! -d ${results} ]; then
        mkdir ${results}
fi

queue="ding-lab"
#queue="long"
group="/khuang"

	sample="/gscmnt/gc3020/dinglab/medseq/Germline/projects/pan8000_germline_clinical/variant_files/201604_PCGP_variants/2015_stJude_germline_nejm_S4_AD_varOnly.vep.vcf"
	output="${results}charged.2015_stJude_germline_nejm_S4_AD_varOnly.vep.tsv"
	command="python CharGer/bin/charger -f ${sample} -o ${output} -O -D -g ${mmGenes} -z ${mmVariants} -H ${hotspot} -l --mac-clinvar-tsv ${clinvar} > ${results}charger.PCGP_AD.out"
	log="${results}charger.PCGP_AD.log"
	echo "bsub -g ${group} -q ${queue} -oo ${log} \"${command}\""

        sample="/gscmnt/gc3020/dinglab/medseq/Germline/projects/pan8000_germline_clinical/variant_files/201604_PCGP_variants/2015_stJude_germline_nejm_S4_AR_varOnly.vep.vcf"
        output="${results}charged.2015_stJude_germline_nejm_S4_AR_varOnly.vep.tsv"
        command="python CharGer/bin/charger -f ${sample} -o ${output} -O -D -g ${mmGenes} -z ${mmVariants} -H ${hotspot} -l --mac-clinvar-tsv ${clinvar} > ${results}charger.PCGP_AR.out"
        log="${results}charger.PCGP_AR.log"
        echo "bsub -g ${group} -q ${queue} -oo ${log} \"${command}\""
