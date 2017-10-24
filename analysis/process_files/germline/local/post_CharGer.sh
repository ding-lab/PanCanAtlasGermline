#!/bin/bash

input="Charged_VEP/"
results="ChargedSample/"
variantBuffer=2000
queue="long"
group="/khuang"
forks=4
mem=20000000

inputDir="/gscmnt/gc3020/dinglab/medseq/Germline/projects/PanCanAtlasGermline/analysis/process_files/germline/local/Charged_VEP/"

if [ ! -d ${results} ]; then
	mkdir ${results}
fi

for i in {1..22} X Y
do
	tsv=${input}charged.PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf.tsv
	vcf=AnnotatedVCFs/anno.PCA.r1.TCGAbarcode.merge.exon.chr${i}.norm.vcf.gz
	out=${results}charged.PCA.r1.TCGAbarcode.merge.exon.chr${i}.vcf.samples.tsv
	runCMD="python combine_CharGer2VCF.py ${tsv} ${vcf} > ${out}"
	log="${out}.log"
	echo "bsub -g ${group} -q ${queue} -n ${forks} -M ${mem} -oo ${log} \"${runCMD}\""
done
